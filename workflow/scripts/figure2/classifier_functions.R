# Shared functions of the label-trained 20-gene classifier (Figure 2D-E). They follow the code of the
# leave-one-study-out computation whose outputs are deposited (run on a compute cluster).
# Sourced by gene_panel_effects.R, primary_label_classifier.R and hio_grafts.R.
#
#   Training labels      in studies whose own object names EC and other terminal EEC cells
#                        (Egozi 2021, Elmentaite 2021, Hickey 2023) only those primary study labels
#                        are used; elsewhere the four-marker classifier's results (high-confidence EC
#                        positive, insufficient evidence negative).
#   Gene ranking         limma moderated t for EC versus other, on donor/sample-by-class mean
#                        log-normalized expression, paired within donor/sample, each donor/sample
#                        weighted 1 / (donors/samples in its dataset), so each of the seven gene-ranking
#                        datasets counts equally (Elmentaite 2020 and 2021 are two datasets of one study). The four classifier markers,
#                        MT- and RP[SL] genes and genes detected (>= 5 % of cells) in fewer than three
#                        datasets are not eligible. No redundancy filter: panel k = top k genes.
#   Classifier           weights = positive median per-dataset effect (normalized to sum one);
#                        score = weighted sum of genes standardized with equal weight per dataset;
#                        threshold = equal-posterior boundary of a two-component equal-variance
#                        Gaussian mixture fitted on training cells meeting detection, equal weight per
#                        dataset; call = >= 25 % of panel genes detected, score >= threshold, and a
#                        within-dataset 25-nearest-neighbor local-consistency gate (q <= 0.05).
suppressPackageStartupMessages(library(Matrix))

MIN_CELLS_PER_UNIT_CLASS <- 5L
TRAINING_DETECTION_FLOOR <- 0.05
MIN_TRAINING_FAMILIES <- 3L
SELECTION_EXCLUSIONS <- c("FEV", "LMX1A", "SLC18A1", "TPH1")
GMM_SPEC <- list(pairs = list(c(0.1, 0.9), c(0.25, 0.75), c(0.4, 0.8)), max_iterations = 1000L, tolerance = 1e-8,
                 variance_floor = 0.0025, minimum_component_weight = 0.02, minimum_component_sd = 0.05,
                 minimum_mean_separation_pooled_sd = 0.5)
GATE <- list(k = 25L, q_threshold = 0.05, highly_variable_genes = 250L, principal_components = 20L,
             detection_minimum = 0.01, detection_maximum = 0.99, minimum_features = 20L)

require_true <- function(value, message) if (!isTRUE(value)) stop(message, call. = FALSE)

# ---- training labels ------------------------------------------------------------------------
training_labels <- function(metadata) {
  presentation <- metadata$presentation_token; author <- metadata$author_group_token
  stopifnot(all(presentation %in% c("HIGH_CONFIDENCE_EC", "LOW_CONFIDENCE_EC", "INSUFFICIENT_EC_EVIDENCE")))
  author_ec <- author %in% c("EXPLICIT_AUTHOR_EC", "EXPLICIT_EC", "EXACT_AUTHOR_EC", "AUTHOR_EC")
  author_other <- author %in% c("EXPLICIT_TERMINAL_NON_EC_EEC", "TERMINAL_NON_EC_EEC", "EXACT_AUTHOR_TERMINAL_NON_EC_EEC", "AUTHOR_TERMINAL_NON_EC_EEC")
  label <- rep("ABSTAIN", nrow(metadata)); reason <- rep("AUTHOR_SILENT", nrow(metadata))
  label[presentation == "HIGH_CONFIDENCE_EC" & !author_other] <- "POSITIVE"   # four-marker classifier labels
  label[presentation == "INSUFFICIENT_EC_EVIDENCE"] <- "NEGATIVE"
  own <- metadata$author_label_evidence_layer %in% "ORIGINAL_AUTHOR_OBJECT"
  positive <- own & author_ec; negative <- own & author_other
  author_study <- metadata$family_id %in% unique(metadata$family_id[positive | negative])
  label[author_study] <- "ABSTAIN"; reason[author_study] <- "AUTHOR_STUDY_CELL_WITHOUT_TERMINAL_AUTHOR_LABEL"
  label[positive] <- "POSITIVE"; reason[positive] <- "ORIGINAL_AUTHOR_EC"
  label[negative] <- "NEGATIVE"; reason[negative] <- "ORIGINAL_AUTHOR_TERMINAL_NON_EC_EEC"
  data.frame(project_cell_id = metadata$project_cell_id, label = label, reason = reason, stringsAsFactors = FALSE)
}

# ---- gene eligibility, unit means and the limma ranking ----------------------------------------
training_feature_mask <- function(log_expression, genes, metadata, train) {
  pass_count <- integer(length(genes))
  for (family in sort(unique(metadata$family_id[train]))) {
    cells <- train & metadata$family_id == family
    pass_count <- pass_count + (Matrix::rowMeans(log_expression[, cells, drop = FALSE] > 0) >= TRAINING_DETECTION_FLOOR)
  }
  pass_count >= MIN_TRAINING_FAMILIES & !(genes %in% SELECTION_EXCLUSIONS) &
    !grepl("^MT-", genes, ignore.case = TRUE) & !grepl("^RP[SL][0-9]", genes, ignore.case = TRUE)
}

unit_class_means <- function(log_expression, metadata, labels, cells, genes, minimum_cells) {
  labels <- labels[match(metadata$project_cell_id, labels$project_cell_id), , drop = FALSE]
  cells <- as.logical(cells) & labels$label %in% c("POSITIVE", "NEGATIVE")
  gene_index <- match(genes, rownames(log_expression)); stopifnot(!anyNA(gene_index))
  grouping <- paste(metadata$family_id, metadata$biological_unit_id, labels$label, sep = "::")
  groups <- split(which(cells), grouping[cells])
  groups <- groups[vapply(groups, length, integer(1L)) >= as.integer(minimum_cells)]
  values <- vapply(groups, function(index) as.numeric(Matrix::rowMeans(log_expression[gene_index, index, drop = FALSE])),
                   numeric(length(genes)))
  if (is.null(dim(values))) values <- matrix(values, nrow = length(genes))
  rownames(values) <- genes
  pieces <- strsplit(names(groups), "::", fixed = TRUE)
  list(values = values, metadata = data.frame(
    family_id = vapply(pieces, function(x) x[[1L]], ""),
    biological_unit_id = vapply(pieces, function(x) paste(x[2:(length(x) - 1L)], collapse = "::"), ""),
    label = vapply(pieces, function(x) x[[length(x)]], ""), stringsAsFactors = FALSE))
}

limma_ec_effect <- function(units) {
  m <- units$metadata
  unit <- factor(paste(m$family_id, m$biological_unit_id, sep = "::"))
  ec <- as.integer(m$label == "POSITIVE")
  units_per_family <- tapply(as.character(unit), m$family_id, function(v) length(unique(v)))
  w <- as.numeric(1 / units_per_family[m$family_id])
  design <- stats::model.matrix(~ 0 + unit + ec)
  fit <- limma::eBayes(limma::lmFit(units$values, design, weights = matrix(w, nrow(units$values), length(w), byrow = TRUE)), trend = TRUE)
  tt <- limma::topTable(fit, coef = "ec", number = Inf, sort.by = "none")
  data.frame(gene_symbol = rownames(units$values), logFC = tt$logFC, t = tt$t, p = tt$P.Value, fdr = tt$adj.P.Val, stringsAsFactors = FALSE)
}

# All-training ranking of eligible genes, best first (columns gene_symbol, limma_t, limma_fdr, ...).
rank_training_genes <- function(log_expression, genes, metadata, labels, training_cells) {
  eligible_genes <- genes[training_feature_mask(log_expression, genes, metadata, training_cells)]
  units <- unit_class_means(log_expression, metadata, labels, training_cells, eligible_genes, MIN_CELLS_PER_UNIT_CLASS)
  fit <- limma_ec_effect(units)
  r <- data.frame(gene_symbol = eligible_genes, limma_logfc = fit$logFC[match(eligible_genes, fit$gene_symbol)],
                  limma_t = fit$t[match(eligible_genes, fit$gene_symbol)], limma_p = fit$p[match(eligible_genes, fit$gene_symbol)],
                  limma_fdr = fit$fdr[match(eligible_genes, fit$gene_symbol)], stringsAsFactors = FALSE)
  stopifnot(!anyNA(r$limma_t))
  r <- r[order(-r$limma_t, r$gene_symbol, method = "radix"), , drop = FALSE]; rownames(r) <- NULL
  r
}

# ---- the classifier --------------------------------------------------------------------------
family_balanced_weights <- function(metadata, cells) {
  family <- metadata$family_id[as.logical(cells)]
  w <- 1 / as.numeric(table(family)[family]); w / sum(w)
}
fit_standardizer <- function(log_expression, metadata, training_cells, genes) {
  index <- match(genes, rownames(log_expression)); stopifnot(!anyNA(index))
  w <- family_balanced_weights(metadata, training_cells); w <- w / sum(w)   # normalized again, as in the deposited computation
  e <- log_expression[index, training_cells, drop = FALSE]
  mean <- as.numeric(e %*% w); second <- as.numeric((e^2) %*% w)
  data.frame(gene_symbol = genes, training_mean = mean, training_sd = sqrt(pmax(second - mean^2, 0)), stringsAsFactors = FALSE)
}
apply_standardizer <- function(log_expression, standardizer) {
  index <- match(standardizer$gene_symbol, rownames(log_expression)); stopifnot(!anyNA(index))
  dense <- as.matrix(log_expression[index, , drop = FALSE])
  sweep(sweep(dense, 1L, standardizer$training_mean, "-"), 1L, standardizer$training_sd, "/")
}
detection_requirement <- function(panel_size, fraction = 0.25) min(panel_size, max(2L, as.integer(ceiling(fraction * panel_size))))

weighted_quantile <- function(values, weights, probability) {
  o <- order(values, seq_along(values), method = "radix")
  v <- values[o]; w <- weights[o]
  as.numeric(stats::approx(cumsum(w) - 0.5 * w, v, xout = probability, rule = 2, ties = "ordered")$y)
}
mixture_log_likelihood <- function(values, weights, mixture, means, variance) {
  ld <- cbind(log(mixture[[1L]]) - 0.5 * log(2 * pi * variance) - 0.5 * (values - means[[1L]])^2 / variance,
              log(mixture[[2L]]) - 0.5 * log(2 * pi * variance) - 0.5 * (values - means[[2L]])^2 / variance)
  row_max <- apply(ld, 1L, max)
  sum(weights * (row_max + log(rowSums(exp(ld - row_max)))))
}
fit_equal_variance_gmm <- function(values, weights, spec = GMM_SPEC) {
  keep <- is.finite(values) & is.finite(weights) & weights > 0
  values <- values[keep]; weights <- weights[keep] / sum(weights[keep])
  overall_variance <- sum(weights * (values - sum(values * weights))^2)
  if (overall_variance <= spec$variance_floor) return(list(status = "UNAVAILABLE", threshold = NA_real_))
  best <- NULL
  for (init_index in seq_along(spec$pairs)) {
    pair <- spec$pairs[[init_index]]
    means <- c(weighted_quantile(values, weights, pair[1]), weighted_quantile(values, weights, pair[2]))
    if (means[[1L]] >= means[[2L]]) next
    mixture <- c(0.5, 0.5); variance <- max(overall_variance, spec$variance_floor)
    previous <- mixture_log_likelihood(values, weights, mixture, means, variance); converged <- FALSE
    for (iteration in seq_len(spec$max_iterations)) {
      ld <- cbind(log(mixture[[1L]]) - 0.5 * log(2 * pi * variance) - 0.5 * (values - means[[1L]])^2 / variance,
                  log(mixture[[2L]]) - 0.5 * log(2 * pi * variance) - 0.5 * (values - means[[2L]])^2 / variance)
      ex <- exp(ld - apply(ld, 1L, max))
      wr <- (ex / rowSums(ex)) * weights
      mass <- colSums(wr); if (any(mass <= 0)) break
      next_mixture <- mass / sum(mass)
      next_means <- colSums(wr * values) / mass
      next_variance <- max(sum(wr * sweep(matrix(values, nrow = length(values), ncol = 2L), 2L, next_means, "-")^2), spec$variance_floor)
      ll <- mixture_log_likelihood(values, weights, next_mixture, next_means, next_variance)
      mixture <- next_mixture; means <- next_means; variance <- next_variance
      if (abs(ll - previous) <= spec$tolerance) { converged <- TRUE; break }
      previous <- ll
    }
    if (!converged) next
    o <- order(means, method = "radix"); means <- means[o]; mixture <- mixture[o]
    cand <- list(ll = mixture_log_likelihood(values, weights, mixture, means, variance), init = init_index,
                 mixture = mixture, means = means, variance = variance)
    if (is.null(best) || cand$ll > best$ll) best <- cand
  }
  unavailable <- list(status = "UNAVAILABLE", threshold = NA_real_)
  if (is.null(best)) return(unavailable)
  sd <- sqrt(best$variance)
  if (any(best$mixture < spec$minimum_component_weight) || sd < spec$minimum_component_sd ||
      (best$means[[2L]] - best$means[[1L]]) / sd < spec$minimum_mean_separation_pooled_sd) return(unavailable)
  thr <- 0.5 * sum(best$means) + best$variance * log(best$mixture[[1L]] / best$mixture[[2L]]) / diff(best$means)
  if (!is.finite(thr) || thr < best$means[[1L]] || thr > best$means[[2L]]) return(unavailable)
  list(status = "PASS", threshold = thr)
}

fit_panel_model <- function(log_expression, metadata, labels, training_cells, panel_genes, detection_fraction = 0.25) {
  standardizer <- fit_standardizer(log_expression, metadata, training_cells, panel_genes)
  units <- unit_class_means(log_expression, metadata, labels, training_cells, panel_genes, 1L)
  families <- sort(unique(units$metadata$family_id))
  family_effects <- vapply(panel_genes, function(gene) {
    stats::median(vapply(families, function(family) {
      block <- units$metadata$family_id == family
      p <- block & units$metadata$label == "POSITIVE"; n <- block & units$metadata$label == "NEGATIVE"
      if (!any(p) || !any(n)) return(NA_real_)
      mean(units$values[gene, p]) - mean(units$values[gene, n])
    }, numeric(1L)), na.rm = TRUE)
  }, numeric(1L))
  raw <- pmax(family_effects, 0); raw[!is.finite(raw)] <- 0
  require_true(any(raw > 0), "panel has no positive training-family effect weights")
  weights <- raw / sum(raw)
  scores <- as.numeric(Matrix::crossprod(weights, apply_standardizer(log_expression, standardizer)))
  detected <- as.numeric(Matrix::colSums(log_expression[match(panel_genes, rownames(log_expression)), , drop = FALSE] > 0))
  requirement <- detection_requirement(length(panel_genes), detection_fraction)
  gmm_cells <- training_cells & detected >= requirement & is.finite(scores)
  gmm <- fit_equal_variance_gmm(scores[gmm_cells], family_balanced_weights(metadata, gmm_cells))
  require_true(identical(gmm$status, "PASS"), "panel mixture unavailable")
  list(panel_genes = panel_genes, gene_weights = weights, standardizer = standardizer,
       detection_requirement = requirement, threshold = gmm$threshold,
       training_families = sort(unique(metadata$family_id[training_cells])))
}
score_panel_model <- function(log_expression, model) {
  index <- match(model$panel_genes, rownames(log_expression)); stopifnot(!anyNA(index))
  list(score = as.numeric(Matrix::crossprod(model$gene_weights, apply_standardizer(log_expression, model$standardizer))),
       detected_genes = as.numeric(Matrix::colSums(log_expression[index, , drop = FALSE] > 0)))
}

# ---- within-dataset neighborhood gate -------------------------------------------------------
family_graph <- function(log_expression, excluded_genes, gate = GATE) {
  block <- log_expression; n <- ncol(block); genes <- rownames(block)
  detection <- as.numeric(Matrix::rowMeans(block > 0))
  sums <- as.numeric(Matrix::rowSums(block)); squares <- as.numeric(Matrix::rowSums(block^2))
  variance <- (squares - sums^2 / n) / (n - 1)
  eligible <- detection >= gate$detection_minimum & detection <= gate$detection_maximum &
    is.finite(variance) & variance > 0 & !(genes %in% excluded_genes)
  stopifnot(sum(eligible) >= gate$minimum_features)
  ordered <- order(-variance, genes, method = "radix"); ordered <- ordered[eligible[ordered]]
  features <- genes[ordered[seq_len(min(gate$highly_variable_genes, length(ordered)))]]
  dense <- t(as.matrix(block[match(features, genes), , drop = FALSE]))
  center <- colMeans(dense)
  psd <- sqrt(colMeans(sweep(dense, 2L, center)^2)); psd[psd <= 0] <- 1
  dense <- sweep(sweep(dense, 2L, center), 2L, psd, "/"); dense[!is.finite(dense)] <- 0
  rank <- min(gate$principal_components, ncol(dense), n - 1L)
  d <- svd(dense, nu = rank, nv = 0L)
  distance <- as.matrix(stats::dist(d$u %*% diag(d$d[seq_len(rank)], rank, rank)))^2
  diag(distance) <- Inf
  k <- min(gate$k, n - 1L)
  list(cells = colnames(block), neighbours = t(apply(distance, 1L, function(row) order(row)[seq_len(k)])), k = k)
}
gated_call <- function(score, detected_genes, requirement, threshold, graph, gate = GATE) {
  stage_one <- detected_genes >= requirement & is.finite(score) & score >= threshold
  nb <- graph$neighbours; total <- length(stage_one); positives <- sum(stage_one)
  hits <- as.integer(rowSums(matrix(stage_one[nb], nrow = nrow(nb))))
  q <- rep(NA_real_, total)
  if (positives > 0L) {
    p <- stats::phyper(hits[stage_one] - 1L, positives - 1L, (total - 1L) - (positives - 1L), ncol(nb), lower.tail = FALSE)
    q[stage_one] <- stats::p.adjust(p, method = "BH")
  }
  list(stage_one = stage_one, call = stage_one & !is.na(q) & q <= gate$q_threshold)
}

binary_auc <- function(response, score) {
  keep <- is.finite(score) & response %in% c(0L, 1L); response <- response[keep]; score <- score[keep]
  pos <- sum(response == 1L); neg <- sum(response == 0L)
  if (pos == 0L || neg == 0L) return(NA_real_)
  r <- rank(score, ties.method = "average")
  (sum(r[response == 1L]) - pos * (pos + 1) / 2) / (pos * neg)
}
