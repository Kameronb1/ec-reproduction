# Functions of the label-trained gene-panel classifier (the 20-gene classifier). Source
# figure2_postvalidation_core.R and figure2_postvalidation_analysis.R before this file.
#
#   - genes are ranked by a study-weighted limma model of donor/sample-by-class means
#     (cp_rank_training_features);
#   - a panel is fitted on the training studies only: per-gene standardization, gene weights and a
#     two-component mixture threshold (cp_fit_panel_model);
#   - a cell is classified as EC in two stages: the detection requirement and a score at or above the
#     threshold, then the four-marker classifier's local-enrichment check inside the cell's own
#     study, from that study's unlabeled cells (cp_two_stage_call);
#   - matched-random null panels are drawn to the real panel's expression, detection and variance
#     bins (cp_draw_matched_random_panels).

cp_read_json <- function(path) {
  require_true(file.exists(path), sprintf("missing JSON: %s", path))
  jsonlite::fromJSON(path, simplifyVector = FALSE)
}

cp_detection_requirement <- function(panel_size, fraction = 0.25) {
  panel_size <- as.integer(panel_size)
  require_true(length(panel_size) == 1L && is.finite(panel_size) && panel_size >= 1L, "panel size must be a positive integer")
  min(panel_size, max(2L, as.integer(ceiling(fraction * panel_size))))
}

cp_align_labels <- function(metadata, labels) {
  index <- match(metadata$project_cell_id, labels$project_cell_id)
  require_true(!anyNA(index), "labels do not cover every metadata cell")
  labels[index, , drop = FALSE]
}

cp_family_balanced_weights <- function(metadata, cells) {
  cells <- as.logical(cells)
  require_true(length(cells) == nrow(metadata) && any(cells), "study weights require selected cells")
  family <- metadata$family_id[cells]
  size <- table(family)
  weights <- 1 / as.numeric(size[family])
  weights / sum(weights)
}

cp_weighted_row_moments <- function(expression, weights) {
  require_true(ncol(expression) == length(weights), "expression and weight lengths differ")
  require_true(all(is.finite(weights)) && all(weights >= 0) && sum(weights) > 0, "invalid expression weights")
  weights <- weights / sum(weights)
  mean <- as.numeric(expression %*% weights)
  second <- as.numeric((expression ^ 2) %*% weights)
  variance <- pmax(second - mean ^ 2, 0)
  list(mean = mean, sd = sqrt(variance))
}

cp_fit_standardizer <- function(log_expression, metadata, training_cells, genes) {
  # NULL rownames on a dgCMatrix read back by readRDS mean Matrix is loaded but not attached;
  # that is reported as such, not as missing genes.
  require_true(!is.null(rownames(log_expression)),
               paste0("expression matrix has NULL rownames: attach Matrix (library(Matrix)) ",
                      "before any panel fit; the genes are not missing."))
  index <- match(genes, rownames(log_expression))
  # The error names the missing genes, so a failure can be checked against the data.
  require_true(!anyNA(index), sprintf("panel genes are missing from expression: %s",
                                      paste(genes[is.na(index)], collapse = ", ")))
  weights <- cp_family_balanced_weights(metadata, training_cells)
  moments <- cp_weighted_row_moments(log_expression[index, training_cells, drop = FALSE], weights)
  require_true(all(is.finite(moments$mean)) && all(is.finite(moments$sd)), "standardization moments are non-finite")
  require_true(all(moments$sd > 0), "at least one panel gene has zero training variance")
  data.frame(gene_symbol = genes, training_mean = moments$mean, training_sd = moments$sd, stringsAsFactors = FALSE)
}

cp_apply_standardizer <- function(log_expression, standardizer) {
  index <- match(standardizer$gene_symbol, rownames(log_expression))
  require_true(!anyNA(index), "standardizer genes are missing from expression")
  dense <- as.matrix(log_expression[index, , drop = FALSE])
  sweep(sweep(dense, 1L, standardizer$training_mean, "-"), 1L, standardizer$training_sd, "/")
}

cp_unit_class_means <- function(log_expression, metadata, labels, cells, genes,
                                minimum_cells_per_unit_class = 5L) {
  labels <- cp_align_labels(metadata, labels)
  cells <- as.logical(cells) & labels$label %in% c("POSITIVE", "NEGATIVE")
  gene_index <- match(genes, rownames(log_expression))
  require_true(!anyNA(gene_index), "unit-class genes are missing")
  grouping <- paste(metadata$family_id, metadata$biological_unit_id, labels$label, sep = "::")
  groups <- split(which(cells), grouping[cells])
  keep <- vapply(groups, length, integer(1L)) >= as.integer(minimum_cells_per_unit_class)
  groups <- groups[keep]
  require_true(length(groups) > 0L, "no biological-unit-class group meets the cell floor")
  values <- vapply(groups, function(index) {
    as.numeric(Matrix::rowMeans(log_expression[gene_index, index, drop = FALSE]))
  }, numeric(length(genes)))
  if (is.null(dim(values))) values <- matrix(values, nrow = length(genes))
  rownames(values) <- genes
  pieces <- strsplit(names(groups), "::", fixed = TRUE)
  group_metadata <- data.frame(
    group_id = names(groups),
    family_id = vapply(pieces, function(x) x[[1L]], character(1L)),
    biological_unit_id = vapply(pieces, function(x) paste(x[2:(length(x) - 1L)], collapse = "::"), character(1L)),
    label = vapply(pieces, function(x) x[[length(x)]], character(1L)),
    cells = vapply(groups, length, integer(1L)),
    stringsAsFactors = FALSE
  )
  list(values = values, metadata = group_metadata)
}

cp_training_gene_summary <- function(log_expression, genes, metadata, labels, training_cells,
                                     excluded_genes, detection_floor, minimum_training_families,
                                     minimum_cells_per_unit_class) {
  eligible <- training_feature_mask(
    log_expression, genes, metadata, training_cells,
    detection_floor, minimum_training_families, excluded_genes
  )
  eligible_genes <- genes[eligible]
  require_true(length(eligible_genes) >= 2L, "fewer than two eligible training genes")
  units <- cp_unit_class_means(
    log_expression, metadata, labels, training_cells, eligible_genes,
    minimum_cells_per_unit_class
  )
  families <- sort(unique(units$metadata$family_id))
  effect_blocks <- lapply(families, function(family) {
    block <- units$metadata$family_id == family
    positive <- block & units$metadata$label == "POSITIVE"
    negative <- block & units$metadata$label == "NEGATIVE"
    if (!any(positive) || !any(negative)) return(NULL)
    data.frame(
      family_id = family,
      gene_symbol = eligible_genes,
      effect = rowMeans(units$values[, positive, drop = FALSE]) -
        rowMeans(units$values[, negative, drop = FALSE]),
      stringsAsFactors = FALSE
    )
  })
  effect_blocks <- effect_blocks[!vapply(effect_blocks, is.null, logical(1L))]
  require_true(length(effect_blocks) >= as.integer(minimum_training_families), "too few training studies have both classes at donor/sample level")
  effects <- do.call(rbind, effect_blocks)
  labels_aligned <- cp_align_labels(metadata, labels)
  detection <- do.call(rbind, lapply(sort(unique(metadata$family_id[training_cells])), function(family) {
    block <- training_cells & metadata$family_id == family & labels_aligned$label %in% c("POSITIVE", "NEGATIVE")
    if (!any(block)) return(NULL)
    data.frame(
      family_id = family,
      gene_symbol = eligible_genes,
      detection_fraction = as.numeric(Matrix::rowMeans(log_expression[eligible, block, drop = FALSE] > 0)),
      stringsAsFactors = FALSE
    )
  }))
  effect_by_gene <- split(effects, effects$gene_symbol)
  detection_by_gene <- split(detection, detection$gene_symbol)
  summary <- do.call(rbind, lapply(eligible_genes, function(gene) {
    effect_block <- effect_by_gene[[gene]]
    detection_block <- detection_by_gene[[gene]]
    data.frame(
      gene_symbol = gene,
      effect_families = length(unique(effect_block$family_id)),
      positive_effect_families = length(unique(effect_block$family_id[effect_block$effect > 0])),
      median_family_effect = stats::median(effect_block$effect, na.rm = TRUE),
      minimum_family_effect = min(effect_block$effect, na.rm = TRUE),
      median_detection_fraction = stats::median(detection_block$detection_fraction, na.rm = TRUE),
      stringsAsFactors = FALSE
    )
  }))
  # Genes are ranked by limma's moderated t for EC versus the other class,
  # fitted on the same donor/sample-by-class means, paired within each donor/sample, with each
  # donor/sample weighted 1 / (donors/samples in its study) so every study contributes equally.
  fit <- cp_limma_ec_effect(units)
  summary$limma_logfc <- fit$logFC[match(summary$gene_symbol, fit$gene_symbol)]
  summary$limma_t <- fit$t[match(summary$gene_symbol, fit$gene_symbol)]
  summary$limma_p <- fit$p[match(summary$gene_symbol, fit$gene_symbol)]
  summary$limma_fdr <- fit$fdr[match(summary$gene_symbol, fit$gene_symbol)]
  require_true(!anyNA(summary$limma_t), "limma left a gene without a moderated t")
  summary <- summary[order(-summary$limma_t, summary$gene_symbol, method = "radix"), , drop = FALSE]
  rownames(summary) <- NULL
  list(summary = summary, effects = effects, detection = detection, unit_class = units)
}

cp_limma_ec_effect <- function(units) {
  m <- units$metadata
  unit <- factor(paste(m$family_id, m$biological_unit_id, sep = "::"))
  ec <- as.integer(m$label == "POSITIVE")
  both <- tapply(ec, unit, function(v) length(unique(v)) == 2L)
  require_true(any(both), "no donor/sample carries both classes, so the paired EC effect is undefined")
  # Study weight: each donor/sample counts 1 / (donors/samples in its study), counted as distinct units.
  units_per_family <- tapply(as.character(unit), m$family_id, function(v) length(unique(v)))
  w <- as.numeric(1 / units_per_family[m$family_id])
  design <- stats::model.matrix(~ 0 + unit + ec)
  fit <- limma::eBayes(limma::lmFit(units$values, design, weights = matrix(w, nrow(units$values), length(w), byrow = TRUE)), trend = TRUE)
  tt <- limma::topTable(fit, coef = "ec", number = Inf, sort.by = "none")
  data.frame(gene_symbol = rownames(units$values), logFC = tt$logFC, t = tt$t, p = tt$P.Value,
             fdr = tt$adj.P.Val, stringsAsFactors = FALSE)
}

cp_reduce_redundancy <- function(ranked_summary, unit_class, maximum_size, correlation_ceiling = 0.85) {
  ordered <- ranked_summary$gene_symbol
  selected <- character()
  decisions <- list()
  for (gene in ordered) {
    reason <- "KEEP"
    maximum_absolute_correlation <- 0
    if (length(selected)) {
      candidate <- as.numeric(unit_class$values[gene, ])
      correlations <- vapply(selected, function(kept) {
        suppressWarnings(stats::cor(candidate, as.numeric(unit_class$values[kept, ]), method = "spearman", use = "pairwise.complete.obs"))
      }, numeric(1L))
      correlations[!is.finite(correlations)] <- 0
      maximum_absolute_correlation <- max(abs(correlations))
      if (maximum_absolute_correlation > correlation_ceiling) reason <- "DROP_REDUNDANT"
    }
    if (identical(reason, "KEEP")) selected <- c(selected, gene)
    decisions[[length(decisions) + 1L]] <- data.frame(
      gene_symbol = gene, decision = reason,
      maximum_absolute_correlation_to_previously_kept = maximum_absolute_correlation,
      selected_rank = if (gene %in% selected) match(gene, selected) else NA_integer_,
      stringsAsFactors = FALSE
    )
    if (length(selected) >= as.integer(maximum_size)) break
  }
  list(genes = selected, decisions = do.call(rbind, decisions))
}

cp_rank_training_features <- function(log_expression, genes, metadata, labels, training_cells,
                                      excluded_genes, parameters, maximum_size = 30L,
                                      candidate_pool = NULL, correlation_ceiling = 0.85) {
  summary <- cp_training_gene_summary(
    log_expression, genes, metadata, labels, training_cells, excluded_genes,
    as.numeric(parameters$training_detection_floor),
    as.integer(parameters$minimum_training_families),
    as.integer(parameters$minimum_cells_per_unit_class)
  )
  ranked <- summary$summary
  if (!is.null(candidate_pool)) {
    ranked <- ranked[ranked$gene_symbol %in% candidate_pool, , drop = FALSE]
    require_true(nrow(ranked) >= 2L, "candidate pool leaves fewer than two eligible training genes")
  }
  reduced_units <- summary$unit_class
  reduced_units$values <- reduced_units$values[ranked$gene_symbol, , drop = FALSE]
  reduced <- cp_reduce_redundancy(ranked, reduced_units, maximum_size, correlation_ceiling)
  list(
    genes = reduced$genes,
    ranking = ranked,
    redundancy = reduced$decisions,
    effects = summary$effects,
    detection = summary$detection,
    unit_class = summary$unit_class
  )
}

# Fits one panel on the training cells: per-gene standardization, gene weights (each gene's median
# EC-minus-other effect across the training studies, negative effects set to 0, normalized to sum
# to 1), the detection requirement and the two-component mixture threshold. null_panel_policy is
# FALSE for every real panel, so a panel without a positive weight or without a mixture fit stops
# with an error. It is TRUE only for matched-random null panels, where both outcomes are legitimate
# draws from the null and are kept: equal weights, and no threshold.
cp_fit_panel_model <- function(log_expression, genes, metadata, labels, training_cells,
                               panel_genes, gmm_spec, detection_fraction = 0.25,
                               null_panel_policy = FALSE) {
  standardizer <- cp_fit_standardizer(log_expression, metadata, training_cells, panel_genes)
  units <- cp_unit_class_means(log_expression, metadata, labels, training_cells, panel_genes, 1L)
  families <- sort(unique(units$metadata$family_id))
  family_effects <- vapply(panel_genes, function(gene) {
    effects <- vapply(families, function(family) {
      block <- units$metadata$family_id == family
      positive <- block & units$metadata$label == "POSITIVE"
      negative <- block & units$metadata$label == "NEGATIVE"
      if (!any(positive) || !any(negative)) return(NA_real_)
      mean(units$values[gene, positive]) - mean(units$values[gene, negative])
    }, numeric(1L))
    stats::median(effects, na.rm = TRUE)
  }, numeric(1L))
  raw_weights <- pmax(family_effects, 0)
  raw_weights[!is.finite(raw_weights)] <- 0
  weight_rule <- "POSITIVE_TRAINING_FAMILY_EFFECT_NORMALIZED_TO_SUM_ONE"
  if (!any(raw_weights > 0)) {
    require_true(isTRUE(null_panel_policy), "panel has no positive training-family effect weights")
    raw_weights <- rep(1, length(panel_genes))
    weight_rule <- "EQUAL_WEIGHTS_FALLBACK__NO_POSITIVE_EFFECT__WORST_CASE_NULL"
  }
  weights <- raw_weights / sum(raw_weights)
  standardized <- cp_apply_standardizer(log_expression, standardizer)
  scores <- as.numeric(Matrix::crossprod(weights, standardized))
  panel_index <- match(panel_genes, rownames(log_expression))
  detected <- as.numeric(Matrix::colSums(log_expression[panel_index, , drop = FALSE] > 0))
  requirement <- cp_detection_requirement(length(panel_genes), detection_fraction)
  gmm_cells <- training_cells & detected >= requirement & is.finite(scores)
  require_true(any(gmm_cells), "no training cells pass panel detection for GMM")
  gmm_weights <- cp_family_balanced_weights(metadata, gmm_cells)
  fit <- fit_equal_variance_gmm(scores[gmm_cells], gmm_weights, gmm_spec)
  if (!identical(fit$status, "PASS")) {
    require_true(isTRUE(null_panel_policy), sprintf("panel GMM unavailable: %s", fit$reason))
    fit$threshold <- NA_real_
  }
  call <- if (is.finite(fit$threshold)) {
    detected >= requirement & scores >= fit$threshold
  } else {
    rep(NA, length(scores))
  }
  list(
    model = list(
      panel_genes = panel_genes,
      gene_weights = weights,
      weight_rule = weight_rule,
      training_family_effects = family_effects,
      standardizer = standardizer,
      detection_requirement = requirement,
      detection_fraction = detection_fraction,
      gmm = fit,
      training_families = sort(unique(metadata$family_id[training_cells]))
    ),
    score = scores,
    detected_genes = detected,
    call = call
  )
}

# ---- Two-stage call -------------------------------------------------------------------
# Stage one: the detection requirement and score >= training threshold minus a relaxation (the
# scripts here pass 0). Stage two: the four-marker classifier's local-enrichment check, applied
# inside the scored study using only that study's unlabeled cells. A neighborhood graph is built
# from the study's 250 most variable genes (detected in 1-99 percent of its cells, excluding the
# genes the caller names: the four EC markers and the panel genes), standardized to population SD,
# reduced to 20 principal components, with k = 25 nearest neighbors by squared Euclidean distance.
# A stage-one candidate passes if its neighborhood holds more candidates than a hypergeometric draw
# from the study-wide candidate prevalence predicts, Benjamini-Hochberg corrected within the study
# and within candidates, at q <= 0.05. The constants are the four-marker classifier's and are not tuned.

cp_gate_spec <- function(contract) {
  gate <- contract$local_enrichment_check
  list(
    k = as.integer(gate$k_nearest_neighbors),
    q_threshold = as.numeric(gate$q_threshold),
    highly_variable_genes = as.integer(gate$highly_variable_genes),
    principal_components = as.integer(gate$principal_components),
    detection_minimum = as.numeric(gate$feature_detection_fraction_minimum),
    detection_maximum = as.numeric(gate$feature_detection_fraction_maximum),
    minimum_features = as.integer(gate$minimum_graph_features)
  )
}

cp_family_graph <- function(log_expression, cells, excluded_genes, gate_spec) {
  cells <- as.logical(cells)
  n <- sum(cells)
  require_true(n >= gate_spec$k + 2L, sprintf("study has %d cells, too few for a %d-neighbor graph", n, gate_spec$k))
  block <- log_expression[, cells, drop = FALSE]
  genes <- rownames(block)
  require_true(!is.null(genes) && !is.null(colnames(block)), "expression must carry gene and cell names")
  detection <- as.numeric(Matrix::rowMeans(block > 0))
  sums <- as.numeric(Matrix::rowSums(block))
  squares <- as.numeric(Matrix::rowSums(block ^ 2))
  variance <- (squares - sums ^ 2 / n) / (n - 1)
  eligible <- detection >= gate_spec$detection_minimum & detection <= gate_spec$detection_maximum &
    is.finite(variance) & variance > 0 & !(genes %in% excluded_genes)
  require_true(sum(eligible) >= gate_spec$minimum_features, sprintf("only %d eligible genes for the study's neighborhood graph", sum(eligible)))
  ordered <- order(-variance, genes, method = "radix")
  ordered <- ordered[eligible[ordered]]
  features <- genes[ordered[seq_len(min(gate_spec$highly_variable_genes, length(ordered)))]]
  dense <- t(as.matrix(block[match(features, genes), , drop = FALSE]))
  centre <- colMeans(dense)
  population_sd <- sqrt(colMeans(sweep(dense, 2L, centre) ^ 2))
  population_sd[population_sd <= 0] <- 1
  dense <- sweep(sweep(dense, 2L, centre), 2L, population_sd, "/")
  dense[!is.finite(dense)] <- 0
  rank <- min(gate_spec$principal_components, ncol(dense), n - 1L)
  decomposition <- svd(dense, nu = rank, nv = 0L)
  scores <- decomposition$u %*% diag(decomposition$d[seq_len(rank)], rank, rank)
  distance <- as.matrix(stats::dist(scores)) ^ 2
  diag(distance) <- Inf
  k <- min(gate_spec$k, n - 1L)
  neighbours <- t(apply(distance, 1L, function(row) order(row)[seq_len(k)]))
  list(
    cells = colnames(block), neighbours = neighbours, k = k,
    features = features, principal_components = rank,
    excluded_genes = sort(unique(excluded_genes)),
    eligible_features = sum(eligible)
  )
}

cp_local_enrichment_gate <- function(candidate, graph, q_threshold) {
  candidate <- as.logical(candidate)
  candidate[is.na(candidate)] <- FALSE
  neighbours <- graph$neighbours
  require_true(length(candidate) == nrow(neighbours), "candidate vector and graph size differ")
  total <- length(candidate)
  positives <- sum(candidate)
  k <- ncol(neighbours)
  hits <- as.integer(rowSums(matrix(candidate[neighbours], nrow = nrow(neighbours))))
  p_value <- rep(NA_real_, total)
  q_value <- rep(NA_real_, total)
  if (positives > 0L) {
    p_value[candidate] <- stats::phyper(hits[candidate] - 1L, positives - 1L, (total - 1L) - (positives - 1L), k, lower.tail = FALSE)
    q_value[candidate] <- stats::p.adjust(p_value[candidate], method = "BH")
  }
  pass <- candidate & !is.na(q_value) & q_value <= q_threshold
  list(pass = pass, hits = hits, p_value = p_value, q_value = q_value, candidates = positives, k = k)
}

cp_two_stage_call <- function(score, detected_genes, detection_requirement, threshold,
                              relaxation, graph, gate_spec) {
  require_true(is.finite(threshold) && is.finite(relaxation) && relaxation >= 0, "two-stage call needs a finite threshold and non-negative relaxation")
  stage_one <- detected_genes >= detection_requirement & is.finite(score) & score >= (threshold - relaxation)
  gate <- cp_local_enrichment_gate(stage_one, graph, gate_spec$q_threshold)
  list(
    stage_one = stage_one, gate_pass = gate$pass, call = gate$pass,
    hits = gate$hits, q_value = gate$q_value, candidates = gate$candidates
  )
}

cp_training_feature_statistics <- function(log_expression, genes, metadata, training_cells,
                                           eligible_genes) {
  index <- match(eligible_genes, genes)
  require_true(!anyNA(index), "eligible genes are absent from expression")
  weights <- cp_family_balanced_weights(metadata, training_cells)
  expression <- log_expression[index, training_cells, drop = FALSE]
  moments <- cp_weighted_row_moments(expression, weights)
  detection <- as.numeric((expression > 0) %*% weights)
  data.frame(
    gene_symbol = eligible_genes,
    mean_expression = moments$mean,
    detection_frequency = detection,
    expression_variance = moments$sd ^ 2,
    stringsAsFactors = FALSE
  )
}

cp_rank_bins <- function(values, bins = 5L) {
  values <- as.numeric(values)
  require_true(all(is.finite(values)), "matching statistics contain non-finite values")
  p <- (rank(values, ties.method = "average") - 0.5) / length(values)
  pmin(as.integer(bins), pmax(1L, as.integer(ceiling(p * bins))))
}

cp_draw_matched_random_panels <- function(log_expression, genes, metadata, training_cells,
                                          real_panel, eligible_genes, replicates, seed,
                                          forbidden_genes = character(), bins = 5L) {
  require_true(all(real_panel %in% eligible_genes), "real panel is not contained in the eligible matching universe")
  statistics <- cp_training_feature_statistics(
    log_expression, genes, metadata, training_cells, eligible_genes
  )
  statistics$mean_bin <- cp_rank_bins(statistics$mean_expression, bins)
  statistics$detection_bin <- cp_rank_bins(statistics$detection_frequency, bins)
  statistics$variance_bin <- cp_rank_bins(statistics$expression_variance, bins)
  rownames(statistics) <- statistics$gene_symbol
  pool <- setdiff(eligible_genes, unique(c(real_panel, forbidden_genes)))
  require_true(length(pool) >= length(real_panel), "random-panel pool is smaller than the real panel")
  panel_rows <- list()
  audit_rows <- list()
  for (replicate in seq_len(as.integer(replicates))) {
    set.seed(as.integer(seed) + replicate)
    available <- pool
    selected <- character()
    for (gene in real_panel) {
      target <- statistics[gene, c("mean_bin", "detection_bin", "variance_bin"), drop = FALSE]
      candidate <- statistics[available, , drop = FALSE]
      distance <- abs(candidate$mean_bin - target$mean_bin) +
        abs(candidate$detection_bin - target$detection_bin) +
        abs(candidate$variance_bin - target$variance_bin)
      minimum <- min(distance)
      choices <- sort(candidate$gene_symbol[distance == minimum])
      choice <- sample(choices, 1L)
      selected <- c(selected, choice)
      available <- setdiff(available, choice)
      audit_rows[[length(audit_rows) + 1L]] <- data.frame(
        replicate = replicate, real_gene = gene, matched_gene = choice,
        adjacent_bin_distance = minimum,
        exact_joint_bin = minimum == 0,
        stringsAsFactors = FALSE
      )
    }
    panel_rows[[replicate]] <- data.frame(
      replicate = replicate, panel_rank = seq_along(selected),
      gene_symbol = selected, stringsAsFactors = FALSE
    )
  }
  list(
    panels = do.call(rbind, panel_rows),
    matching_audit = do.call(rbind, audit_rows),
    training_statistics = statistics
  )
}

cp_downsample_counts <- function(counts, fraction, seed) {
  require_true(inherits(counts, "sparseMatrix"), "downsampling requires a Matrix sparseMatrix")
  require_true(fraction > 0 && fraction < 1, "downsampling fraction must be in (0,1)")
  require_true(all(is.finite(counts@x)) && all(counts@x >= 0) && all(counts@x == floor(counts@x)), "downsampling requires nonnegative integral sparse counts")
  set.seed(as.integer(seed))
  result <- counts
  result@x <- as.numeric(stats::rbinom(length(counts@x), size = counts@x, prob = fraction))
  Matrix::drop0(result)
}
