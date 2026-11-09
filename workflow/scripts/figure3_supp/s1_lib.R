# Shared functions for the Figure S15 scripts (s1_marker_limma.R, s1_marker_stability.R). They rank genes and score
# panels the way the Figure 2 classifier does: genes are ranked by limma's moderated t for SC-SIC versus other SC
# endocrine cells, fitted on sample-by-class mean expression, paired within samples, each sample weighted
# 1 / (samples in its study).

suppressPackageStartupMessages({ library(Matrix); library(limma) })

require_true <- function(value, message) if (!isTRUE(value)) stop(message, call. = FALSE)

# Figure 2's fixed numerical settings (the three values these analyses use)
PARAMETERS <- list(training_detection_floor = 0.05, minimum_training_families = 3L, minimum_cells_per_unit_class = 5L)
# genes never ranked: the four canonical EC marker genes
EXCLUDED_MARKERS <- c("FEV", "LMX1A", "SLC18A1", "TPH1")

# Read the SC endocrine count package: ln(1 + count / full-object cell total x 1e6), genes x cells.
read_sc_endocrine <- function(meta_path, genes_path, indptr_path, indices_path, data_path) {
  genes <- readLines(genes_path)
  meta <- read.delim(meta_path, stringsAsFactors = FALSE, colClasses = "character")
  meta$total_umis <- as.numeric(meta$total_umis)
  nG <- length(genes); nC <- nrow(meta)
  rd <- function(f, n) readBin(f, what = "integer", n = n, size = 4L, endian = "little")
  p <- rd(indptr_path, nC + 1L); nnz <- p[nC + 1L]
  i <- rd(indices_path, nnz); x <- as.numeric(rd(data_path, nnz))
  require_true(length(i) == nnz && length(x) == nnz, "sparse arrays truncated")
  x <- log1p(x / rep(meta$total_umis, times = diff(p)) * 1e6)
  require_true(all(is.finite(x)), "non-finite normalized value")
  expr <- new("dgCMatrix", i = i, p = p, x = x, Dim = c(nG, nC), Dimnames = list(genes, meta$project_cell_id))
  rm(i, x); invisible(gc())
  list(expr = expr, genes = genes, meta = meta)
}

training_feature_mask <- function(log_expression, genes, metadata, train, detection_floor, minimum_training_families, excluded_genes) {
  families <- sort(unique(metadata$family_id[train]))
  pass_count <- integer(length(genes))
  for (family in families) {
    cells <- train & metadata$family_id == family
    pass_count <- pass_count + (Matrix::rowMeans(log_expression[, cells, drop = FALSE] > 0) >= detection_floor)
  }
  eligible <- pass_count >= as.integer(minimum_training_families)
  eligible <- eligible & !(genes %in% excluded_genes)
  eligible <- eligible & !grepl("^MT-", genes, ignore.case = TRUE)
  eligible <- eligible & !grepl("^RP[SL][0-9]", genes, ignore.case = TRUE)
  eligible
}

binary_auc <- function(response, score) {
  keep <- is.finite(score) & response %in% c(0L, 1L)
  response <- response[keep]
  score <- score[keep]
  positives <- sum(response == 1L)
  negatives <- sum(response == 0L)
  if (positives == 0L || negatives == 0L) return(NA_real_)
  ranks <- rank(score, ties.method = "average")
  (sum(ranks[response == 1L]) - positives * (positives + 1) / 2) / (positives * negatives)
}

cp_align_labels <- function(metadata, labels) {
  index <- match(metadata$project_cell_id, labels$project_cell_id)
  require_true(!anyNA(index), "labels do not cover every metadata cell")
  labels[index, , drop = FALSE]
}

cp_family_balanced_weights <- function(metadata, cells) {
  cells <- as.logical(cells)
  family <- metadata$family_id[cells]
  size <- table(family)
  weights <- 1 / as.numeric(size[family])
  weights / sum(weights)
}

cp_weighted_row_moments <- function(expression, weights) {
  weights <- weights / sum(weights)
  mean <- as.numeric(expression %*% weights)
  second <- as.numeric((expression ^ 2) %*% weights)
  variance <- pmax(second - mean ^ 2, 0)
  list(mean = mean, sd = sqrt(variance))
}

cp_fit_standardizer <- function(log_expression, metadata, training_cells, genes) {
  index <- match(genes, rownames(log_expression))
  require_true(!anyNA(index), "panel genes are missing from expression")
  weights <- cp_family_balanced_weights(metadata, training_cells)
  moments <- cp_weighted_row_moments(log_expression[index, training_cells, drop = FALSE], weights)
  require_true(all(moments$sd > 0), "at least one panel gene has zero training variance")
  data.frame(gene_symbol = genes, training_mean = moments$mean, training_sd = moments$sd, stringsAsFactors = FALSE)
}

cp_apply_standardizer <- function(log_expression, standardizer) {
  index <- match(standardizer$gene_symbol, rownames(log_expression))
  dense <- as.matrix(log_expression[index, , drop = FALSE])
  sweep(sweep(dense, 1L, standardizer$training_mean, "-"), 1L, standardizer$training_sd, "/")
}

cp_unit_class_means <- function(log_expression, metadata, labels, cells, genes, minimum_cells_per_unit_class = 5L) {
  labels <- cp_align_labels(metadata, labels)
  cells <- as.logical(cells) & labels$label %in% c("POSITIVE", "NEGATIVE")
  gene_index <- match(genes, rownames(log_expression))
  require_true(!anyNA(gene_index), "unit-class genes are missing")
  grouping <- paste(metadata$family_id, metadata$biological_unit_id, labels$label, sep = "::")
  groups <- split(which(cells), grouping[cells])
  keep <- vapply(groups, length, integer(1L)) >= as.integer(minimum_cells_per_unit_class)
  groups <- groups[keep]
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

cp_limma_ec_effect <- function(units) {
  m <- units$metadata
  unit <- factor(paste(m$family_id, m$biological_unit_id, sep = "::"))
  ec <- as.integer(m$label == "POSITIVE")
  units_per_family <- tapply(as.character(unit), m$family_id, function(v) length(unique(v)))
  w <- as.numeric(1 / units_per_family[m$family_id])
  design <- stats::model.matrix(~ 0 + unit + ec)
  fit <- limma::eBayes(limma::lmFit(units$values, design, weights = matrix(w, nrow(units$values), length(w), byrow = TRUE)), trend = TRUE)
  tt <- limma::topTable(fit, coef = "ec", number = Inf, sort.by = "none")
  data.frame(gene_symbol = rownames(units$values), logFC = tt$logFC, t = tt$t, p = tt$P.Value,
             fdr = tt$adj.P.Val, stringsAsFactors = FALSE)
}

cp_training_gene_summary <- function(log_expression, genes, metadata, labels, training_cells,
                                     excluded_genes, detection_floor, minimum_training_families,
                                     minimum_cells_per_unit_class) {
  eligible <- training_feature_mask(log_expression, genes, metadata, training_cells,
                                    detection_floor, minimum_training_families, excluded_genes)
  eligible_genes <- genes[eligible]
  units <- cp_unit_class_means(log_expression, metadata, labels, training_cells, eligible_genes, minimum_cells_per_unit_class)
  families <- sort(unique(units$metadata$family_id))
  effect_blocks <- lapply(families, function(family) {
    block <- units$metadata$family_id == family
    positive <- block & units$metadata$label == "POSITIVE"
    negative <- block & units$metadata$label == "NEGATIVE"
    if (!any(positive) || !any(negative)) return(NULL)
    data.frame(family_id = family, gene_symbol = eligible_genes,
               effect = rowMeans(units$values[, positive, drop = FALSE]) - rowMeans(units$values[, negative, drop = FALSE]),
               stringsAsFactors = FALSE)
  })
  effect_blocks <- effect_blocks[!vapply(effect_blocks, is.null, logical(1L))]
  effects <- do.call(rbind, effect_blocks)
  labels_aligned <- cp_align_labels(metadata, labels)
  detection <- do.call(rbind, lapply(sort(unique(metadata$family_id[training_cells])), function(family) {
    block <- training_cells & metadata$family_id == family & labels_aligned$label %in% c("POSITIVE", "NEGATIVE")
    if (!any(block)) return(NULL)
    data.frame(family_id = family, gene_symbol = eligible_genes,
               detection_fraction = as.numeric(Matrix::rowMeans(log_expression[eligible, block, drop = FALSE] > 0)),
               stringsAsFactors = FALSE)
  }))
  effect_by_gene <- split(effects, effects$gene_symbol)
  detection_by_gene <- split(detection, detection$gene_symbol)
  summary <- do.call(rbind, lapply(eligible_genes, function(gene) {
    effect_block <- effect_by_gene[[gene]]
    detection_block <- detection_by_gene[[gene]]
    data.frame(gene_symbol = gene,
               effect_families = length(unique(effect_block$family_id)),
               positive_effect_families = length(unique(effect_block$family_id[effect_block$effect > 0])),
               median_family_effect = stats::median(effect_block$effect, na.rm = TRUE),
               minimum_family_effect = min(effect_block$effect, na.rm = TRUE),
               median_detection_fraction = stats::median(detection_block$detection_fraction, na.rm = TRUE),
               stringsAsFactors = FALSE)
  }))
  fit <- cp_limma_ec_effect(units)
  summary$limma_logfc <- fit$logFC[match(summary$gene_symbol, fit$gene_symbol)]
  summary$limma_t <- fit$t[match(summary$gene_symbol, fit$gene_symbol)]
  summary$limma_p <- fit$p[match(summary$gene_symbol, fit$gene_symbol)]
  summary$limma_fdr <- fit$fdr[match(summary$gene_symbol, fit$gene_symbol)]
  summary <- summary[order(-summary$limma_t, summary$gene_symbol, method = "radix"), , drop = FALSE]
  rownames(summary) <- NULL
  list(summary = summary, effects = effects, detection = detection, unit_class = units)
}

# Gene ranking as in Figure 2, with no redundancy filter between correlated genes: the ranking is the limma order.
rank_genes <- function(log_expression, genes, metadata, labels, training_cells, excluded_genes = EXCLUDED_MARKERS) {
  s <- cp_training_gene_summary(log_expression, genes, metadata, labels, training_cells, excluded_genes,
                                PARAMETERS$training_detection_floor, PARAMETERS$minimum_training_families,
                                PARAMETERS$minimum_cells_per_unit_class)
  list(ranking = s$summary, unit_class = s$unit_class)
}

cp_rank_bins <- function(values, bins = 5L) {
  values <- as.numeric(values)
  p <- (rank(values, ties.method = "average") - 0.5) / length(values)
  pmin(as.integer(bins), pmax(1L, as.integer(ceiling(p * bins))))
}

cp_training_feature_statistics <- function(log_expression, genes, metadata, training_cells, eligible_genes) {
  index <- match(eligible_genes, genes)
  weights <- cp_family_balanced_weights(metadata, training_cells)
  expression <- log_expression[index, training_cells, drop = FALSE]
  moments <- cp_weighted_row_moments(expression, weights)
  detection <- as.numeric((expression > 0) %*% weights)
  data.frame(gene_symbol = eligible_genes, mean_expression = moments$mean, detection_frequency = detection,
             expression_variance = moments$sd ^ 2, stringsAsFactors = FALSE)
}

# Random gene sets matched to the real panel on training-study mean, detection and variance (quintile bins).
cp_draw_matched_random_panels <- function(log_expression, genes, metadata, training_cells, real_panel, eligible_genes,
                                          replicates, seed, forbidden_genes = character(), bins = 5L) {
  statistics <- cp_training_feature_statistics(log_expression, genes, metadata, training_cells, eligible_genes)
  statistics$mean_bin <- cp_rank_bins(statistics$mean_expression, bins)
  statistics$detection_bin <- cp_rank_bins(statistics$detection_frequency, bins)
  statistics$variance_bin <- cp_rank_bins(statistics$expression_variance, bins)
  rownames(statistics) <- statistics$gene_symbol
  pool <- setdiff(eligible_genes, unique(c(real_panel, forbidden_genes)))
  panel_rows <- list()
  for (replicate in seq_len(as.integer(replicates))) {
    set.seed(as.integer(seed) + replicate)
    available <- pool
    selected <- character()
    for (gene in real_panel) {
      target <- statistics[gene, c("mean_bin", "detection_bin", "variance_bin"), drop = FALSE]
      candidate <- statistics[available, , drop = FALSE]
      distance <- abs(candidate$mean_bin - target$mean_bin) + abs(candidate$detection_bin - target$detection_bin) +
        abs(candidate$variance_bin - target$variance_bin)
      choices <- sort(candidate$gene_symbol[distance == min(distance)])
      choice <- sample(choices, 1L)
      selected <- c(selected, choice)
      available <- setdiff(available, choice)
    }
    panel_rows[[replicate]] <- data.frame(replicate = replicate, panel_rank = seq_along(selected), gene_symbol = selected,
                                          stringsAsFactors = FALSE)
  }
  list(panels = do.call(rbind, panel_rows))
}

# Per-gene quantities of a panel fit on the training cells: standardizer, median training-study effect, standardized values.
cp_prepare_panel_fit <- function(log_expression, metadata, labels, training_cells, panel_genes) {
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
  list(panel_genes = panel_genes, standardizer = standardizer, family_effects = family_effects,
       standardized = cp_apply_standardizer(log_expression, standardizer))
}
