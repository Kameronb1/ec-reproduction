# Helper functions for s2_classifier_swap.R (Figure S6): the parts of the Figure 2 classifier code it uses (the
# candidate-panel fit without threshold relaxation) and the four-marker classifier's neighborhood check. Their
# counterparts in workflow/scripts/figure2_supp/ are:
#   require_true, abortf                         (classifier_lib/vendor/base/R/figure2_postvalidation_core.R)
#   weighted_quantile_contract, mixture_log_likelihood, fit_equal_variance_gmm   (classifier_lib/vendor/base/R/figure2_postvalidation_analysis.R)
#   cp_* panel fit, scoring and two-stage call   (classifier_lib/R/candidate_panel_core.R)
#   ec_graph_veto                                (ec_classifier_lib.R; the four-marker classifier's neighborhood check)
# Only the functions the calculation uses are kept; nothing here reads or writes files.

abortf <- function(format, ...) {
  stop(sprintf(format, ...), call. = FALSE)
}

require_true <- function(value, message) {
  if (!isTRUE(value)) abortf("%s", message)
  invisible(TRUE)
}

weighted_quantile_contract <- function(values, weights, probability) {
  require_true(probability >= 0 && probability <= 1, "weighted quantile probability is outside [0,1]")
  keep <- is.finite(values) & is.finite(weights) & weights > 0
  require_true(any(keep) && all(weights[keep] >= 0), "weighted quantile has no valid positive weights")
  values <- as.numeric(values[keep])
  weights <- as.numeric(weights[keep]) / sum(weights[keep])
  order_index <- order(values, seq_along(values), method = "radix")
  ordered_values <- values[order_index]
  ordered_weights <- weights[order_index]
  midpoint_cdf <- cumsum(ordered_weights) - 0.5 * ordered_weights
  as.numeric(stats::approx(midpoint_cdf, ordered_values, xout = probability, rule = 2, ties = "ordered")$y)
}

mixture_log_likelihood <- function(values, weights, mixture, means, variance) {
  log_density <- cbind(
    log(mixture[[1L]]) - 0.5 * log(2 * pi * variance) - 0.5 * (values - means[[1L]])^2 / variance,
    log(mixture[[2L]]) - 0.5 * log(2 * pi * variance) - 0.5 * (values - means[[2L]])^2 / variance
  )
  row_max <- apply(log_density, 1L, max)
  sum(weights * (row_max + log(rowSums(exp(log_density - row_max)))))
}

fit_equal_variance_gmm <- function(values, weights, spec) {
  values <- as.numeric(values)
  weights <- as.numeric(weights)
  keep <- is.finite(values) & is.finite(weights) & weights > 0
  require_true(any(keep) && all(weights[keep] >= 0), "GMM has no valid positive weights")
  values <- values[keep]
  weights <- weights[keep] / sum(weights[keep])
  overall_mean <- sum(values * weights)
  overall_variance <- sum(weights * (values - overall_mean)^2)
  unavailable <- function(reason) list(status = "UNAVAILABLE", reason = reason)
  if (overall_variance <= as.numeric(spec$variance_floor)) return(unavailable("OVERALL_VARIANCE_BELOW_FLOOR"))
  pairs <- spec$deterministic_initial_quantile_pairs
  best <- NULL
  for (init_index in seq_along(pairs)) {
    pair <- as.numeric(unlist(pairs[[init_index]]))
    means <- c(weighted_quantile_contract(values, weights, pair[[1L]]), weighted_quantile_contract(values, weights, pair[[2L]]))
    if (means[[1L]] >= means[[2L]]) next
    mixture <- c(0.5, 0.5)
    variance <- max(overall_variance, as.numeric(spec$variance_floor))
    previous <- mixture_log_likelihood(values, weights, mixture, means, variance)
    converged <- FALSE
    log_likelihood <- previous
    iterations <- 0L
    for (iteration in seq_len(as.integer(spec$max_iterations))) {
      iterations <- iteration
      log_density <- cbind(
        log(mixture[[1L]]) - 0.5 * log(2 * pi * variance) - 0.5 * (values - means[[1L]])^2 / variance,
        log(mixture[[2L]]) - 0.5 * log(2 * pi * variance) - 0.5 * (values - means[[2L]])^2 / variance
      )
      row_max <- apply(log_density, 1L, max)
      exponent <- exp(log_density - row_max)
      responsibilities <- exponent / rowSums(exponent)
      weighted_responsibility <- responsibilities * weights
      component_mass <- colSums(weighted_responsibility)
      if (any(component_mass <= 0)) break
      next_mixture <- component_mass / sum(component_mass)
      next_means <- colSums(weighted_responsibility * values) / component_mass
      centered <- sweep(matrix(values, nrow = length(values), ncol = 2L), 2L, next_means, "-")
      next_variance <- sum(weighted_responsibility * centered^2)
      next_variance <- max(next_variance, as.numeric(spec$variance_floor))
      log_likelihood <- mixture_log_likelihood(values, weights, next_mixture, next_means, next_variance)
      mixture <- next_mixture
      means <- next_means
      variance <- next_variance
      if (abs(log_likelihood - previous) <= as.numeric(spec$tolerance)) {
        converged <- TRUE
        break
      }
      previous <- log_likelihood
    }
    if (!converged) next
    component_order <- order(means, method = "radix")
    means <- means[component_order]
    mixture <- mixture[component_order]
    log_likelihood <- mixture_log_likelihood(values, weights, mixture, means, variance)
    candidate <- list(log_likelihood = log_likelihood, initialization_index = init_index - 1L, mixture = mixture, means = means, variance = variance, iterations = iterations)
    if (is.null(best) || candidate$log_likelihood > best$log_likelihood ||
        (identical(candidate$log_likelihood, best$log_likelihood) && candidate$initialization_index < best$initialization_index)) best <- candidate
  }
  if (is.null(best)) return(unavailable("NO_DETERMINISTIC_START_CONVERGED"))
  shared_sd <- sqrt(best$variance)
  if (any(best$mixture < as.numeric(spec$minimum_component_weight))) return(unavailable("COMPONENT_WEIGHT_BELOW_FLOOR"))
  if (shared_sd < as.numeric(spec$minimum_component_sd)) return(unavailable("COMPONENT_SD_BELOW_FLOOR"))
  separation <- (best$means[[2L]] - best$means[[1L]]) / shared_sd
  if (separation < as.numeric(spec$minimum_mean_separation_pooled_sd)) return(unavailable("COMPONENT_MEANS_INSUFFICIENTLY_SEPARATED"))
  threshold <- 0.5 * sum(best$means) + best$variance * log(best$mixture[[1L]] / best$mixture[[2L]]) / diff(best$means)
  if (!is.finite(threshold) || threshold < best$means[[1L]] || threshold > best$means[[2L]]) return(unavailable("NO_EQUAL_POSTERIOR_ROOT_BETWEEN_MEANS"))
  list(
    status = "PASS", reason = "", threshold = threshold,
    lower_mean = best$means[[1L]], upper_mean = best$means[[2L]],
    lower_weight = best$mixture[[1L]], upper_weight = best$mixture[[2L]],
    shared_sd = shared_sd, separation_shared_sd = separation,
    weighted_log_likelihood = best$log_likelihood,
    initialization_index = best$initialization_index, iterations = best$iterations
  )
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
  require_true(length(cells) == nrow(metadata) && any(cells), "family weights require selected cells")
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
  # Missing row names mean the Matrix package is not loaded (a setup problem), not that genes are missing from the data.
  require_true(!is.null(rownames(log_expression)),
               paste0("expression matrix has no row names: the Matrix package is not attached, so S4 ",
                      "dispatch failed. This is a setup problem, not missing genes. ",
                      "Load the Matrix package before any panel fit."))
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

cp_score_panel_model <- function(log_expression, model) {
  standardized <- cp_apply_standardizer(log_expression, model$standardizer)
  score <- as.numeric(Matrix::crossprod(model$gene_weights, standardized))
  index <- match(model$panel_genes, rownames(log_expression))
  require_true(!anyNA(index), "trained panel genes are missing from evaluation expression")
  detected <- as.numeric(Matrix::colSums(log_expression[index, , drop = FALSE] > 0))
  call <- detected >= model$detection_requirement & score >= model$gmm$threshold
  list(score = score, detected_genes = detected, call = call)
}

cp_fit_panel_model <- function(log_expression, genes, metadata, labels, training_cells,
                               panel_genes, gmm_spec, detection_fraction = 0.25,
                               null_panel_policy = FALSE, prepared = NULL) {
  if (is.null(prepared)) {
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
  } else {
    require_true(identical(prepared$training_cells,training_cells) &&
                 identical(prepared$cell_ids,metadata$project_cell_id),"panel cache population differs")
    ix <- match(panel_genes,prepared$panel_genes)
    require_true(!anyNA(ix),"panel cache lacks requested genes")
    standardizer <- prepared$standardizer[ix,,drop=FALSE];rownames(standardizer)<-NULL
    family_effects <- prepared$family_effects[ix]
  }
  raw_weights <- pmax(family_effects, 0)
  raw_weights[!is.finite(raw_weights)] <- 0
  weight_rule <- "POSITIVE_TRAINING_FAMILY_EFFECT_NORMALIZED_TO_SUM_ONE"
  if (!any(raw_weights > 0)) {
    require_true(isTRUE(null_panel_policy), "panel has no positive training-family effect weights")
    raw_weights <- rep(1, length(panel_genes))
    weight_rule <- "EQUAL_WEIGHTS_FALLBACK__NO_POSITIVE_EFFECT__WORST_CASE_NULL"
  }
  weights <- raw_weights / sum(raw_weights)
  standardized <- if(is.null(prepared)) cp_apply_standardizer(log_expression, standardizer) else prepared$standardized[ix,,drop=FALSE]
  scores <- as.numeric(Matrix::crossprod(weights, standardized))
  panel_index <- match(panel_genes, rownames(log_expression))
  detected <- if(is.null(prepared)) as.numeric(Matrix::colSums(log_expression[panel_index, , drop = FALSE] > 0)) else as.numeric(Matrix::colSums(prepared$detected[ix,,drop=FALSE]))
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

EC_VETO_Q             <- 0.05   # a fixed tuning threshold, not false-discovery-rate control

ec_graph_veto <- function(candidate, graphs, veto_q = EC_VETO_Q) {
  local_q <- rep(NA_real_, length(candidate))
  for (g in graphs) {
    if (is.null(g)) next
    is_cand_a <- candidate[g$cols]
    S <- sum(is_cand_a)
    if (S < 2L) next
    hits <- rowSums(matrix(is_cand_a[g$nn], nrow = g$n_a, ncol = g$k))
    foc <- which(is_cand_a)
    p <- phyper(hits[foc] - 1L, S - 1L, (g$n_a - 1L) - (S - 1L), g$k, lower.tail = FALSE)
    local_q[g$cols[foc]] <- p.adjust(p, method = "BH")
  }
  candidate & !is.na(local_q) & local_q <= veto_q
}
