# The numerical parameters, the two-component Gaussian mixture and the ROC AUC of the 20-gene
# classifier. Sourcing this file reads no data and fits no model.

validate_numerical_contract <- function(spec) {
  parameters <- spec$parameters
  require_true(as.integer(parameters$minimum_cells_per_unit_class) >= 1L, "minimum cells per unit-class must be positive")
  require_true(parameters$training_detection_floor > 0 && parameters$training_detection_floor < 1, "training detection floor must be in (0,1)")
  require_true(as.integer(parameters$minimum_training_families) >= 2L, "at least two training studies are required")
  invisible(TRUE)
}

read_numerical_contract <- function(path) {
  require_true(requireNamespace("jsonlite", quietly = TRUE), "jsonlite is required to read the numerical parameter file")
  spec <- jsonlite::fromJSON(path, simplifyVector = FALSE)
  validate_numerical_contract(spec)
  spec
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
