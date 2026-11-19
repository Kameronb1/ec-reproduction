# Shared functions of the four-marker EC classifier: the two-component Gaussian mixture that sets
# the score cutoff, and the within-group nearest-neighbor graph used for the local-consistency check.
# Sourced by the analysis scripts of this folder; nothing here reads or writes files.


EC_GMM_VAR_FLOOR      <- 0.0025
EC_GMM_INIT_QUANTILES <- list(c(0.1, 0.9), c(0.25, 0.75), c(0.4, 0.8))
EC_GMM_TOL            <- 1e-8
EC_GMM_MAX_ITER       <- 1000L

ec_gmm_wquantile <- function(x, w, q) {
  o <- order(x)                       # stable on ties
  xs <- x[o]; ws <- w[o]
  cdf <- cumsum(ws) - 0.5 * ws
  approx(cdf, xs, xout = q, rule = 2)$y
}

ec_gmm_loglik <- function(x, w, mu, mix, v) {
  log_dens <- cbind(log(mix[1]) - 0.5 * log(2 * pi * v) - 0.5 * (x - mu[1])^2 / v,
                    log(mix[2]) - 0.5 * log(2 * pi * v) - 0.5 * (x - mu[2])^2 / v)
  m <- pmax(log_dens[, 1], log_dens[, 2])
  sum(w * (m + log(rowSums(exp(log_dens - m)))))
}

ec_gmm_fit <- function(x, w = rep(1 / length(x), length(x))) {
  stopifnot(length(x) == length(w), all(is.finite(x)))
  overall_mean <- sum(w * x)
  overall_var  <- sum(w * (x - overall_mean)^2)
  if (overall_var <= EC_GMM_VAR_FLOOR) stop("score variance below the floor; mixture refused")

  best <- NULL
  for (init_index in seq_along(EC_GMM_INIT_QUANTILES)) {
    qs <- EC_GMM_INIT_QUANTILES[[init_index]]
    mu <- c(ec_gmm_wquantile(x, w, qs[1]), ec_gmm_wquantile(x, w, qs[2]))
    if (mu[1] >= mu[2]) next
    mix <- c(0.5, 0.5)
    v <- max(overall_var, EC_GMM_VAR_FLOOR)
    prev_ll <- NA_real_; converged <- FALSE; iters <- 0L

    for (it in seq_len(EC_GMM_MAX_ITER)) {
      log_dens <- cbind(log(mix[1]) - 0.5 * log(2 * pi * v) - 0.5 * (x - mu[1])^2 / v,
                        log(mix[2]) - 0.5 * log(2 * pi * v) - 0.5 * (x - mu[2])^2 / v)
      m <- pmax(log_dens[, 1], log_dens[, 2])
      e <- exp(log_dens - m)
      resp <- e / rowSums(e)
      wr <- w * resp
      mass <- colSums(wr)
      mix <- mass / sum(mass)
      mu <- c(sum(wr[, 1] * x), sum(wr[, 2] * x)) / mass
      v <- max(sum(wr[, 1] * (x - mu[1])^2) + sum(wr[, 2] * (x - mu[2])^2), EC_GMM_VAR_FLOOR)

      ll <- ec_gmm_loglik(x, w, mu, mix, v)     # recomputed with the NEW parameters
      iters <- it
      if (!is.na(prev_ll) && abs(ll - prev_ll) <= EC_GMM_TOL) { converged <- TRUE; break }
      prev_ll <- ll
    }
    if (!converged) next

    o <- order(mu)                              # components ascending by mean
    mu <- mu[o]; mix <- mix[o]
    cand <- list(mu = mu, mix = mix, var = v, loglik = ec_gmm_loglik(x, w, mu, mix, v),
                 init_index = init_index - 1L, iterations = iters)
    # best of three by log-likelihood; ties keep the EARLIER initialization
    if (is.null(best) || cand$loglik > best$loglik) best <- cand
  }
  if (is.null(best)) stop("no initialization converged")

  sd_shared <- sqrt(best$var)
  if (min(best$mix) < 0.02)                              stop("degenerate mixing weight")
  if (sd_shared < 0.05)                                  stop("degenerate component width")
  if ((best$mu[2] - best$mu[1]) / sd_shared < 0.5)       stop("components not separated")

  # Equal-posterior boundary: where a cell is equally likely to belong to either component.
  thr <- 0.5 * (best$mu[1] + best$mu[2]) +
    best$var * log(best$mix[1] / best$mix[2]) / (best$mu[2] - best$mu[1])
  if (thr < best$mu[1] || thr > best$mu[2]) stop("threshold outside the component means")

  list(mu_low = best$mu[1], mu_high = best$mu[2],
       prop_low = best$mix[1], prop_high = best$mix[2],
       var = best$var, sd = sd_shared, loglik = best$loglik,
       init_index = best$init_index, iterations = best$iterations,
       threshold = thr)
}


EC_GRAPH_NEIGHBORS   <- 25L
EC_GRAPH_COMPONENTS   <- 20L
EC_GRAPH_TOP_HVG      <- 250L
EC_GRAPH_DETECT_RANGE <- c(0.01, 0.99)
EC_VETO_Q             <- 0.05   # a fixed cutoff on the local check's q-values, not false-discovery-rate control

ec_graph_build_one <- function(X, cols, excluded_rows = integer(0)) {
  n_a <- length(cols)
  if (n_a <= EC_GRAPH_NEIGHBORS) return(NULL)
  A <- X[, cols, drop = FALSE]

  detect_frac <- rowSums(A > 0) / n_a
  rmeans <- rowMeans(A)
  rvars <- rowSums((A - rmeans)^2) / n_a
  eligible <- which(detect_frac >= EC_GRAPH_DETECT_RANGE[1] &
                    detect_frac <= EC_GRAPH_DETECT_RANGE[2] & rvars > 0)
  eligible <- setdiff(eligible, excluded_rows)
  stopifnot(length(eligible) > 0)
  chosen <- eligible[order(rvars[eligible], decreasing = TRUE)][
    seq_len(min(EC_GRAPH_TOP_HVG, length(eligible)))]

  W <- t(A[chosen, , drop = FALSE])                       # cells x genes
  W <- scale(W, center = TRUE, scale = apply(W, 2, function(v) sqrt(mean((v - mean(v))^2))))
  W[!is.finite(W)] <- 0

  k_comp <- min(EC_GRAPH_COMPONENTS, ncol(W), n_a - 1L)
  sv <- svd(W, nu = k_comp, nv = 0)
  P <- sv$u[, seq_len(k_comp), drop = FALSE] %*% diag(sv$d[seq_len(k_comp)], nrow = k_comp)

  D <- as.matrix(dist(P))^2
  diag(D) <- Inf
  k <- min(EC_GRAPH_NEIGHBORS, n_a - 1L)
  nn <- t(apply(D, 1, function(row) order(row)[seq_len(k)]))   # within-group indices
  list(cols = cols, n_a = n_a, k = k, nn = nn, hvgs = length(chosen), components = k_comp,
       P = P)
}

ec_graphs_build <- function(X, arm, excluded_rows = integer(0)) {
  arms <- sort(unique(arm))
  g <- lapply(arms, function(a) ec_graph_build_one(X, which(arm == a), excluded_rows))
  names(g) <- arms
  g
}

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
