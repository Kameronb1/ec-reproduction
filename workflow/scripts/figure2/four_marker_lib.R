# Shared functions of the label-free four-marker EC classifier (Figure 2A-C): the two-component
# equal-variance Gaussian mixture that sets the cutoff, and the within-group 25-nearest-neighbor
# local-consistency test that a candidate must pass to be classified as high-confidence EC.
# Sourced by four_marker_fit.R and heldout_marker_calls.R.
LOCAL_Q_MAX <- 0.05                  # a candidate is locally consistent at BH q <= this value

# ---- two-component equal-variance Gaussian mixture --------------------------------------------
gmm_wquantile <- function(x, w, q) {
  o <- order(x); xs <- x[o]; ws <- w[o]
  approx(cumsum(ws) - 0.5 * ws, xs, xout = q, rule = 2)$y
}
gmm_loglik <- function(x, w, mu, mix, v) {
  ld <- cbind(log(mix[1]) - 0.5 * log(2 * pi * v) - 0.5 * (x - mu[1])^2 / v,
              log(mix[2]) - 0.5 * log(2 * pi * v) - 0.5 * (x - mu[2])^2 / v)
  m <- pmax(ld[, 1], ld[, 2])
  sum(w * (m + log(rowSums(exp(ld - m)))))
}
gmm_fit <- function(x, w = rep(1 / length(x), length(x))) {
  floor_v <- 0.0025
  overall_var <- sum(w * (x - sum(w * x))^2)
  best <- NULL
  for (i in 1:3) {
    qs <- list(c(0.1, 0.9), c(0.25, 0.75), c(0.4, 0.8))[[i]]
    mu <- c(gmm_wquantile(x, w, qs[1]), gmm_wquantile(x, w, qs[2]))
    if (mu[1] >= mu[2]) next
    mix <- c(0.5, 0.5); v <- max(overall_var, floor_v); prev <- NA_real_; ok <- FALSE
    for (it in 1:1000) {
      ld <- cbind(log(mix[1]) - 0.5 * log(2 * pi * v) - 0.5 * (x - mu[1])^2 / v,
                  log(mix[2]) - 0.5 * log(2 * pi * v) - 0.5 * (x - mu[2])^2 / v)
      e <- exp(ld - pmax(ld[, 1], ld[, 2]))
      wr <- w * (e / rowSums(e))
      mass <- colSums(wr)
      mix <- mass / sum(mass)
      mu <- c(sum(wr[, 1] * x), sum(wr[, 2] * x)) / mass
      v <- max(sum(wr[, 1] * (x - mu[1])^2) + sum(wr[, 2] * (x - mu[2])^2), floor_v)
      ll <- gmm_loglik(x, w, mu, mix, v)
      if (!is.na(prev) && abs(ll - prev) <= 1e-8) { ok <- TRUE; break }
      prev <- ll
    }
    if (!ok) next
    o <- order(mu); mu <- mu[o]; mix <- mix[o]
    cand <- list(mu = mu, mix = mix, var = v, loglik = gmm_loglik(x, w, mu, mix, v))
    if (is.null(best) || cand$loglik > best$loglik) best <- cand
  }
  stopifnot(!is.null(best), min(best$mix) >= 0.02, sqrt(best$var) >= 0.05,
            (best$mu[2] - best$mu[1]) / sqrt(best$var) >= 0.5)
  thr <- 0.5 * (best$mu[1] + best$mu[2]) + best$var * log(best$mix[1] / best$mix[2]) / (best$mu[2] - best$mu[1])
  list(mu_low = best$mu[1], mu_high = best$mu[2], prop_low = best$mix[1], prop_high = best$mix[2],
       var = best$var, sd = sqrt(best$var), threshold = thr)
}

# ---- within-group neighborhood graph and local-consistency test ------------------------------
graph_one <- function(X, cols) {
  n <- length(cols); if (n <= 25L) return(NULL)
  A <- X[, cols, drop = FALSE]
  det <- rowSums(A > 0) / n; rm <- rowMeans(A); rv <- rowSums((A - rm)^2) / n
  eligible <- which(det >= 0.01 & det <= 0.99 & rv > 0)
  chosen <- eligible[order(rv[eligible], decreasing = TRUE)][seq_len(min(250L, length(eligible)))]
  W <- t(A[chosen, , drop = FALSE])
  W <- scale(W, center = TRUE, scale = apply(W, 2, function(v) sqrt(mean((v - mean(v))^2))))
  W[!is.finite(W)] <- 0
  kc <- min(20L, ncol(W), n - 1L)
  sv <- svd(W, nu = kc, nv = 0)
  D <- as.matrix(dist(sv$u[, seq_len(kc), drop = FALSE] %*% diag(sv$d[seq_len(kc)], nrow = kc)))^2
  diag(D) <- Inf
  k <- min(25L, n - 1L)
  list(cols = cols, n = n, k = k, nn = t(apply(D, 1, function(r) order(r)[seq_len(k)])))
}
locally_consistent <- function(candidate, X, group) {
  q <- rep(NA_real_, length(candidate))
  for (g in sort(unique(group))) {
    G <- graph_one(X, which(group == g)); if (is.null(G)) next   # a group of <= 25 cells has no graph
    cand <- candidate[G$cols]; S <- sum(cand); if (S < 2L) next
    hits <- rowSums(matrix(cand[G$nn], nrow = G$n, ncol = G$k)); foc <- which(cand)
    p <- phyper(hits[foc] - 1L, S - 1L, (G$n - 1L) - (S - 1L), G$k, lower.tail = FALSE)
    q[G$cols[foc]] <- p.adjust(p, method = "BH")
  }
  candidate & !is.na(q) & q <= LOCAL_Q_MAX
}
