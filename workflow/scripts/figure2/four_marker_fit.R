# Figure 2A: fit of the label-free four-marker EC classifier on the fitting studies.
#
# What it computes
#   1. Keeps the 1,199 primary-tissue fitting cells (model_context == "primary_intestine").
#   2. Standardizes TPH1, SLC18A1, FEV and LMX1A (and DDC, TAC1, which decide nothing) by their
#      mean and population SD over those cells. The EC marker score is the mean of the four
#      marker z-scores; "markers detected" counts the four markers with a nonzero count.
#   3. Fits a two-component equal-variance Gaussian mixture to the score (three deterministic
#      starts, best by log-likelihood); the cutoff is the equal-posterior boundary.
#   4. Panel A tables: cells per 0.1-wide score bin (bins anchored on the cutoff) by markers
#      detected, and the fitted mixture and its two components in cells per bin.
#   5. Applies the fitted constants unchanged to the hPSC-derived HIO cells (not in the fit) and
#      classifies them as the classifier does everywhere: score >= cutoff, >= 2 of 4 markers detected,
#      and a within-group 25-nearest-neighbor graph holding more candidates than chance
#      (hypergeometric, BH q <= 0.05). These results feed the HIO rows of panels B and C.
#
# Inputs : matrix (log1p-CPM, float64, cells x features, C order), features, cells, meta
# Outputs: hist (A_histogram.tsv), curves (A_mixture_curves.tsv), fit (JSON of fitted
#          constants, used to score the test studies), hio (four-marker results for the HIO cells)
#   Rscript four_marker_fit.R matrix=... features=... cells=... meta=... hist=... curves=... fit=... hio=...
# Optional markers= (comma-separated) and min_detected= refit with other marker genes; the Figure 2 supplement
# fits the TPH1-only classifier with markers=TPH1 min_detected=1 (same cells, standardization, mixture and check).
args <- commandArgs(trailingOnly = TRUE)
a <- setNames(sub("^[^=]*=", "", args), sub("=.*$", "", args))
suppressPackageStartupMessages(library(jsonlite))
source(file.path(dirname(sub("^--file=", "", grep("^--file=", commandArgs(), value = TRUE))), "four_marker_lib.R"))

EC_CORE <- if ("markers" %in% names(a)) strsplit(a[["markers"]], ",")[[1]] else c("TPH1", "SLC18A1", "FEV", "LMX1A")
STANDARDIZED_GENES <- c(EC_CORE, "DDC", "TAC1")
MIN_MARKERS_DETECTED <- if ("min_detected" %in% names(a)) as.integer(a[["min_detected"]]) else 2L
BIN_WIDTH <- 0.10

# ---- load -------------------------------------------------------------------------------------
features <- read.delim(a[["features"]], stringsAsFactors = FALSE)
cells <- read.delim(a[["cells"]], stringsAsFactors = FALSE)
con <- file(a[["matrix"]], "rb")
X <- readBin(con, "double", n = as.numeric(nrow(cells)) * nrow(features), size = 8, endian = "little")
close(con)
stopifnot(length(X) == as.numeric(nrow(cells)) * nrow(features))
dim(X) <- c(nrow(features), nrow(cells))                       # features x cells
meta <- read.delim(gzfile(a[["meta"]]), stringsAsFactors = FALSE)
meta <- meta[!duplicated(meta$capture_qualified_key), ]
meta <- meta[match(cells$actual_cell_key, meta$capture_qualified_key), ]
stopifnot(!anyNA(meta$capture_qualified_key))
rows <- vapply(STANDARDIZED_GENES, function(s) { r <- which(features$feature_symbol == s); stopifnot(length(r) == 1L); r }, 1L)

# ---- 1-3. fit on the 1,199 primary-tissue cells ------------------------------------------------
fit_cells <- meta$model_context == "primary_intestine"
Xf <- X[, fit_cells, drop = FALSE]
stopifnot(ncol(Xf) == 1199L)
consts <- t(vapply(rows, function(r) { v <- Xf[r, ]; c(mean = mean(v), sd = sqrt(mean((v - mean(v))^2))) }, numeric(2)))
rownames(consts) <- STANDARDIZED_GENES
Z <- vapply(STANDARDIZED_GENES, function(s) (Xf[rows[[s]], ] - consts[s, "mean"]) / consts[s, "sd"], numeric(ncol(Xf)))
score <- rowMeans(Z[, EC_CORE, drop = FALSE])
detected <- colSums(Xf[rows[EC_CORE], , drop = FALSE] > 0)
gmm <- gmm_fit(score)
cutoff <- gmm$threshold
cat(sprintf("fit: %d cells, cutoff %.17g, %d at or above\n", length(score), cutoff, sum(score >= cutoff)))

fit_out <- list(n_cells = length(score), ec_core_genes = EC_CORE,
                standardisation = data.frame(gene = STANDARDIZED_GENES, mean = consts[, "mean"], sd = consts[, "sd"]),
                ec_core_cutoff = cutoff, min_markers_detected = MIN_MARKERS_DETECTED, q_threshold = LOCAL_Q_MAX,
                gaussian_mixture = gmm)
write_json(fit_out, a[["fit"]], digits = I(17), auto_unbox = TRUE, pretty = TRUE)
back <- fromJSON(a[["fit"]])
stopifnot(identical(back$ec_core_cutoff, cutoff), identical(back$standardisation$mean, unname(consts[, "mean"])),
          identical(back$standardisation$sd, unname(consts[, "sd"])))

# ---- 4. panel A tables -------------------------------------------------------------------------
# The per-cell scores are rounded to 15 significant digits (the precision of a text table) before
# binning, so the bins match the published panel A tables.
s15 <- as.numeric(sprintf("%.15g", score))
breaks <- c(rev(seq(cutoff, min(score) - BIN_WIDTH, by = -BIN_WIDTH)),
            seq(cutoff + BIN_WIDTH, max(score) + BIN_WIDTH, by = BIN_WIDTH))
mid <- breaks[-length(breaks)] + BIN_WIDTH / 2
bin <- cut(s15, breaks = breaks, right = FALSE, include.lowest = TRUE, labels = FALSE)
stopifnot(!anyNA(bin))
A <- aggregate(list(cells = rep(1L, length(bin))), list(bin = bin, markers = detected), sum)
A$x <- as.numeric(sprintf("%.15g", mid))[A$bin]
A$markers <- factor(A$markers, levels = 4:0)
xx <- seq(min(s15) - .15, max(s15) + .15, length.out = 600)
lo <- length(s15) * BIN_WIDTH * (1 - gmm$prop_high) * dnorm(xx, gmm$mu_low, gmm$sd)
hi <- length(s15) * BIN_WIDTH * gmm$prop_high * dnorm(xx, gmm$mu_high, gmm$sd)
curves <- data.frame(x = rep(xx, 3), y = c(lo + hi, lo, hi),
                     curve = factor(rep(c("Mixture", "Lower", "Upper"), each = length(xx)), levels = c("Mixture", "Lower", "Upper")))
write.table(A, a[["hist"]], sep = "\t", row.names = FALSE, quote = TRUE)
write.table(curves, a[["curves"]], sep = "\t", row.names = FALSE, quote = TRUE)

# ---- 5. the HIO cells, scored with the fitted constants ----------------------------------------
hio <- meta$model_context == "hpsc_derived_hio"
Xh <- X[, hio, drop = FALSE]; rm(X, Xf); invisible(gc(verbose = FALSE))
stopifnot(ncol(Xh) == 755L)
Zh <- vapply(STANDARDIZED_GENES, function(s) (Xh[rows[[s]], ] - consts[s, "mean"]) / consts[s, "sd"], numeric(ncol(Xh)))
h_score <- rowMeans(Zh[, EC_CORE, drop = FALSE])
h_det <- colSums(Xh[rows[EC_CORE], , drop = FALSE] > 0)
h_cand <- h_score >= cutoff & h_det >= MIN_MARKERS_DETECTED
h_high <- locally_consistent(h_cand, Xh, meta$arm_id[hio])
cls <- ifelse(h_high, "HIGH_CONFIDENCE_EC", ifelse(h_cand, "LOW_CONFIDENCE_EC", "INSUFFICIENT_EC_EVIDENCE"))
out <- data.frame(cell_key = cells$actual_cell_key[hio], arm_id = meta$arm_id[hio],
                  gex_total_counts = as.numeric(meta$cell_total_counts[hio]), ec_core_score = h_score,
                  markers_detected = h_det, presentation_class = cls)
write.table(out, a[["hio"]], sep = "\t", row.names = FALSE, quote = FALSE)
print(table(out$arm_id, out$presentation_class))
