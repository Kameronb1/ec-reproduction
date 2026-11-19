# Figure S3A-B, McCauley 2023 HIO graft rows: the four-marker EC classifier's cutoff sweep and
# read-thinning test for the HIO graft cells, computed as heldout.R and four_marker_depth.R compute
# them for the test studies. Published labels: the authors' Endocrine_1 = EC, Endocrine_2 =
# other terminal EEC. The classes are those of marker_fit.R (fitted constants; within-group graph on
# the group's own cells in the fitting-cell matrix).
#
#   S3A  cutoff sweep (same grid plus the refits' display positions) and the five transported
#          leave-one-fitting-group-out refits (fitted without the seven pediatric Crohn's disease donors,
#          roster), drawn at the positions computed by heldout.R
#   S3B  binomial read thinning (f = 0.5 / 0.25, 3 draws; McCauley = study 5 in the seed formula).
#          Counts are recovered from the log1p-CPM matrix; counts on genes outside the matrix are the
#          remainder of the library total and are thinned with the same probability.
#
# Inputs : matrix, features, cells, meta, roster, fit (fitted_constants.json), stability, hio_classes, hio_cells,
#          refits (transported_refits.tsv from heldout.R)
# Outputs: sweep, refits_out, depth
args <- commandArgs(trailingOnly = TRUE)
a <- setNames(sub("^[^=]*=", "", args), sub("=.*$", "", args))
HERE <- dirname(sub("^--file=", "", grep("^--file=", commandArgs(), value = TRUE)))
source(file.path(HERE, "ec_classifier_lib.R"))
suppressPackageStartupMessages({library(Matrix); library(jsonlite)})

HIO <- "SF_GSE214852_HIO_GRAFT"; HIO_FAM <- "gse214852_hio"; HIO_ARM <- "GSE214852_WT"
HIO_LABEL <- "McCauley 2023 · HIO graft"
gz <- function(p) read.delim(gzfile(p), stringsAsFactors = FALSE)
wq <- function(x, p, quote) write.table(x, if (grepl("[.]gz$", p)) gzfile(p) else p, sep = "\t", row.names = FALSE, quote = quote)
EC_CORE <- c("TPH1", "SLC18A1", "FEV", "LMX1A")
live <- fromJSON(a[["fit"]], simplifyVector = FALSE)
CUT <- as.numeric(live$ec_core_cutoff); MINM <- as.integer(live$min_markers_detected); VQ <- as.numeric(live$neighbourhood$q_threshold)
std <- do.call(rbind, lapply(live$standardisation, as.data.frame)); rownames(std) <- std$gene
MU <- std[EC_CORE, "mean"]; SD <- std[EC_CORE, "sd"]; names(MU) <- names(SD) <- EC_CORE
binary_auc <- function(y, s) { r <- rank(s); n1 <- sum(y == 1L); n0 <- sum(y == 0L); (sum(r[y == 1L]) - n1 * (n1 + 1) / 2) / (n1 * n0) }
thin <- function(X, f, seed) { set.seed(as.integer(seed)); Y <- X; Y@x <- as.numeric(rbinom(length(X@x), size = X@x, prob = f)); drop0(Y) }
CLASS_KEEP <- c(EXPLICIT_EC = "Primary-labeled EC", TERMINAL_NON_EC_EEC = "Other terminal EEC", BROAD_EEC_PARENT_UNRESOLVED = "Unresolved EEC label")
SWEEP <- sort(unique(c(seq(-0.90, 1.50, by = 0.02), CUT)))

# ---- the five leave-one-fitting-group-out refits --------------------------------------------------
stab <- read.delim(a[["stability"]], stringsAsFactors = FALSE); loo <- stab[stab$design == "Leave one arm out", ]
dg <- read.delim(gzfile(a[["meta"]]), stringsAsFactors = FALSE, colClasses = "character")
dg <- dg[dg$model_context == "primary_intestine", ]
ROSTER <- jsonlite::fromJSON(a[["roster"]])                         # the seven Crohn's donors are left out
crohn <- dg$study_family == ROSTER$family_id & dg$biological_unit %in% ROSTER$excluded_units
stopifnot(sum(crohn) == 35L)
dg <- dg[!crohn, ]
stopifnot(nrow(dg) == 1164L)
d_counts <- vapply(paste0("raw_count__", EC_CORE), function(cc) as.numeric(dg[[cc]]), numeric(nrow(dg))); colnames(d_counts) <- EC_CORE
dL <- log1p(d_counts * (1e6 / as.numeric(dg$cell_total_counts))); d_arm <- dg$arm_id
pop_sd <- function(m, mu) sqrt(colMeans(sweep(m, 2, mu, "-")^2))
fit_fold <- function(fit) { mu <- colMeans(dL[fit, , drop = FALSE]); sdv <- pop_sd(dL[fit, , drop = FALSE], mu)
  s <- rowMeans(sweep(sweep(dL, 2, mu, "-"), 2, sdv, "/")); list(mu = mu, sd = sdv, threshold = ec_gmm_fit(s[fit])$threshold) }
FOLDS <- lapply(seq_len(nrow(loo)), function(i) fit_fold(d_arm != loo$held_out[i])); names(FOLDS) <- loo$held_out
REF <- read.delim(a[["refits"]], stringsAsFactors = FALSE)
POS <- unique(REF[c("held_out", "equiv_cutoff", "candidates")]); stopifnot(nrow(POS) == 5L)
SWEEP_ALL <- sort(unique(c(SWEEP, POS$equiv_cutoff)))

# ---- McCauley cells: the fitting-cell matrix columns scored by marker_fit.R ---------------------------
h <- gz(a[["hio_classes"]]); h <- h[h$arm_id == HIO_ARM, ]; stopifnot(nrow(h) == 598L)
ex <- read.delim(a[["hio_cells"]], stringsAsFactors = FALSE)
ex <- ex[match(h$cell_key, ex$figure2_cell_id), ]; stopifnot(!anyNA(ex$figure2_cell_id))
h$author_label_class <- ifelse(ex$author_label == "Endocrine_1", "EXPLICIT_EC", "TERMINAL_NON_EC_EEC")
feat <- read.delim(a[["features"]], stringsAsFactors = FALSE); mcells <- read.delim(a[["cells"]], stringsAsFactors = FALSE)
n_all <- nrow(mcells); n_feat <- nrow(feat)
con <- file(a[["matrix"]], "rb"); Xd <- readBin(con, "double", n = as.numeric(n_all) * n_feat, size = 8, endian = "little"); close(con)
dim(Xd) <- c(n_feat, n_all)
ci <- match(h$cell_key, mcells$actual_cell_key); stopifnot(!anyNA(ci))
Xh <- Xd[, ci, drop = FALSE]; rm(Xd); invisible(gc(verbose = FALSE))
colnames(Xh) <- h$cell_key
gh <- vapply(EC_CORE, function(s) { r <- which(feat$feature_symbol == s); stopifnot(length(r) == 1L); r }, integer(1))
harm <- rep(HIO_ARM, ncol(Xh)); hg <- ec_graphs_build(Xh, harm, integer(0))
stopifnot(identical(unname(ec_graph_veto(as.logical(h$candidate), hg, VQ)), as.logical(h$locally_consistent)))

# ---- S3A ---------------------------------------------------------------------------------------
lch <- Xh[gh, , drop = FALSE]; rownames(lch) <- EC_CORE
sw <- list(); tr <- list()
for (cut in SWEEP_ALL) {
  cand <- (h$ec_core_score >= cut) & (h$markers_detected >= MINM)
  cons <- if (any(cand)) ec_graph_veto(cand, hg, VQ) else cand
  for (k in names(CLASS_KEEP)) { i <- h$author_label_class == k; if (!any(i)) next
    sw[[length(sw) + 1L]] <- data.frame(family = HIO_FAM, study = HIO_LABEL, cutoff = cut, class = unname(CLASS_KEEP[[k]]), n = sum(i),
                                        pre = sum(cand & i), post = sum(cons & i), stringsAsFactors = FALSE) }
}
sw <- do.call(rbind, sw); sw$pre_rate <- 100 * sw$pre / sw$n; sw$post_rate <- 100 * sw$post / sw$n
for (g in names(FOLDS)) {
  fo <- FOLDS[[g]]; s_fold <- colMeans((lch - fo$mu[EC_CORE]) / fo$sd[EC_CORE])
  cand_f <- (s_fold >= fo$threshold) & (h$markers_detected >= MINM)
  cons_f <- if (any(cand_f)) ec_graph_veto(cand_f, hg, VQ) else cand_f
  for (k in names(CLASS_KEEP)) { i <- h$author_label_class == k; if (!any(i)) next
    tr[[length(tr) + 1L]] <- data.frame(held_out = g, family = HIO_FAM, study = HIO_LABEL, raw_threshold = fo$threshold, class = unname(CLASS_KEEP[[k]]),
                                        n = sum(i), pre = sum(cand_f & i), post = sum(cons_f & i), stringsAsFactors = FALSE) }
}
tr <- do.call(rbind, tr); tr$post_rate <- 100 * tr$post / tr$n
tr$equiv_cutoff <- POS$equiv_cutoff[match(tr$held_out, POS$held_out)]; tr$candidates <- POS$candidates[match(tr$held_out, POS$held_out)]
tr$curve_rate <- mapply(function(cl, ct) sw$post_rate[sw$class == cl & sw$cutoff == ct], tr$class, tr$equiv_cutoff)
SWEEP_ORDER <- CLASS_KEEP[c("EXPLICIT_EC", "BROAD_EEC_PARENT_UNRESOLVED", "TERMINAL_NON_EC_EEC")]   # published row order
wq(sw[order(match(sw$class, SWEEP_ORDER), sw$cutoff), ], a[["sweep"]], FALSE)
wq(tr[c("held_out", "family", "study", "raw_threshold", "class", "n", "pre", "post", "post_rate", "curve_rate", "equiv_cutoff", "candidates")],
   a[["refits_out"]], FALSE)

# ---- S3B ---------------------------------------------------------------------------------------
tot_h <- as.numeric(h$gex_total_counts)
raw <- expm1(Xh) * rep(tot_h / 1e6, each = nrow(Xh)); Cn <- round(raw)
stopifnot(max(abs(raw - Cn)) < 1e-4, all(Cn >= 0)); rm(raw)
Cn <- Matrix(Cn, sparse = TRUE); rownames(Cn) <- feat$feature_symbol; colnames(Cn) <- h$cell_key
remainder <- round(tot_h - colSums(Cn))
lab <- h$author_label_class; ev <- lab %in% c("EXPLICIT_EC", "TERMINAL_NON_EC_EEC"); resp <- as.integer(lab[ev] == "EXPLICIT_EC")
runs <- rbind(data.frame(fraction = 1, replicate = 0L), expand.grid(fraction = c(.5, .25), replicate = 1:3))
rows <- list()
for (i in seq_len(nrow(runs))) {
  fr <- runs$fraction[i]; rp <- runs$replicate[i]
  Xt <- if (fr == 1) Cn else thin(Cn, fr, 20261008L + rp + 100L * 5L + round(1000 * fr))
  ext <- if (fr == 1) remainder else rbinom(length(remainder), size = remainder, prob = fr)   # same RNG stream, after thin()
  tot <- colSums(Xt) + ext
  L <- log1p(sweep(Xt, 2, tot, "/") * 1e6)
  score <- colMeans((as.matrix(L[gh, , drop = FALSE]) - MU) / SD); det <- as.integer(colSums(Xt[gh, , drop = FALSE] > 0))
  cand <- score >= CUT & det >= MINM
  call <- if (any(cand)) ec_graph_veto(cand, ec_graphs_build(L, harm, integer(0)), VQ) else cand
  cl <- call[ev]
  rows[[length(rows) + 1L]] <- data.frame(study = HIO, fraction = fr, replicate = rp,
    roc_auc = binary_auc(resp, score[ev]), sensitivity = mean(cl[resp == 1L]), specificity = mean(!cl[resp == 0L]),
    ec_called = sum(cl[resp == 1L]), ec_n = sum(resp == 1L), other_called = sum(cl[resp == 0L]), other_n = sum(resp == 0L),
    median_total_counts = median(tot))
}
wq(do.call(rbind, rows), a[["depth"]], TRUE)
