# Figure S3A: the four-marker EC classifier applied unchanged to the three intestinal test studies with
# published labels (Elmentaite 2021, Hickey 2023, Egozi 2021).
#
# What it computes, per study
#   - log1p counts per million from raw counts (Gene Expression features only); four-marker score
#     with the fitted constants of marker_fit.R; candidate = score >= cutoff and >= 2 of 4 markers;
#     within-group 25-neighbor graph on all genes; local-consistency check (q <= 0.05).
#   - Share of each published-label class classified as EC over cutoffs -0.90..1.50 (step 0.02,
#     plus the cutoff used), with the graph reused and the check recomputed at each cutoff.
#   - Each leave-one-fitting-group-out refit (standardization and cutoff, from the fitting cells without
#     the seven pediatric Crohn's disease donors, roster) is applied to the test-study cells and placed at
#     the highest cutoff on the full fit's score scale that gives the same candidate count, the three test
#     studies pooled. The sweep also holds the full fit evaluated exactly at those positions.
#
# Inputs : h5ad_dir (one <study>.h5ad per test study), fit (fitted_constants.json from marker_fit.R),
#          stability (from cutoff_stability.R), meta (raw marker counts of the fitting cells), roster
# Outputs: scoring (every test-study cell), sweep, refits
args <- commandArgs(trailingOnly = TRUE)
a <- setNames(sub("^[^=]*=", "", args), sub("=.*$", "", args))
HERE <- dirname(sub("^--file=", "", grep("^--file=", commandArgs(), value = TRUE)))
source(file.path(HERE, "ec_classifier_lib.R")); source(file.path(HERE, "read_h5ad.R"))
suppressPackageStartupMessages({library(Matrix); library(jsonlite)})

EC_CORE <- c("TPH1", "SLC18A1", "FEV", "LMX1A")
STANDARDISED_GENES <- c(EC_CORE, "DDC", "TAC1")
CLASS_HIGH <- "HIGH_CONFIDENCE_EC"; CLASS_LOW <- "LOW_CONFIDENCE_EC"; CLASS_NONE <- "INSUFFICIENT_EC_EVIDENCE"
FAMILIES <- c(gut_cell_atlas = "Gut Cell Atlas 2021", hickey = "Hickey 2023 adult EEC", egozi = "Egozi 2021")
EXPECTED_ARMS <- list(gut_cell_atlas = c("E-MTAB-9536" = 780L, "E-MTAB-9543" = 117L), hickey = c("HICKEY_2023_DRYAD" = 3533L),
                      egozi = c("EGOZI_2021_ZENODO" = 784L))
STUDY_LABEL <- c(gut_cell_atlas = "Elmentaite 2021 · fetal + adult intestine", hickey = "Hickey 2023 · adult intestine",
                 egozi = "Egozi 2021 · fetal + neonatal intestine")
SWEEP_FAMILIES <- c("gut_cell_atlas", "hickey", "egozi")
CLASS_KEEP <- c(EXPLICIT_EC = "Primary-labeled EC", TERMINAL_NON_EC_EEC = "Other terminal EEC",
                BROAD_EEC_PARENT_UNRESOLVED = "Unresolved EEC label")
cpm <- function(counts, total) log1p(sweep(counts, 2, total, "/") * 1e6)

fitj <- jsonlite::fromJSON(a[["fit"]], simplifyVector = FALSE)
std <- do.call(rbind, lapply(fitj$standardisation, function(r) data.frame(gene = r$gene, mean = as.numeric(r$mean), sd = as.numeric(r$sd))))
rownames(std) <- std$gene
LIVE_MU <- setNames(std[STANDARDISED_GENES, "mean"], STANDARDISED_GENES)
LIVE_SD <- setNames(std[STANDARDISED_GENES, "sd"], STANDARDISED_GENES)
LIVE_CUT <- as.numeric(fitj$ec_core_cutoff); MIN_MARKERS <- as.integer(fitj$min_markers_detected)
VETO_Q <- as.numeric(fitj$neighbourhood$q_threshold)

# ---- 1. score and classify every test-study cell; keep each study's graphs ----------------------
per_cell <- list(); fam_data <- list()
for (fam in names(FAMILIES)) {
  sce <- read_h5ad_as_sce(file.path(a[["h5ad_dir"]], paste0(fam, ".h5ad")))
  X <- SummarizedExperiment::assay(sce, "X")
  cd <- as.data.frame(SingleCellExperiment::colData(sce), stringsAsFactors = FALSE)
  rd <- as.data.frame(SummarizedExperiment::rowData(sce), stringsAsFactors = FALSE)
  stopifnot(inherits(X, "dgCMatrix"), sum(X@x != round(X@x)) == 0)   # raw counts
  ftype <- intersect(c("feature_type", "original_feature_type"), colnames(rd))
  keep <- if (length(ftype)) as.character(rd[[ftype[1]]]) == "Gene Expression" else rep(TRUE, nrow(X))
  X <- X[keep, , drop = FALSE]; rd <- rd[keep, , drop = FALSE]
  symbols <- if ("gene_symbol" %in% colnames(rd)) as.character(rd$gene_symbol) else rownames(X)
  total <- Matrix::colSums(X)
  stopifnot(all(total > 0))
  gene_rows <- match(STANDARDISED_GENES, symbols)
  stopifnot(!anyNA(gene_rows), !any(STANDARDISED_GENES %in% symbols[duplicated(symbols)]))
  counts6 <- as.matrix(X[gene_rows, , drop = FALSE]); rownames(counts6) <- STANDARDISED_GENES
  detected <- as.integer(Matrix::colSums(X[gene_rows[seq_along(EC_CORE)], , drop = FALSE] > 0))
  arm <- as.character(cd$object_dataset_id)
  want <- EXPECTED_ARMS[[fam]]; got <- table(arm)
  stopifnot(identical(sort(names(got)), sort(names(want))), all(as.integer(got[names(want)]) == as.integer(want)))

  L6 <- cpm(counts6, total)
  z <- (L6 - LIVE_MU[rownames(L6)]) / LIVE_SD[rownames(L6)]
  score <- colMeans(z[EC_CORE, , drop = FALSE])
  cand <- (score >= LIVE_CUT) & (detected >= MIN_MARKERS)
  full_lcpm <- cpm(X, total)
  graphs <- ec_graphs_build(full_lcpm, arm, excluded_rows = integer(0))
  lcpm_core <- as.matrix(full_lcpm[match(EC_CORE, symbols), , drop = FALSE]); rownames(lcpm_core) <- EC_CORE
  rm(full_lcpm); invisible(gc(verbose = FALSE))
  stopifnot(!any(vapply(graphs, is.null, logical(1))))
  cons <- ec_graph_veto(cand, graphs, VETO_Q)
  presentation <- rep(CLASS_NONE, ncol(X)); presentation[cand] <- CLASS_LOW; presentation[cons] <- CLASS_HIGH
  reason <- rep("did not meet the EC_CORE score cut-off", ncol(X))
  reason[(score >= LIVE_CUT) & !(detected >= MIN_MARKERS)] <- "score met, fewer than 2 of 4 markers detected"
  reason[cand & !cons] <- "score and floor met, not locally consistent"
  reason[cons] <- "score, floor and local consistency all met"
  per_cell[[fam]] <- data.frame(family = fam, family_label = unname(FAMILIES[[fam]]), cell_key = colnames(X), arm_id = arm,
    independent_unit_id = as.character(cd$independent_unit_id), donor_id = as.character(cd$donor_id),
    author_label_class = as.character(cd$author_label_class), gex_total_counts = as.numeric(total),
    markers_detected = detected, ec_core_score = unname(score), score_pass = unname(score >= LIVE_CUT),
    floor_pass = detected >= MIN_MARKERS, candidate = unname(cand), locally_consistent = unname(cons),
    presentation_class = presentation, class_reason = reason, cpm_order = "divide_then_multiply",
    stringsAsFactors = FALSE)
  fam_data[[fam]] <- list(cell_key = colnames(X), arm = arm, graphs = graphs, lcpm_core = lcpm_core)
  rm(sce, X); invisible(gc(verbose = FALSE))
}
write.table(do.call(rbind, per_cell), gzfile(a[["scoring"]]), sep = "\t", row.names = FALSE, quote = FALSE)
cells <- read.delim(gzfile(a[["scoring"]]), stringsAsFactors = FALSE)   # downstream steps read the table as written
stopifnot(nrow(cells) == 5214L)

# ---- 2. cutoff sweep and the transported leave-one-group-out refits -------------------------
stab <- read.delim(a[["stability"]], stringsAsFactors = FALSE)
SWEEP <- sort(unique(c(seq(-0.90, 1.50, by = 0.02), LIVE_CUT)))
loo_rows <- stab[stab$design == "Leave one arm out", ]

dg <- read.delim(gzfile(a[["meta"]]), stringsAsFactors = FALSE, colClasses = "character")
dg <- dg[dg$model_context == "primary_intestine", ]
ROSTER <- jsonlite::fromJSON(a[["roster"]])                         # the seven Crohn's donors are left out
crohn <- dg$study_family == ROSTER$family_id & dg$biological_unit %in% ROSTER$excluded_units
stopifnot(sum(crohn) == 35L)
dg <- dg[!crohn, ]
stopifnot(nrow(dg) == 1164L)
d_counts <- vapply(paste0("raw_count__", EC_CORE), function(cc) as.numeric(dg[[cc]]), numeric(nrow(dg)))
colnames(d_counts) <- EC_CORE
dL <- log1p(d_counts * (1e6 / as.numeric(dg$cell_total_counts)))
d_arm <- dg$arm_id
pop_sd <- function(m, mu) sqrt(colMeans(sweep(m, 2, mu, "-")^2))
fit_fold <- function(fit) {
  mu <- colMeans(dL[fit, , drop = FALSE]); sdv <- pop_sd(dL[fit, , drop = FALSE], mu)
  s <- rowMeans(sweep(sweep(dL, 2, mu, "-"), 2, sdv, "/"))
  list(mu = mu, sd = sdv, threshold = ec_gmm_fit(s[fit])$threshold, fit_cells = sum(fit))
}
FOLDS <- lapply(loo_rows$held_out, function(g) fit_fold(d_arm != g))
names(FOLDS) <- loo_rows$held_out
stopifnot(max(abs(vapply(FOLDS, `[[`, 0, "threshold") - loo_rows$threshold)) < 1e-9)

sweep_rows <- list(); tr_rows <- list(); tr_tot <- list()
for (fam in SWEEP_FAMILIES) {
  fd <- fam_data[[fam]]; graphs <- fd$graphs; lcpm_core <- fd$lcpm_core
  ff <- cells[cells$family == fam, , drop = FALSE]
  ff <- ff[match(fd$cell_key, ff$cell_key), , drop = FALSE]
  for (g in names(FOLDS)) {
    fo <- FOLDS[[g]]
    s_fold <- colMeans((lcpm_core - fo$mu[EC_CORE]) / fo$sd[EC_CORE])
    cand_f <- (s_fold >= fo$threshold) & (ff$markers_detected >= MIN_MARKERS)
    cons_f <- if (any(cand_f)) ec_graph_veto(cand_f, graphs, VETO_Q) else cand_f
    tr_tot[[length(tr_tot) + 1L]] <- data.frame(family = fam, held_out = g, candidates = sum(cand_f), calls = sum(cons_f))
    for (k in names(CLASS_KEEP)) {
      i <- ff$author_label_class == k
      if (!any(i)) next
      tr_rows[[length(tr_rows) + 1L]] <- data.frame(family = fam, study = unname(STUDY_LABEL[[fam]]), held_out = g,
        raw_threshold = fo$threshold, class = unname(CLASS_KEEP[[k]]), n = sum(i), pre = sum(cand_f & i),
        post = sum(cons_f & i), stringsAsFactors = FALSE)
    }
  }
  cls <- ff$author_label_class
  for (cut in SWEEP) {
    cand <- (ff$ec_core_score >= cut) & (ff$markers_detected >= MIN_MARKERS)
    cons <- if (any(cand)) ec_graph_veto(cand, graphs, VETO_Q) else cand
    for (k in names(CLASS_KEEP)) {
      i <- cls == k
      if (!any(i)) next
      sweep_rows[[length(sweep_rows) + 1L]] <- data.frame(family = fam, study = unname(STUDY_LABEL[[fam]]), cutoff = cut,
        class = unname(CLASS_KEEP[[k]]), n = sum(i), pre = sum(cand & i), post = sum(cons & i), stringsAsFactors = FALSE)
    }
  }
}
sw <- do.call(rbind, sweep_rows)
sw$pre_rate <- 100 * sw$pre / sw$n; sw$post_rate <- 100 * sw$post / sw$n
tr <- do.call(rbind, tr_rows); trt <- do.call(rbind, tr_tot)
tr$post_rate <- 100 * tr$post / tr$n

# display position of a refit: the highest cutoff on the full fit's scale giving the same candidate count, the three
# test studies pooled
ok <- cells$markers_detected >= MIN_MARKERS
s_ok <- sort(cells$ec_core_score[ok], decreasing = TRUE)
equiv_cut <- function(n) {
  c_n <- s_ok[n]
  if (sum(cells$ec_core_score >= c_n & ok) != n) stop("score ties prevent an exact equivalent cut-off at ", n)
  c_n
}
fold_pos <- do.call(rbind, lapply(split(trt, trt$held_out), function(d) {
  n <- sum(d$candidates)
  data.frame(held_out = d$held_out[1], candidates = n, calls = sum(d$calls), equiv_cutoff = equiv_cut(n), stringsAsFactors = FALSE)
}))
fold_pos$raw_threshold <- loo_rows$threshold[match(fold_pos$held_out, loo_rows$held_out)]
fold_pos <- fold_pos[order(fold_pos$equiv_cutoff), ]

# the full-fit classifier evaluated exactly at every refit position
exact_rows <- list()
for (fam in SWEEP_FAMILIES) {
  fd <- fam_data[[fam]]
  ff <- cells[cells$family == fam, , drop = FALSE]; ff <- ff[match(fd$cell_key, ff$cell_key), , drop = FALSE]
  for (g in names(FOLDS)) {
    cut <- fold_pos$equiv_cutoff[match(g, fold_pos$held_out)]
    cand <- ff$ec_core_score >= cut & ff$markers_detected >= MIN_MARKERS
    cons <- if (any(cand)) ec_graph_veto(cand, fd$graphs, VETO_Q) else cand
    for (k in names(CLASS_KEEP)) {
      i <- ff$author_label_class == k
      if (!any(i)) next
      exact_rows[[length(exact_rows) + 1L]] <- data.frame(family = fam, study = unname(STUDY_LABEL[[fam]]), cutoff = cut,
        class = unname(CLASS_KEEP[[k]]), n = sum(i), pre = sum(cand & i), post = sum(cons & i), stringsAsFactors = FALSE)
    }
  }
}
exact <- unique(do.call(rbind, exact_rows))
exact$pre_rate <- 100 * exact$pre / exact$n; exact$post_rate <- 100 * exact$post / exact$n
sw <- unique(rbind(sw, exact))
SWEEP_ORDER <- CLASS_KEEP[c("EXPLICIT_EC", "BROAD_EEC_PARENT_UNRESOLVED", "TERMINAL_NON_EC_EEC")]   # published row order
sw <- sw[order(sw$family, match(sw$class, SWEEP_ORDER), sw$cutoff), ]
stopifnot(!anyDuplicated(sw[, c("family", "class", "cutoff")]))
tr$curve_rate <- mapply(function(fam, k, cut) sw$post_rate[sw$family == fam & sw$class == k & sw$cutoff == cut],
                        tr$family, tr$class, fold_pos$equiv_cutoff[match(tr$held_out, fold_pos$held_out)])
wt <- function(x, path) write.table(x, if (grepl("[.]gz$", path)) gzfile(path) else path, sep = "\t", row.names = FALSE, quote = FALSE)
wt(merge(tr, fold_pos[, c("held_out", "equiv_cutoff", "candidates")], by = "held_out"), a[["refits"]])
wt(sw, a[["sweep"]])
print(fold_pos)
