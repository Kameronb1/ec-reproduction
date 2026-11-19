# Figures S2 and S3: the four-marker EC classifier fitted on the fitting cells, then applied unchanged to
# the two HIO groups. Writes the fitted constants that every other use of the classifier reads.
#
# What it computes
#   1. Keeps the primary-tissue fitting cells (model_context == "primary_intestine") and leaves out the
#      cells of the seven pediatric Crohn's disease donors of Elmentaite 2020 (the donor list in roster):
#      1,164 cells.
#   2. Standardizes TPH1, SLC18A1, FEV, LMX1A (and DDC, TAC1 as diagnostics) by their mean and
#      population SD over those cells; the four-marker score is the mean z-score of the four markers.
#   3. Fits a two-component equal-variance Gaussian mixture to the score; the cutoff is the
#      equal-posterior boundary.
#   4. Candidate = score >= cutoff and >= 2 of 4 markers detected. A candidate is classified as EC when
#      its within-group 25-nearest-neighbor graph holds more candidates than chance (hypergeometric,
#      BH q <= 0.05); otherwise it is above the cutoff but not classified.
#      The graph is run under four gene-exclusion sets for comparison; "none" is the one used.
#   5. Applies the fixed constants and cutoff to the 755 HIO cells (model_context ==
#      "hpsc_derived_hio"); only their within-group graph is built from their own cells.
#
# Inputs : matrix (log1p-CPM, float64, cells x features, C order), features, cells, meta, roster
# Outputs: fit_out (fitted_constants.json), cells_out (class of every fitting cell), counts_out (cells per class),
#          hio_cells_out (class of every HIO cell)
#   Rscript marker_fit.R matrix=... features=... cells=... meta=... roster=... fit_out=... cells_out=... counts_out=... hio_cells_out=...
args <- commandArgs(trailingOnly = TRUE)
a <- setNames(sub("^[^=]*=", "", args), sub("=.*$", "", args))
suppressPackageStartupMessages(library(jsonlite))
source(file.path(dirname(sub("^--file=", "", grep("^--file=", commandArgs(), value = TRUE))), "ec_classifier_lib.R"))

EC_CORE <- c("TPH1", "SLC18A1", "FEV", "LMX1A")
EC_SUPPORT <- c("DDC", "TAC1")                     # diagnostic only
STANDARDISED_GENES <- c(EC_CORE, EC_SUPPORT)
# 27 enteroendocrine and islet hormone and marker genes: the widest of the four gene-exclusion sets compared in step 4
SHIPPED_EXCLUSIONS <- c("ARX", "CCK", "CHGA", "CHGB", "DDC", "FEV", "GCG", "GHRL", "GIP", "HHEX",
  "IAPP", "INS", "LMX1A", "MLN", "NKX6-1", "NTS", "PCSK1", "PCSK2", "PDX1",
  "PPY", "PYY", "SCG5", "SCT", "SLC18A1", "SST", "TAC1", "TPH1")
EXCLUSION_SETS <- list(none = character(0), core4 = EC_CORE, core6 = c(EC_CORE, EC_SUPPORT),
                       shipped27 = SHIPPED_EXCLUSIONS)
MIN_MARKERS_DETECTED <- 2L
VETO_Q <- 0.05
CLASS_HIGH <- "HIGH_CONFIDENCE_EC"; CLASS_LOW <- "LOW_CONFIDENCE_EC"; CLASS_NONE <- "INSUFFICIENT_EC_EVIDENCE"

# ---- load ----------------------------------------------------------------------------------
features <- read.delim(a[["features"]], stringsAsFactors = FALSE)
cells_all <- read.delim(a[["cells"]], stringsAsFactors = FALSE)
n_all <- nrow(cells_all); n_feat <- nrow(features)
con <- file(a[["matrix"]], "rb")
X_all <- readBin(con, "double", n = as.numeric(n_all) * n_feat, size = 8, endian = "little")
close(con)
stopifnot(length(X_all) == as.numeric(n_all) * n_feat)
dim(X_all) <- c(n_feat, n_all)                    # features x cells

meta_all <- read.delim(gzfile(a[["meta"]]), stringsAsFactors = FALSE)
meta_all <- meta_all[!duplicated(meta_all$capture_qualified_key), ]
meta_all <- meta_all[match(cells_all$actual_cell_key, meta_all$capture_qualified_key), , drop = FALSE]
stopifnot(!anyNA(meta_all$arm_id))
# the seven pediatric Crohn's disease donors of Elmentaite 2020 (35 cells) are left out before anything is fitted
roster <- fromJSON(a[["roster"]])
crohn <- meta_all$study_family == roster$family_id & meta_all$biological_unit %in% roster$excluded_units
stopifnot(sum(crohn) == 35L, all(meta_all$model_context[crohn] == "primary_intestine"))
X_all <- X_all[, !crohn, drop = FALSE]; meta_all <- meta_all[!crohn, , drop = FALSE]
cells_all <- cells_all[!crohn, , drop = FALSE]; n_all <- nrow(cells_all)
invisible(gc(verbose = FALSE))

row_of <- function(sym) {
  r <- which(features$feature_symbol == sym)
  if (length(r) != 1L) stop("expected exactly one row for ", sym, ", found ", length(r))
  r
}
gene_rows <- vapply(STANDARDISED_GENES, row_of, integer(1))

# ---- 1. fitting cells ------------------------------------------------------------------------
keep <- meta_all$model_context == "primary_intestine"
X <- X_all[, keep, drop = FALSE]
cells <- cells_all[keep, , drop = FALSE]
arm <- meta_all$arm_id[keep]; biological_unit <- meta_all$biological_unit[keep]
gex_total <- as.numeric(meta_all$cell_total_counts[keep])
n_cells <- sum(keep)
stopifnot(n_cells == 1164L, all(gex_total > 0))

consts <- t(vapply(gene_rows, function(row) { v <- X[row, ]; c(mean = mean(v), sd = sqrt(mean((v - mean(v))^2))) }, numeric(2)))
rownames(consts) <- STANDARDISED_GENES
Z <- vapply(STANDARDISED_GENES, function(s) (X[gene_rows[[s]], ] - consts[s, "mean"]) / consts[s, "sd"], numeric(n_cells))
colnames(Z) <- STANDARDISED_GENES
ec_core_score <- rowMeans(Z[, EC_CORE, drop = FALSE])
detected <- colSums(X[gene_rows[EC_CORE], , drop = FALSE] > 0)
floor_pass <- detected >= MIN_MARKERS_DETECTED

gmm <- ec_gmm_fit(ec_core_score)
cutoff <- gmm$threshold
score_pass <- ec_core_score >= cutoff
candidate <- score_pass & floor_pass

ddc_threshold <- as.numeric(median(Z[, "DDC"]))
tac1_threshold <- as.numeric(median(Z[, "TAC1"]))
serotonergic <- (X[gene_rows[["DDC"]], ] > 0 & Z[, "DDC"] >= ddc_threshold) |
                (X[gene_rows[["TAC1"]], ] > 0 & Z[, "TAC1"] >= tac1_threshold)

# ---- 2. within-group neighborhood check, under each exclusion set ---------------------------
run_neighborhood <- function(excluded_symbols) {
  excluded_rows <- which(features$feature_symbol %in% excluded_symbols)
  local_q <- rep(NA_real_, n_cells); hits_all <- rep(NA_integer_, n_cells)
  for (g in ec_graphs_build(X, arm, excluded_rows)) {
    is_cand_a <- candidate[g$cols]; S <- sum(is_cand_a)
    hits <- rowSums(matrix(is_cand_a[g$nn], nrow = g$n_a, ncol = g$k))
    hits_all[g$cols] <- hits
    if (S >= 2L) {
      foc <- which(is_cand_a)
      p <- phyper(hits[foc] - 1L, S - 1L, (g$n_a - 1L) - (S - 1L), g$k, lower.tail = FALSE)
      local_q[g$cols[foc]] <- p.adjust(p, method = "BH")
    }
  }
  list(q = local_q, hits = hits_all, consistent = candidate & !is.na(local_q) & local_q <= VETO_Q,
       n_excluded_symbols = length(excluded_symbols), n_excluded_rows = length(excluded_rows))
}
runs <- lapply(EXCLUSION_SETS, run_neighborhood)
primary <- runs[["none"]]
locally_consistent <- primary$consistent

n_cand <- sum(candidate)
exclusion_comparison <- do.call(rbind, lapply(names(EXCLUSION_SETS), function(nm) {
  r <- runs[[nm]]
  both_pass <- sum(candidate & primary$consistent & r$consistent)
  both_veto <- sum(candidate & !primary$consistent & !r$consistent)
  only_pri <- sum(candidate & !primary$consistent & r$consistent)
  only_alt <- sum(candidate & primary$consistent & !r$consistent)
  q_both <- !is.na(primary$q) & !is.na(r$q)
  data.frame(exclusion_set = nm, is_primary = identical(nm, "none"),
             excluded_symbols = r$n_excluded_symbols, excluded_matrix_rows = r$n_excluded_rows,
             candidates = n_cand, vetoed = sum(candidate & !r$consistent), high_confidence = sum(r$consistent),
             both_pass = both_pass, both_veto = both_veto,
             vetoed_only_by_primary = only_pri, vetoed_only_by_this = only_alt,
             disagreements = only_pri + only_alt, agreement_pct = 100 * (both_pass + both_veto) / n_cand,
             q_spearman_vs_primary = suppressWarnings(cor(primary$q[q_both], r$q[q_both], method = "spearman")),
             stringsAsFactors = FALSE)
}))

# ---- 3. three classes ------------------------------------------------------------------------
classes_of <- function(consistent) { p <- rep(CLASS_NONE, n_cells); p[candidate] <- CLASS_LOW; p[consistent] <- CLASS_HIGH; p }
presentation <- classes_of(locally_consistent)
reason <- rep("did not meet the EC_CORE score cut-off", n_cells)
reason[score_pass & !floor_pass] <- "score met, fewer than 2 of 4 markers detected"
reason[candidate & !locally_consistent] <- "score and floor met, not locally consistent"
reason[locally_consistent] <- "score, floor and local consistency all met"
class_by_set <- lapply(runs, function(r) classes_of(r$consistent))
exclusion_comparison$class_changed_vs_primary <-
  vapply(class_by_set, function(p) sum(p != presentation), integer(1))[names(EXCLUSION_SETS)]
counts <- table(factor(presentation, levels = c(CLASS_HIGH, CLASS_LOW, CLASS_NONE)))

per_cell <- data.frame(cell_key = cells$actual_cell_key, arm_id = arm, biological_unit = biological_unit,
  gex_total_counts = gex_total, ec_core_score = ec_core_score, markers_detected = detected,
  score_pass = score_pass, floor_pass = floor_pass, candidate = candidate,
  serotonergic_diagnostic = serotonergic, neighbour_candidate_hits = primary$hits, local_q = primary$q,
  locally_consistent = locally_consistent, presentation_class = presentation, class_reason = reason,
  stringsAsFactors = FALSE)
for (nm in setdiff(names(EXCLUSION_SETS), "none")) {
  per_cell[[paste0("local_q__", nm)]] <- runs[[nm]]$q
  per_cell[[paste0("locally_consistent__", nm)]] <- runs[[nm]]$consistent
  per_cell[[paste0("presentation_class__", nm)]] <- class_by_set[[nm]]
}
write.table(per_cell, gzfile(a[["cells_out"]]), sep = "\t", row.names = FALSE, quote = FALSE)
write.table(data.frame(presentation_class = names(counts), n = as.integer(counts), fraction = as.numeric(counts) / n_cells),
            a[["counts_out"]], sep = "\t", row.names = FALSE, quote = FALSE)

fitted <- list(n_cells = n_cells, ec_core_genes = EC_CORE, ec_support_genes = EC_SUPPORT,
  standardisation = lapply(STANDARDISED_GENES, function(g) list(gene = g, mean = consts[g, "mean"], sd = consts[g, "sd"])),
  ec_core_cutoff = cutoff, gaussian_mixture = gmm, min_markers_detected = MIN_MARKERS_DETECTED,
  ddc_threshold = ddc_threshold, tac1_threshold = tac1_threshold,
  neighbourhood = list(neighbours = EC_GRAPH_NEIGHBORS, components = EC_GRAPH_COMPONENTS, top_hvgs = EC_GRAPH_TOP_HVG,
                       detected_fraction_range = EC_GRAPH_DETECT_RANGE, excluded_genes = character(0),
                       excluded_gene_set_name = "none", scope = "WITHIN_ANALYTIC_ARM", q_threshold = VETO_Q),
  classes = c(CLASS_HIGH, CLASS_LOW, CLASS_NONE),
  exclusion_set_comparison = exclusion_comparison)
write(jsonlite::toJSON(fitted, auto_unbox = TRUE, pretty = TRUE, digits = 17), a[["fit_out"]])

# ---- 4. the two HIO groups, scored with the constants as written ------------------------------
fit <- jsonlite::fromJSON(a[["fit_out"]])
sd_tab <- fit$standardisation
hc <- as.matrix(sd_tab[match(STANDARDISED_GENES, sd_tab$gene), c("mean", "sd")]); rownames(hc) <- STANDARDISED_GENES
hcut <- fit$ec_core_cutoff
hk <- meta_all$model_context == "hpsc_derived_hio"
XH <- X_all[, hk, drop = FALSE]; mh <- meta_all[hk, , drop = FALSE]; ch <- cells_all[hk, , drop = FALSE]
rm(X_all, X); invisible(gc(verbose = FALSE))
nh <- sum(hk); harm <- mh$arm_id; hdepth <- as.numeric(mh$cell_total_counts)
stopifnot(nh == 755L, identical(sort(unique(harm)), c("E-MTAB-10268", "GSE214852_WT")))
ZH <- vapply(STANDARDISED_GENES, function(s) (XH[gene_rows[[s]], ] - hc[s, "mean"]) / hc[s, "sd"], numeric(nh))
colnames(ZH) <- STANDARDISED_GENES
h_score <- rowMeans(ZH[, EC_CORE, drop = FALSE])
h_det <- colSums(XH[gene_rows[EC_CORE], , drop = FALSE] > 0)
h_floor <- h_det >= MIN_MARKERS_DETECTED; h_pass <- h_score >= hcut; h_cand <- h_pass & h_floor
h_sero <- (XH[gene_rows[["DDC"]], ] > 0 & ZH[, "DDC"] >= fit$ddc_threshold) |
          (XH[gene_rows[["TAC1"]], ] > 0 & ZH[, "TAC1"] >= fit$tac1_threshold)
h_q <- rep(NA_real_, nh); h_hits <- rep(NA_integer_, nh)
for (g in ec_graphs_build(XH, harm, integer(0))) {
  stopifnot(!is.null(g))
  is_cand_a <- h_cand[g$cols]; S <- sum(is_cand_a)
  hits <- rowSums(matrix(is_cand_a[g$nn], nrow = g$n_a, ncol = g$k)); h_hits[g$cols] <- hits
  if (S >= 2L) {
    foc <- which(is_cand_a)
    p <- phyper(hits[foc] - 1L, S - 1L, (g$n_a - 1L) - (S - 1L), g$k, lower.tail = FALSE)
    h_q[g$cols[foc]] <- p.adjust(p, method = "BH")
  }
}
h_cons <- h_cand & !is.na(h_q) & h_q <= fit$neighbourhood$q_threshold
h_class <- rep(CLASS_NONE, nh); h_class[h_cand] <- CLASS_LOW; h_class[h_cons] <- CLASS_HIGH
h_reason <- rep("did not meet the EC_CORE score cut-off", nh)
h_reason[h_pass & !h_floor] <- "score met, fewer than 2 of 4 markers detected"
h_reason[h_cand & !h_cons] <- "score and floor met, not locally consistent"
h_reason[h_cons] <- "score, floor and local consistency all met"
hp <- data.frame(cell_key = ch$actual_cell_key, arm_id = harm, study_family = mh$study_family,
  model_context = "hpsc_derived_hio", scoring_role = "TRANSPORTED_NOT_IN_FIT", biological_unit = mh$biological_unit,
  gex_total_counts = hdepth, ec_core_score = h_score, markers_detected = h_det, score_pass = h_pass,
  floor_pass = h_floor, candidate = h_cand, serotonergic_diagnostic = h_sero, neighbour_candidate_hits = h_hits,
  local_q = h_q, locally_consistent = h_cons, presentation_class = h_class, class_reason = h_reason,
  stringsAsFactors = FALSE)
stopifnot(length(intersect(hp$cell_key, per_cell$cell_key)) == 0L)
write.table(hp, gzfile(a[["hio_cells_out"]]), sep = "\t", row.names = FALSE, quote = FALSE)
cat("marker_fit:", n_cells, "cells, cutoff", format(cutoff, digits = 17), "; classes", paste(as.integer(counts), collapse = "/"), "\n")
