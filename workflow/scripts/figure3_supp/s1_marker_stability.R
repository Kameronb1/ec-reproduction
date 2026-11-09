# Figure S15B-D: how stable the 20 SC-SIC genes of S15A are (the checks Figure 2's supplement runs on its 20 genes).
#   B  each SC-islet study left out of the ranking in turn: rank of every gene in each of the six top-20 lists
#   C  (label swaps) is computed by s1_label_swaps.R from the sample-by-class means this script saves (swaps_input).
#   D  against random genes of matched expression: each study left out in turn, the 20 genes and 200 random 20-gene sets
#      (matched on training-study mean, detection and variance bins) are weighted on the other five studies (positive
#      median study effect, normalized to sum one; equal weights if none is positive) and scored on standardized
#      expression; ROC AUC of SC-SIC vs other endocrine cells in the left-out study.
# Same corpus and ranking as s1_marker_limma.R. Seed 20261010.
# Inputs:  lib (s1_lib.R); the SC endocrine count package: meta, genes, indptr, indices, data
# Outputs: a (a_fold_ranks.tsv.gz, panel B), c (c_matched_random.tsv.gz, panel D), swaps_input (.rds for s1_label_swaps.R)
a <- commandArgs(TRUE); A <- setNames(sub("^[^=]*=", "", a), sub("=.*", "", a))
source(A[["lib"]])
EXCL <- EXCLUDED_MARKERS
TOP_N <- 20L; SEED <- 20261010L; REPS <- 200L
log_step <- function(...) cat(format(Sys.time(), "%H:%M:%S"), ..., "\n")

sc <- read_sc_endocrine(A[["meta"]], A[["genes"]], A[["indptr"]], A[["indices"]], A[["data"]])
expr <- sc$expr; genes <- sc$genes; meta <- sc$meta; nC <- nrow(meta)
pos <- meta$target == "SC_SIC"; neg <- meta$target %in% c("SC_BETA", "SC_ALPHA", "SC_DELTA", "SC_EPSILON", "SC_GAMMA")
stopifnot(all(pos | neg))
labels <- data.frame(project_cell_id = meta$project_cell_id, label = ifelse(pos, "POSITIVE", "NEGATIVE"), stringsAsFactors = FALSE)
STUDIES <- sort(unique(meta$family_id))
stopifnot(length(STUDIES) == 6L)
all_cells <- rep(TRUE, nC)

# the 20 genes of S15A (all six studies)
full <- rank_genes(expr, genes, meta, labels, all_cells)
D_GENES <- head(full$ranking$gene_symbol, TOP_N)

# ---- B: leave-one-study-out rankings --------------------------------------------------------------------------------------
log_step("folds")
ELIG <- list(); FOLD_RANK <- list()
for (s in STUDIES) {
  tr <- meta$family_id != s
  FOLD_RANK[[s]] <- rank_genes(expr, genes, meta, labels, tr)$ranking$gene_symbol
  ELIG[[s]] <- genes[training_feature_mask(expr, genes, meta, tr, PARAMETERS$training_detection_floor,
                                           PARAMETERS$minimum_training_families, EXCL)]
}
top_sets <- lapply(FOLD_RANK, head, TOP_N)
OBS <- length(Reduce(intersect, top_sets))
show <- unique(c(D_GENES, unlist(top_sets)))
A_TAB <- do.call(rbind, lapply(STUDIES, function(s) data.frame(left_out = s, gene = show, rank = match(show, FOLD_RANK[[s]]))))
A_TAB$in_d <- A_TAB$gene %in% D_GENES

# ---- C (computed by s1_label_swaps.R): what the swaps need ------------------------------------------------------------
UG <- sort(unique(unlist(ELIG)))
units <- cp_unit_class_means(expr, meta, labels, all_cells, UG, PARAMETERS$minimum_cells_per_unit_class)
saveRDS(list(units = units, ELIG = ELIG, fold_top30 = lapply(FOLD_RANK, head, 30L), studies = STUDIES, observed = OBS), A[["swaps_input"]])

# ---- D: the 20 genes against matched random sets, each study left out -----------------------------------------------------
log_step("matched random")
eligible <- genes[training_feature_mask(expr, genes, meta, all_cells, PARAMETERS$training_detection_floor,
                                        PARAMETERS$minimum_training_families, EXCL)]
stopifnot(all(D_GENES %in% eligible))
response <- function(te) as.integer(labels$label[te] == "POSITIVE")
# weights and scores as Figure 2's panel fit: positive median training-study effects normalized to sum one
weights_of <- function(effects, null_policy) {
  w <- pmax(effects, 0); w[!is.finite(w)] <- 0
  fallback <- !any(w > 0)
  if (fallback) { if (!null_policy) stop("no positive weight"); w <- rep(1, length(effects)) }
  list(w = w / sum(w), fallback = fallback)
}
auc_real <- function(panel, s) {   # the real panel: standardized on the left-out study's cells only, as the full fit scores them
  tr <- meta$family_id != s; te <- !tr
  prep <- cp_prepare_panel_fit(expr, meta, labels, tr, panel)
  w <- weights_of(prep$family_effects, FALSE)$w
  binary_auc(response(te), as.numeric(Matrix::crossprod(w, cp_apply_standardizer(expr[, te, drop = FALSE], prep$standardizer))))
}
auc_prepared <- function(prep, panel, te, null_policy) {
  ix <- match(panel, prep$panel_genes); stopifnot(!anyNA(ix))
  w <- weights_of(prep$family_effects[ix], null_policy)
  c(auc = binary_auc(response(te), as.numeric(Matrix::crossprod(w$w, prep$standardized[ix, te, drop = FALSE]))), failed = 0,
    equal_weights = as.numeric(w$fallback))
}
C_TAB <- do.call(rbind, lapply(seq_along(STUDIES), function(j) { s <- STUDIES[j]
  tr <- meta$family_id != s; te <- which(!tr)
  draws <- cp_draw_matched_random_panels(expr, genes, meta, tr, D_GENES, eligible, REPS, seed = SEED + j * 10000L, forbidden_genes = EXCL)
  real <- auc_real(D_GENES, s)
  rnd <- do.call(rbind, lapply(split(seq_len(REPS), ceiling(seq_len(REPS) / 25)), function(rr) {
    sets <- lapply(rr, function(r) draws$panels$gene_symbol[draws$panels$replicate == r])
    prep <- cp_prepare_panel_fit(expr, meta, labels, tr, sort(unique(unlist(sets))))
    out <- t(vapply(sets, function(gs) auc_prepared(prep, gs, te, TRUE), c(auc = 0, failed = 0, equal_weights = 0))); rm(prep); invisible(gc()); out }))
  log_step(sprintf("  %s: real AUC %.3f; random median %.3f", s, real, median(rnd[, "auc"], na.rm = TRUE)))
  rbind(data.frame(study = s, replicate = 0L, panel = "REAL", roc_auc = real, fit_failed = 0, equal_weights = 0),
        data.frame(study = s, replicate = seq_len(REPS), panel = "MATCHED_RANDOM", roc_auc = rnd[, "auc"], fit_failed = rnd[, "failed"],
                   equal_weights = rnd[, "equal_weights"])) }))

wt <- function(d, f) { h <- gzfile(f, "wt"); write.table(d, h, sep = "\t", quote = FALSE, row.names = FALSE, na = "NA"); close(h) }
wt(A_TAB, A[["a"]]); wt(C_TAB, A[["c"]])
