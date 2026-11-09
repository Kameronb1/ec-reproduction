# Figure S15C: is the stability of the 20 SC-SIC genes beyond chance? The statistic is the number of genes in all six
# leave-one-study-out top-20 lists (s1_marker_stability.R). It is recomputed with SC-SIC and other labels swapped
#   (a) for whole studies: all 2^6 = 64 assignments (p = exceedances / 64, identity included), and
#   (b) within samples: 200 random draws, each paired sample flipped with probability 1/2, seed 20261010
#       (p = (1 + exceedances) / 201).
# A swap exchanges which sample-by-class mean is SC-SIC, which is what swapping the cells' labels does to them; each fold
# is re-ranked with the same limma fit as S15A. A fold's ranking depends only on the labels of the studies it trains on,
# so identical label sets are fitted once and reused (same fit, same result). The fits are independent; workers=N runs
# them on N R processes (results do not depend on N).
# Inputs:  lib (s1_lib.R); swaps_input (.rds written by s1_marker_stability.R: sample-by-class means, eligible genes and
#          real top-30 per fold); workers (optional, default 1)
# Output:  b (b_label_swaps.tsv.gz, drawn in S15C)
a <- commandArgs(TRUE); A <- setNames(sub("^[^=]*=", "", a), sub("=.*", "", a))
source(A[["lib"]])
TOP_N <- 20L; SEED <- 20261010L; REPS <- 200L
WORKERS <- if ("workers" %in% names(A)) as.integer(A[["workers"]]) else 1L
X <- readRDS(A[["swaps_input"]])
units <- X$units; ELIG <- X$ELIG; STUDIES <- X$studies; OBS <- X$observed
um <- units$metadata

fast_rank <- function(s, lab) {
  k <- um$family_id != s
  u <- list(values = units$values[ELIG[[s]], k, drop = FALSE], metadata = transform(um[k, ], label = lab[k]))
  f <- cp_limma_ec_effect(u)
  f$gene_symbol[order(-f$t, f$gene_symbol, method = "radix")]
}
flip <- function(lab, sel) ifelse(sel, ifelse(lab == "POSITIVE", "NEGATIVE", "POSITIVE"), lab)

# the label sets: the primary study labels, the 64 whole-study swaps (expand.grid order), then the 200 within-sample draws
grid <- as.matrix(expand.grid(rep(list(c(FALSE, TRUE)), length(STUDIES)))); colnames(grid) <- STUDIES
LAB_STUDY <- lapply(seq_len(nrow(grid)), function(r) flip(um$label, um$family_id %in% STUDIES[grid[r, ]]))
unit_id <- paste(um$family_id, um$biological_unit_id, sep = "::")
paired <- names(which(tapply(um$label, unit_id, function(v) setequal(v, c("POSITIVE", "NEGATIVE")))))
set.seed(SEED)
LAB_UNIT <- lapply(seq_len(REPS), function(r) { pick <- paired[stats::runif(length(paired)) < 0.5]; flip(um$label, unit_id %in% pick) })
ALL <- c(list(um$label), LAB_STUDY, LAB_UNIT)

# every distinct (fold, training labels) fit
key <- function(s, lab) paste(s, paste(lab[um$family_id != s], collapse = ","))
keys <- character(); jobs <- list()
for (lab in ALL) for (s in STUDIES) { k <- key(s, lab); if (!(k %in% keys)) { keys <- c(keys, k); jobs[[length(jobs) + 1L]] <- list(s = s, lab = lab) } }
cat(length(jobs), "distinct fits\n")
run_job <- function(j) head(fast_rank(jobs[[j]]$s, jobs[[j]]$lab), 30L)
if (WORKERS > 1L) {
  cl <- parallel::makeCluster(WORKERS)
  parallel::clusterExport(cl, c("units", "ELIG", "um", "jobs", "fast_rank", "cp_limma_ec_effect"))
  parallel::clusterEvalQ(cl, suppressPackageStartupMessages(library(limma)))
  TOPS <- parallel::parLapply(cl, seq_along(jobs), run_job); parallel::stopCluster(cl)
} else TOPS <- lapply(seq_along(jobs), run_job)
names(TOPS) <- keys

for (s in STUDIES) stopifnot(identical(TOPS[[key(s, um$label)]], X$fold_top30[[s]]))   # real labels reproduce the folds
stat <- function(lab) length(Reduce(intersect, lapply(STUDIES, function(s) head(TOPS[[key(s, lab)]], TOP_N))))
stopifnot(stat(um$label) == OBS)
NULL_STUDY <- vapply(LAB_STUDY, stat, 0L)
NULL_UNIT <- vapply(LAB_UNIT, stat, 0L)
B_TAB <- rbind(data.frame(scope = "STUDY_BLOCK", draw = seq_len(nrow(grid)), genes_in_all = NULL_STUDY),
               data.frame(scope = "WITHIN_SAMPLE", draw = seq_len(REPS), genes_in_all = NULL_UNIT))
cat(sprintf("observed %d genes in all six lists; study-block p = %.4f; within-sample p = %.4f\n", OBS,
            sum(NULL_STUDY >= OBS) / length(NULL_STUDY), (1 + sum(NULL_UNIT >= OBS)) / (REPS + 1)))
h <- gzfile(A[["b"]], "wt"); write.table(B_TAB, h, sep = "\t", quote = FALSE, row.names = FALSE, na = "NA"); close(h)
