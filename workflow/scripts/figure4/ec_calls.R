# Figure 4A (each row's earlier and later samples) and 4C left (per-cell scores), and the share of cells classified as
# EC per row. Rows: before and after the change (in vitro -> graft, or earlier -> later in
# vitro). Three per-cell classifications: the label-free four-marker EC classifier ungated (FOUR_MARKER_UNGATED), gated
# by its own neighborhood veto (FOUR_MARKER_GATED), and the SC-SIC classifier with its local-enrichment gate
# (CLASSIFIER; CLASSIFIER_STAGE_ONE is its ungated first stage).
#
# Per sample and population the percentage of cells classified as EC is computed, then each row end is the mean over
# that end's samples. The Balboa 2022 3-month grafts (row B-0-M3) are pooled: their three samples count as one sample,
# whose percentage is over all their cells (SRR13805905's cells on uncorrected counts: SoupX could not estimate its
# contamination).
# Three tables:
#   rows_called        all eleven rows on uncorrected counts (percentages per sample rounded to 3 decimals): the share
#                      classified as EC at each end, and the samples of each row end (4A)
#   rows_called_soupx  the eight graft rows on SoupX-corrected counts, beside the uncorrected values
#   C_cells            Figure 4C left: one line per SC-SIC and SC-beta cell of each row's earlier ("from") and later
#                      ("to") samples, with its four-marker score, gated four-marker result, SC-SIC classifier score
#                      minus the threshold of its fold, classifier result, and a weight 1 / (cells of its sample) /
#                      (samples at its end), so that weighted shares are the means over samples above. Graft rows come from the
#                      SoupX-corrected table, the culture rows from the uncorrected one (Balboa 2022 week 6 from its
#                      GSE167880 arm, Augsornworawat 2023 months 6 and 12 from the long-culture arm).
#
# Inputs: calls=<per-cell classifier results, uncorrected counts>
#         calls_soupx=<per-cell classifier results, SoupX-corrected counts, graft arms>
#         comparisons=<row definitions, whole-cell studies>  comparisons_nuclei=<row definitions, Augsornworawat 2023>
# Outputs: rows=<rows_called.tsv>  rows_soupx=<rows_called_soupx.tsv>  cells=<C_cells.tsv.gz>
args <- commandArgs(trailingOnly = TRUE)
A <- setNames(sub("^[^=]*=", "", args), sub("=.*$", "", args))
rt <- function(p) read.delim(if (grepl("\\.gz$", p)) gzfile(p) else p, check.names = FALSE, stringsAsFactors = FALSE)
wt <- function(d, path) { con <- file(path, "wb"); on.exit(close(con)); write.table(d, con, sep = "\t", quote = FALSE, row.names = FALSE, na = "NA", eol = "\n") }

ROWS <- c("A-HUES8", "A-WS4", "B-0-M1", "B-0-M3", "B-3-M6-CTRL", "B-3-M6-STZ", "H-10w", "P-TXP", "B-w0-w6", "P-w2-m6", "P-w2-m12")
GRAFT_ROWS <- ROWS[1:8]; SOUPX_ROWS <- GRAFT_ROWS
cmp <- rbind(rt(A[["comparisons"]]), rt(A[["comparisons_nuclei"]]))[, c("comparison", "ref", "test")]
cmp <- cmp[match(ROWS, cmp$comparison), ]; stopifnot(!anyNA(cmp$comparison))
REF <- setNames(strsplit(cmp$ref, ";"), cmp$comparison); TEST <- setNames(strsplit(cmp$test, ";"), cmp$comparison)
# the pooled end: one sample named by its members joined with "+"
M3 <- TEST[["B-0-M3"]]; stopifnot(identical(M3, c("SRR13805903", "SRR13805904", "SRR13805905")))
TEST[["B-0-M3"]] <- paste(M3, collapse = "+")

# ---- uncorrected counts: per sample x population, then rows ----
DRAWN_ARMS <- c("GSE151117", "GSE139535", "GSE167880", "GSE185036", "GSE199636", "GSE199636_LONGCULTURE")
pc_all <- rt(A[["calls"]]); pc_all$sample[pc_all$arm == "GSE167880" & pc_all$sample %in% M3] <- TEST[["B-0-M3"]]
pc <- pc_all[pc_all$arm %in% DRAWN_ARMS, ]
key <- paste(pc$arm, pc$sample, pc$population, sep = "\r"); sp <- split(seq_len(nrow(pc)), key)
pct <- function(v) round(100 * sum(v) / length(v), 3)
bysp <- do.call(rbind, lapply(sp, function(ix) data.frame(arm = pc$arm[ix[1]], sample = pc$sample[ix[1]], population = pc$population[ix[1]],
  fold_left_out = pc$fold_left_out[ix[1]], cells = length(ix),
  four_marker_called_pct = pct(pc$four_marker_candidate[ix]), four_marker_gated_pct = pct(pc$four_marker_gated[ix]),
  stage_one_called_pct = pct(pc$stage_one[ix]), classifier_called_pct = pct(pc$classifier_call[ix]), stringsAsFactors = FALSE)))
ARM_OF <- setNames(bysp$arm, bysp$sample)
RULES <- c(FOUR_MARKER_UNGATED = "four_marker_called_pct", FOUR_MARKER_GATED = "four_marker_gated_pct",
           CLASSIFIER_STAGE_ONE = "stage_one_called_pct", CLASSIFIER = "classifier_called_pct")
armv <- function(col, s, p) { v <- bysp[[col]][match(paste(s, p), paste(bysp$sample, bysp$population))]; stopifnot(!anyNA(v)); mean(v) }
rows <- do.call(rbind, lapply(ROWS, function(id) do.call(rbind, lapply(c("SC_SIC", "SC_BETA"), function(p) do.call(rbind, lapply(names(RULES), function(r)
  data.frame(id = id, population = p, rule = r, from = armv(RULES[[r]], REF[[id]], p), to = armv(RULES[[r]], TEST[[id]], p),
    gated = r %in% c("FOUR_MARKER_GATED", "CLASSIFIER"),
    fold_left_out = unique(bysp$fold_left_out[bysp$sample %in% c(REF[[id]], TEST[[id]])]),
    ref_arm = paste(unique(ARM_OF[REF[[id]]]), collapse = ";"), test_arm = paste(unique(ARM_OF[TEST[[id]]]), collapse = ";"),
    ref_samples = paste(REF[[id]], collapse = ";"), test_samples = paste(TEST[[id]], collapse = ";"),
    ref_cells = sum(bysp$cells[bysp$sample %in% REF[[id]] & bysp$population == p]),
    test_cells = sum(bysp$cells[bysp$sample %in% TEST[[id]] & bysp$population == p]), stringsAsFactors = FALSE)))))))
rows$change <- rows$to - rows$from
wt(rows, A[["rows"]])
rows <- rt(A[["rows"]])   # read back as written (3 decimals): the graft table's uncorrected columns carry these stored values, as the published table does

# ---- SoupX-corrected counts: graft rows (unrounded percentages) ----
ps <- rt(A[["calls_soupx"]]); ps$sample[ps$arm == "GSE167880" & ps$sample %in% M3] <- TEST[["B-0-M3"]]
COLS <- c("four_marker_candidate", "four_marker_gated", "classifier_call")
BS <- aggregate(ps[COLS], ps[c("sample", "population")], function(v) 100 * mean(v))
RULECOL <- c(FOUR_MARKER_UNGATED = "four_marker_candidate", FOUR_MARKER_GATED = "four_marker_gated", CLASSIFIER = "classifier_call")
sxv <- function(col, s, p) { v <- BS[[col]][match(paste(s, p), paste(BS$sample, BS$population))]; stopifnot(!anyNA(v)); mean(v) }
OUT <- do.call(rbind, lapply(SOUPX_ROWS, function(id) do.call(rbind, lapply(c("SC_SIC", "SC_BETA"), function(p) do.call(rbind, lapply(names(RULECOL), function(r) {
  k <- rows$id == id & rows$population == p & rows$rule == r
  data.frame(id = id, population = p, rule = r, uncorrected_from = rows$from[k], uncorrected_to = rows$to[k],
             from = sxv(RULECOL[[r]], REF[[id]], p), to = sxv(RULECOL[[r]], TEST[[id]], p),
             ref_samples = rows$ref_samples[k], test_samples = rows$test_samples[k]) }))))))
OUT$change <- OUT$to - OUT$from; OUT$uncorrected_change <- OUT$uncorrected_to - OUT$uncorrected_from
OUT$direction_same <- sign(round(OUT$change, 6)) == sign(round(OUT$uncorrected_change, 6))
wt(OUT, A[["rows_soupx"]])

# ---- C left: per-cell scores of each row's earlier and later cells ----
THR <- tapply(pc_all$threshold, pc_all$fold_left_out, function(x) { stopifnot(length(unique(x)) == 1L); x[1] })   # classifier cut per fold
CULT_ARM <- c("B-w0-w6" = "GSE167880", "P-w2-m6" = "GSE199636_LONGCULTURE", "P-w2-m12" = "GSE199636_LONGCULTURE")
ENDS <- list(from = REF, to = TEST)
cells_for <- function(id, end, p) {
  s <- ENDS[[end]][[id]]
  x <- if (id %in% SOUPX_ROWS) ps[ps$sample %in% s & ps$population == p, ] else {
    a <- pc_all[pc_all$sample %in% s & pc_all$population == p, ]; a[a$arm == CULT_ARM[[id]], ] }
  stopifnot(setequal(unique(x$sample), s), !anyDuplicated(x$cell_id))
  n <- table(x$sample)
  data.frame(id = id, population = p, end = end, sample = x$sample, four_marker_score = x$four_marker_score,
             four_marker_gated = as.logical(x$four_marker_gated), classifier_margin = x$classifier_score - THR[x$fold_left_out],
             classifier_call = as.logical(x$classifier_call), weight = 1 / as.numeric(n[x$sample]) / length(s), stringsAsFactors = FALSE) }
CC <- do.call(rbind, lapply(ROWS, function(id) do.call(rbind, lapply(c("SC_SIC", "SC_BETA"), function(p)
  do.call(rbind, lapply(c("from", "to"), function(e) cells_for(id, e, p)))))))
# check: the weighted shares of cells classified as EC (gated results) equal the row shares above (graft rows: the
# SoupX table)
called <- function(id, p, r, e) { d <- if (id %in% SOUPX_ROWS) OUT else rows; d[[e]][d$id == id & d$population == p & d$rule == r] }
dev <- unlist(lapply(split(CC, paste(CC$id, CC$population, CC$end)), function(x) c(
  100 * sum(x$weight * x$four_marker_gated) - called(x$id[1], x$population[1], "FOUR_MARKER_GATED", x$end[1]),
  100 * sum(x$weight * x$classifier_call) - called(x$id[1], x$population[1], "CLASSIFIER", x$end[1]))))
stopifnot(length(dev) == 88L, max(abs(dev)) < 1e-2)   # the uncorrected shares are stored at 3 decimals
con <- gzfile(A[["cells"]], "wb"); write.table(CC, con, sep = "\t", quote = FALSE, row.names = FALSE, eol = "\n"); close(con)
