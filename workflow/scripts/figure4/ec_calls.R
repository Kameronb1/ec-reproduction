# Figure 4A (each row's start and end samples) and 4B (cells classified as EC by the label-free four-marker EC
# classifier, start and end of each row).
#
# Rows compare SC-islets before and after a change: in vitro -> graft (eight rows) or earlier -> later in vitro culture
# (three rows). The Balboa 2022 3-month grafts (row B-0-M3) are pooled: their three samples count as one sample.
#
# Each SC-SIC and SC-beta cell is scored with the four-marker classifier's constants (from the classifier unit): per
# marker, z = (log(1 + count / library size * 1e6) - mean) / sd, the score is the mean of the four z, and a cell is
# classified as EC when its score is at or above the cutoff and at least two of the four markers are detected. Graft rows
# use SoupX-corrected counts at both ends; the culture rows use uncorrected counts (Balboa 2022 week 6 from its GSE167880
# arm, Augsornworawat 2023 months 6 and 12 from the long-culture arm).
#
# Outputs:
#   rows      row_samples.tsv: the samples at each row end (a pooled end is its samples joined with "+")
#   cells     C_cells.tsv.gz: one line per SC-SIC and SC-beta cell at each row end, with its four-marker score, a
#             weight 1 / (cells of its sample) / (samples at its end), and its call (score at or above the cutoff and at
#             least two markers detected)
#   shares    four_marker_row_shares.tsv: per row end and cell type, the weighted share of cells with a score at or above
#             the cutoff (pct; the mean of the samples' shares) and the share classified as EC (pct_two_marker_rule).
#             Computed from C_cells.tsv.gz as written.
#
# Inputs: markers_soupx, markers_uncorrected = per-cell marker counts and library size; constants = fitted_constants.json
#         of the four-marker classifier; comparisons, comparisons_nuclei = row definitions
#
#   Rscript ec_calls.R markers_soupx=... markers_uncorrected=... constants=... comparisons=... comparisons_nuclei=... \
#           rows=... cells=... shares=...
args <- commandArgs(trailingOnly = TRUE)
A <- setNames(sub("^[^=]*=", "", args), sub("=.*$", "", args))
suppressPackageStartupMessages(library(jsonlite))
rt <- function(p, ...) read.delim(if (grepl("\\.gz$", p)) gzfile(p) else p, check.names = FALSE, stringsAsFactors = FALSE, ...)
wt <- function(d, path) { con <- if (grepl("\\.gz$", path)) gzfile(path, "wb") else file(path, "wb"); on.exit(close(con))
  write.table(d, con, sep = "\t", quote = FALSE, row.names = FALSE, na = "NA", eol = "\n") }

# ---- rows ----
ROWS <- c("A-HUES8", "A-WS4", "B-0-M1", "B-0-M3", "B-3-M6-CTRL", "B-3-M6-STZ", "H-10w", "P-TXP", "B-w0-w6", "P-w2-m6", "P-w2-m12")
GRAFT_ROWS <- ROWS[1:8]
CULT_ARM <- c("B-w0-w6" = "GSE167880", "P-w2-m6" = "GSE199636_LONGCULTURE", "P-w2-m12" = "GSE199636_LONGCULTURE")
cmp <- rbind(rt(A[["comparisons"]]), rt(A[["comparisons_nuclei"]]))[, c("comparison", "ref", "test")]
cmp <- cmp[match(ROWS, cmp$comparison), ]; stopifnot(!anyNA(cmp$comparison))
REF <- setNames(strsplit(cmp$ref, ";"), cmp$comparison); TEST <- setNames(strsplit(cmp$test, ";"), cmp$comparison)
M3 <- TEST[["B-0-M3"]]; stopifnot(identical(M3, c("SRR13805903", "SRR13805904", "SRR13805905")))
TEST[["B-0-M3"]] <- paste(M3, collapse = "+")
wt(data.frame(id = ROWS, ref_samples = vapply(REF[ROWS], paste, "", collapse = ";"), test_samples = vapply(TEST[ROWS], paste, "", collapse = ";")), A[["rows"]])

# ---- four-marker score and call of every cell ----
K <- fromJSON(A[["constants"]]); CORE4 <- K$ec_core_genes
stopifnot(identical(CORE4, c("TPH1", "SLC18A1", "FEV", "LMX1A")), K$min_markers_detected == 2)
st <- K$standardisation; st <- st[match(CORE4, st$gene), ]; stopifnot(!anyNA(st$gene))
score_cells <- function(path) {
  x <- rt(path, colClasses = c(cell_id = "character", sample = "character"))
  z <- lapply(seq_along(CORE4), function(j) (log1p(x[[CORE4[j]]] / x$library_size * 1e6) - st$mean[j]) / st$sd[j])
  x$four_marker_score <- (((z[[1]] + z[[2]]) + z[[3]]) + z[[4]]) / 4
  x$four_marker_markers <- rowSums(as.matrix(x[, CORE4]) > 0)
  x$four_marker_candidate <- x$four_marker_score >= K$ec_core_cutoff & x$four_marker_markers >= K$min_markers_detected
  x$sample[x$sample %in% M3] <- TEST[["B-0-M3"]]
  x }
ps <- score_cells(A[["markers_soupx"]]); pu <- score_cells(A[["markers_uncorrected"]])

# ---- per-cell table of each row's start ("from") and end ("to") ----
ENDS <- list(from = REF, to = TEST)
cells_for <- function(id, end, p) {
  s <- ENDS[[end]][[id]]
  x <- if (id %in% GRAFT_ROWS) ps[ps$sample %in% s & ps$population == p, ] else {
    a <- pu[pu$sample %in% s & pu$population == p, ]; a[a$arm == CULT_ARM[[id]], ] }
  stopifnot(setequal(unique(x$sample), s), !anyDuplicated(x$cell_id))
  n <- table(x$sample)
  data.frame(id = id, population = p, end = end, sample = x$sample, four_marker_score = x$four_marker_score,
             weight = 1 / as.numeric(n[x$sample]) / length(s), four_marker_candidate = x$four_marker_candidate, stringsAsFactors = FALSE) }
CC <- do.call(rbind, lapply(ROWS, function(id) do.call(rbind, lapply(c("SC_SIC", "SC_BETA"), function(p)
  do.call(rbind, lapply(c("from", "to"), function(e) cells_for(id, e, p)))))))
wt(CC, A[["cells"]])

# ---- shares per row end, from the table as written ----
CC <- rt(A[["cells"]])
grp <- split(CC, paste(CC$id, CC$population, CC$end))
stopifnot(all(vapply(grp, function(x) abs(sum(x$weight) - 1) < 1e-9, TRUE)))
SH <- do.call(rbind, lapply(grp, function(x) data.frame(id = x$id[1], population = x$population[1], end = x$end[1], cells = nrow(x),
  samples = length(unique(x$sample)), pct = 100 * sum(x$weight * (x$four_marker_score >= K$ec_core_cutoff)),
  pct_two_marker_rule = 100 * sum(x$weight * x$four_marker_candidate), stringsAsFactors = FALSE)))
wt(SH, A[["shares"]])
