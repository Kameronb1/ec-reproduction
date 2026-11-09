# Figure S5C (TPH1-only classifier): which cells the TPH1-only classifier classifies as EC, beside the four-marker
# and 20-gene classifiers.
#
# The TPH1-only classifier is the four-marker classifier refitted with TPH1 alone (figure2/four_marker_fit.R
# markers=TPH1 min_detected=1: the same 1,199 fitting cells, standardization, two-component mixture and neighborhood
# check; a cell needs at least one TPH1 read). It is scored on the test studies by figure2/heldout_marker_calls.R and
# on the McCauley 2023 HIO grafts by figure2/four_marker_fit.R, as the four-marker classifier is.
#
# What it computes
#   1. Per donor, sample, day (Singh 2024) or HIO sample group: EC cells (primary study label) classified as EC by
#      the TPH1-only classifier, on the units of calls_by_donor.tsv and hio_calls_by_organoid.tsv.
#   2. Per study: other terminal EEC cells (primary study label) classified as EC by each of the three classifiers.
#   Check: the four-marker and 20-gene counts recomputed here on the same units equal the per-unit tables drawn in S5C.
#
# Inputs : cells4, cells1 (test-study per-cell classes, four-marker and TPH1-only; figure2/heldout_marker_calls.R)
#          hio4, hio1 (HIO per-cell classes, four-marker and TPH1-only; figure2/four_marker_fit.R)
#          predictions, hio_predictions (20-gene classifier per cell; classifier.R)
#          by_donor, by_organoid (calls_by_donor.tsv, hio_calls_by_organoid.tsv)
# Outputs: units (study, unit, n_rule, called_tph1), other (study, n, called_rule, called_tph1, called_cls)
#   Rscript tph1_calls.R cells4=... cells1=... hio4=... hio1=... predictions=... hio_predictions=... by_donor=...
#           by_organoid=... units=... other=...
args <- commandArgs(trailingOnly = TRUE)
a <- setNames(sub("^[^=]*=", "", args), sub("=.*$", "", args))
rd <- function(p) read.delim(if (grepl("\\.gz$", p)) gzfile(p) else p, stringsAsFactors = FALSE, encoding = "UTF-8")

FAM <- c(gut_cell_atlas = "SF_ELMENTAITE2021_GUT_ATLAS", hickey = "SF_HICKEY2023_ADULT_EEC",
         egozi = "SF_EGOZI2021_DEVELOPMENTAL_EEC", gse238274 = "SF_INTESTINAL_GSE238274")
SINGH <- "SF_INTESTINAL_GSE238274"; HIO <- "SF_GSE214852_HIO_GRAFT"
day <- function(id) sub("^GSE238274__(D[0-9]+):.*$", "\\1", id)

# test studies: one row per cell with its unit (as per_sample.R) and whether each of the three classifiers classifies it as EC
c4 <- rd(a[["cells4"]]); c1 <- rd(a[["cells1"]])
stopifnot(identical(c4$project_cell_id, c1$project_cell_id))
pr <- rd(a[["predictions"]]); pr <- pr[pr$panel_size == 20, ]
k <- match(c4$project_cell_id, pr$project_cell_id); stopifnot(!anyNA(k))
tis <- data.frame(study = unname(FAM[c4$family]), cls = c4$author_label_class,
                  unit = sub(" \\([^)]*\\)$", "", pr$biological_unit_id[k]),
                  rule = as.logical(c4$called), tph1 = as.logical(c1$called), cls20 = as.logical(pr$call[k]),
                  stringsAsFactors = FALSE)
tis$unit[tis$study == SINGH] <- day(c4$project_cell_id[tis$study == SINGH])

# HIO grafts: cells with a primary study label (EC or other EEC)
h4 <- rd(a[["hio4"]]); h1 <- rd(a[["hio1"]])
hp <- rd(a[["hio_predictions"]]); hp <- hp[hp$panel_size == 20, ]
j4 <- match(hp$project_cell_id, h4$cell_key); j1 <- match(hp$project_cell_id, h1$cell_key); stopifnot(!anyNA(j4), !anyNA(j1))
hio <- data.frame(study = HIO, cls = ifelse(hp$primary_label == "POSITIVE", "EXPLICIT_EC", "TERMINAL_NON_EC_EEC"),
                  unit = hp$biological_unit_id, rule = h4$presentation_class[j4] == "HIGH_CONFIDENCE_EC",
                  tph1 = h1$presentation_class[j1] == "HIGH_CONFIDENCE_EC", cls20 = as.logical(hp$call),
                  stringsAsFactors = FALSE)
cells <- rbind(tis, hio)

# 1. EC cells per unit
ec <- cells[cells$cls == "EXPLICIT_EC", ]
u <- aggregate(cbind(n = 1, rule, tph1, cls20) ~ study + unit, data = ec, FUN = sum)
pub <- rbind(rd(a[["by_donor"]]), rd(a[["by_organoid"]]))
m <- match(paste(pub$study, pub$unit), paste(u$study, u$unit)); stopifnot(!anyNA(m), nrow(u) == nrow(pub))
stopifnot(all(u$n[m] == pub$n_rule), all(u$rule[m] == pub$called_rule), all(u$cls20[m] == pub$called_cls))  # check
units <- data.frame(study = pub$study, unit = pub$unit, n_rule = pub$n_rule, called_tph1 = u$tph1[m])

# 2. other terminal EEC cells per study
o <- aggregate(cbind(n = 1, called_rule = rule, called_tph1 = tph1, called_cls = cls20) ~ study,
               data = cells[cells$cls == "TERMINAL_NON_EC_EEC", ], FUN = sum)
o <- o[order(match(o$study, c(unname(FAM), HIO))), ]
stopifnot(nrow(o) == 5L)

write.table(units, a[["units"]], sep = "\t", quote = TRUE, row.names = FALSE)
write.table(o, a[["other"]], sep = "\t", quote = TRUE, row.names = FALSE)
