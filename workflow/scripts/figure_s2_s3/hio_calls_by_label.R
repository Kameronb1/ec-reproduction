# Figure S2B, McCauley 2023 HIO grafts (GSE214852): cells classified as EC by the four-marker classifier within
# each published label, per HIO sample group.
#
# Published labels: the authors' Endocrine_1 cluster is taken as EC and Endocrine_2 as other terminal EEC.
# Units are the two HIO sample groups wt1 (libraries wt1.1 and wt1.2 combined) and wt2; they are groups of
# cells from one cell line, not independent organoids.
#
# Inputs : cells (the GSE214852 endocrine export), hio (four-marker classes of the HIO cells, from four_marker_fit.R)
# Outputs: bc (hio_BC_classifier_calls_by_label.tsv)
#   Rscript hio_calls_by_label.R cells=... hio=... bc=...
args <- commandArgs(trailingOnly = TRUE)
a <- setNames(sub("^[^=]*=", "", args), sub("=.*$", "", args))

cells <- read.delim(a[["cells"]], stringsAsFactors = FALSE)
stopifnot(nrow(cells) == 598L)
cells$sample <- sub("[.][0-9]+$", "", cells$sample)            # wt1.1 + wt1.2 -> wt1
stopifnot(setequal(cells$sample, c("wt1", "wt2")))
hio <- read.delim(a[["hio"]], stringsAsFactors = FALSE)
hio <- hio[match(cells$figure2_cell_id, hio$cell_key), ]
stopifnot(!anyNA(hio$cell_key), all(hio$arm_id == "GSE214852_WT"),
          all(hio$presentation_class %in% c("HIGH_CONFIDENCE_EC", "INSUFFICIENT_EC_EVIDENCE")))
rule_label <- ifelse(hio$presentation_class == "HIGH_CONFIDENCE_EC", "POSITIVE", "NEGATIVE")
stopifnot(all(cells$author_label %in% c("Endocrine_1", "Endocrine_2")))

bc <- aggregate(cbind(n = 1, rule_called = rule_label == "POSITIVE") ~ sample + author_label,
                data = data.frame(sample = cells$sample, author_label = cells$author_label, rule_label), FUN = sum)
write.table(bc, a[["bc"]], sep = "\t", quote = TRUE, row.names = FALSE)
