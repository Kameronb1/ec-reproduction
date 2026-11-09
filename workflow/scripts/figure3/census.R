# Figure 3A: how many primary-study-labeled cells of each SC-islet population every study contributes.
#
# Sums the per-sample census of primary study labels into one row per study (and per protocol for the
# inDrops study Veres 2019, GSE114412) and five populations: SC-SIC, SC-beta, SC-alpha, SC-delta and every
# other label. Also counts the distinct samples per row. The panel draws the six 10x studies (block
# "main"); Veres 2019, which is not used to select or fit the SC-SIC genes, is kept in the table only.
#
# Inputs:  census=  per-sample census (dataset, sample, unit, protocol, endpoint, target, n_cells, ...)
# Outputs: out=     long table block / part / dataset / protocol / row_id / key / value
#
#   Rscript census.R census=primary_label_population_census.tsv out=panel_a_census.tsv

args <- commandArgs(trailingOnly = TRUE)
arg <- setNames(sub("^[^=]*=", "", args), sub("=.*", "", args))

HELD_OUT <- "GSE114412"
# drawing order of the six 10x studies, then Veres 2019's three protocols
ORDER <- c("GSE151117", "GSE167880", "GSE202497", "GSE226346", "GSE185036", "GSE199636")
# The key text is the published table's wording, so the output matches it byte for byte;
# OTHER means every other primary study label.
CLASSES <- c(SC_SIC = "SC-SIC", SC_BETA = "SC-beta", SC_ALPHA = "SC-alpha", SC_DELTA = "SC-delta",
             OTHER = "other primary study labels (mostly not mature endocrine)")

census <- read.delim(arg[["census"]], stringsAsFactors = FALSE, colClasses = "character")
census$n_cells <- as.integer(census$n_cells)
if (anyNA(census$n_cells)) stop("a census row has no cell count")
if (!setequal(unique(census$dataset), c(ORDER, HELD_OUT))) stop("the census does not hold exactly the seven studies")
if (!all(census$target %in% names(CLASSES))) stop("the census holds a population with no class")

census$protocol <- ifelse(census$dataset == HELD_OUT, census$protocol, "")
census$row_id <- paste(census$dataset, census$protocol, sep = "__")
rows <- c(paste0(ORDER, "__"), paste(HELD_OUT, c("v8", "x1", "x2"), sep = "__"))
if (!setequal(rows, unique(census$row_id))) stop("unexpected study / protocol rows")

part <- function(id) census[census$row_id == id, ]
cells <- do.call(rbind, lapply(rows, function(id) {
  z <- part(id)
  data.frame(block = if (z$dataset[1] == HELD_OUT) "held_out" else "main", part = "cells", dataset = z$dataset[1],
             protocol = z$protocol[1], row_id = id, key = unname(CLASSES),
             value = vapply(names(CLASSES), function(k) sum(z$n_cells[z$target == k]), integer(1)))
}))
samples <- do.call(rbind, lapply(rows, function(id) {
  z <- part(id)
  data.frame(block = if (z$dataset[1] == HELD_OUT) "held_out" else "main", part = "samples", dataset = z$dataset[1],
             protocol = z$protocol[1], row_id = id, key = "distinct_sample_id", value = length(unique(z$sample)))
}))
marks <- rbind(cells, samples)
rownames(marks) <- NULL
write.table(marks, arg[["out"]], sep = "\t", quote = TRUE, row.names = FALSE, na = "")
