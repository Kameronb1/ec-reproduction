# Figure S5C (tissue rows): EC cells (by primary study label) classified as EC per donor or sample by the
# four-marker EC classifier and by the 20-gene classifier, in the four test studies. Singh 2024 is
# one stem-cell line, so its rows are differentiation days (read from the cell ID). The two McCauley
# 2023 graft rows come from classifier.R's by-organoid table.
#
# Inputs : scoring (four-marker classes of the test-study cells, heldout.R), predictions (20-gene classifier, classifier.R)
# Output : out (study, unit, n_rule, called_rule, called_cls)
args <- commandArgs(trailingOnly = TRUE)
a <- setNames(sub("^[^=]*=", "", args), sub("=.*$", "", args))

FAM <- c(gut_cell_atlas = "SF_ELMENTAITE2021_GUT_ATLAS", hickey = "SF_HICKEY2023_ADULT_EEC",
         egozi = "SF_EGOZI2021_DEVELOPMENTAL_EEC", gse238274 = "SF_INTESTINAL_GSE238274")
STUDY_IDS <- unname(FAM); SINGH <- "SF_INTESTINAL_GSE238274"
day <- function(id) sub("^GSE238274__(D[0-9]+):.*$", "\\1", id)

# four-marker classifier: EC cells per unit and how many are high-confidence EC (classified as EC)
sc <- read.delim(gzfile(a[["scoring"]]), stringsAsFactors = FALSE)
sc <- sc[sc$author_label_class == "EXPLICIT_EC", ]
sc$study <- FAM[sc$family]; sc$unit <- sub(" \\([^)]*\\)$", "", sc$independent_unit_id)
sc$unit[sc$study == SINGH] <- day(sc$cell_key[sc$study == SINGH])
sc$called <- sc$presentation_class == "HIGH_CONFIDENCE_EC"
ru <- aggregate(cbind(called = called, n = 1) ~ study + unit, data = transform(sc, n = 1), FUN = sum)

# 20-gene classifier at panel size 20
pr <- read.delim(gzfile(a[["predictions"]]), stringsAsFactors = FALSE)
pr <- pr[pr$panel_size == 20 & pr$primary_label == "POSITIVE" & pr$study %in% STUDY_IDS, ]
pr$call <- as.logical(pr$call); pr$unit <- sub(" \\([^)]*\\)$", "", pr$biological_unit_id)
pr$unit[pr$study == SINGH] <- day(pr$project_cell_id[pr$study == SINGH])
stopifnot(all(pr$unit[pr$study == SINGH] %in% c("D3", "D4", "D5", "D6")))
cu <- aggregate(cbind(called = call, n = 1) ~ study + unit, data = transform(pr, n = 1), FUN = sum)

du <- merge(cu, ru, by = c("study", "unit"), suffixes = c("_cls", "_rule"))
stopifnot(nrow(du) == nrow(ru), nrow(du) == nrow(cu), all(du$n_cls == du$n_rule))
du <- du[order(match(du$study, STUDY_IDS), ifelse(du$study == SINGH, as.numeric(sub("D", "", du$unit)), -du$n_rule)), ]
stopifnot(identical(as.vector(table(du$study)[STUDY_IDS]), c(10L, 9L, 6L, 4L)))
write.table(du[c("study", "unit", "n_rule", "called_rule", "called_cls")], a[["out"]], sep = "\t", quote = TRUE, row.names = FALSE)
