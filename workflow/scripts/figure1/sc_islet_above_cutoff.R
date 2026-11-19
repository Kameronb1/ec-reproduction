# Figure 1E: the share of each SC-islet study's SC-SIC and SC-beta cells whose four-marker score is above the cutoff of
# the label-free four-marker EC classifier (fitted on primary intestinal cells, Supplementary Figure S2A).
# The scores and the cutoff are results of the classifier (Supplementary Figures S2 and S3); the cutoff must equal the
# fitted constants' ec_core_cutoff. Share = round(100 * mean(score > cutoff)), whole percent.
# Inputs:  scores     study / class / score of every scored mature endocrine cell of the six SC-islet studies (gzipped)
#          cutoff     the cutoff (one column, cutoff)
#          constants  the classifier's fitted constants (JSON)
# Output:  out        study / sc_sic_pct / sc_beta_pct, rows in the panel's order, top to bottom
#   Rscript sc_islet_above_cutoff.R scores=<...tsv.gz> cutoff=<...tsv> constants=<...json> out=<sc_islet_above_cutoff.tsv>
args <- commandArgs(trailingOnly = TRUE)
arg <- setNames(sub("^[^=]*=", "", args), sub("=.*$", "", args))
stopifnot(all(c("scores", "cutoff", "constants", "out") %in% names(arg)))

ORDER <- c("Augsornworawat 2020", "Balboa 2022", "Du 2022", "Augsornworawat 2023 (nuclei)", "Zhu 2023", "Siehler 2024")
cells <- read.delim(gzfile(arg[["scores"]]), stringsAsFactors = FALSE)
cutoff <- read.delim(arg[["cutoff"]], stringsAsFactors = FALSE)$cutoff
fitted <- jsonlite::fromJSON(arg[["constants"]])$ec_core_cutoff
if (length(cutoff) != 1 || !isTRUE(all.equal(cutoff, fitted, tolerance = 1e-15)))
  stop("the cutoff differs from the fitted constants' ec_core_cutoff")
if (!setequal(unique(cells$study), ORDER)) stop("the score table's studies are not the six SC-islet studies")
if (anyNA(cells$score)) stop("a scored cell has no score")

above <- function(cls) vapply(ORDER, function(s) round(100 * mean(cells$score[cells$study == s & cells$class == cls] > cutoff)), 0)
out <- data.frame(study = ORDER, sc_sic_pct = above("SC_SIC"), sc_beta_pct = above("SC_BETA"))
write.table(out, arg[["out"]], sep = "\t", quote = FALSE, row.names = FALSE)
