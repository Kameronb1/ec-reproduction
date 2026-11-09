# Figure 3B: the label-free four-marker EC classifier (fitted on primary intestinal tissue in Figure 2, applied
# unchanged) on primary-study-labeled mature endocrine cells of the six 10x SC-islet studies, one row per study.
#
# Keeps every scored cell of the six 10x studies (the inDrops study Veres 2019, GSE114412, is left out), in file
# order, with its study name and primary-study class (SC_SIC, SC_BETA, SC_ALPHA, ...) and its EC marker score
# (mean z-score of TPH1, SLC18A1, FEV, LMX1A with Figure 2's fixed constants). The panel draws each study's
# score distribution for SC-SIC, SC-beta and the other endocrine cells, the dashed line at Figure 2A's cutoff,
# and, per study, the share of SC-SIC and of SC-beta cells scoring above the cutoff (rounded to whole percent).
#
# Inputs:  scores=  per-cell scores (dataset, protocol, target, score, called, ...)
#          fit=     Figure 2A's fitted constants (ec_core_cutoff)
#          census=  per-sample census of primary study labels (check: SC-SIC cells per study equal Figure 3A's)
# Outputs: cells=   study / class / score, gzipped (scores written with 15 significant digits, 17 where 15 do
#                   not read back as the same number)
#          cutoff=  the cutoff
#          above=   study / class / n_cells / n_above / pct (rows in the panel's study order, most SC-SIC cells first)
#
#   Rscript score_by_study.R scores=... fit=... census=... cells=... cutoff=... above=...

args <- commandArgs(trailingOnly = TRUE)
arg <- setNames(sub("^[^=]*=", "", args), sub("=.*", "", args))
suppressPackageStartupMessages(library(jsonlite))

STUDY <- c(GSE151117 = "Augsornworawat 2020", GSE167880 = "Balboa 2022", GSE185036 = "Du 2022",
           GSE199636 = "Augsornworawat 2023 (nuclei)", GSE202497 = "Zhu 2023", GSE226346 = "Siehler 2024")
HELD_OUT <- "GSE114412"

num <- function(x) { s <- sprintf("%.15g", x); bad <- !is.na(x) & suppressWarnings(as.numeric(s)) != x
  s[bad] <- sprintf("%.17g", x[bad]); s[is.na(x)] <- "NA"; s }

cutoff <- fromJSON(arg[["fit"]])$ec_core_cutoff
if (!is.numeric(cutoff) || length(cutoff) != 1) stop("no single ec_core_cutoff in the fitted constants")

d <- read.delim(gzfile(arg[["scores"]]), stringsAsFactors = FALSE)
if (!setequal(unique(d$dataset), c(names(STUDY), HELD_OUT))) stop("the scores table holds other studies")
d <- d[d$dataset %in% names(STUDY), ]
cells <- data.frame(study = unname(STUDY[d$dataset]), class = d$target, score = as.numeric(d$score))
if (anyNA(cells$score)) stop("a scored cell has no score")
if (any(cells$score == cutoff)) stop("a score equals the cutoff")

census <- read.delim(arg[["census"]], stringsAsFactors = FALSE)
census <- census[census$dataset %in% names(STUDY) & census$target == "SC_SIC", ]
sic_census <- tapply(census$n_cells, unname(STUDY[census$dataset]), sum)
sic_scored <- table(cells$study[cells$class == "SC_SIC"])
if (!setequal(names(sic_census), names(sic_scored)) || any(sic_scored[names(sic_census)] != sic_census))
  stop("SC-SIC cells per study differ from the census (Figure 3A)")

out <- gzfile(arg[["cells"]], "wb")
write.table(transform(cells, score = num(score)), out, sep = "\t", quote = FALSE, row.names = FALSE, eol = "\n")
close(out)
out <- file(arg[["cutoff"]], "wb"); writeLines(c("cutoff", num(cutoff)), out); close(out)

# share above the cutoff, as the panel computes it: round(100 * mean(score > cutoff)) per study and class
order_studies <- names(sort(sic_scored, decreasing = TRUE))
above <- do.call(rbind, lapply(order_studies, function(st) do.call(rbind, lapply(c("SC_SIC", "SC_BETA"), function(cl) {
  s <- cells$score[cells$study == st & cells$class == cl]
  data.frame(study = st, class = cl, n_cells = length(s), n_above = sum(s > cutoff), pct = round(100 * mean(s > cutoff)))
}))))
out <- file(arg[["above"]], "wb")
write.table(above, out, sep = "\t", quote = FALSE, row.names = FALSE, eol = "\n")
close(out)
