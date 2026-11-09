# Figure S17E: are CDX-family motif sites open in each study's EC chromatin, and how much CDX2 RNA does that study's EC
# population carry?
#
# One row per cohort with chromatin data. Motif sites (z): CDX1/CDX2/HOXA10 motif-family openness, EC minus sibling
# populations, mean chromVAR z over the family's three motifs (one motif per member) and the siblings. Augsornworawat 2023
# (GSE199636), Hickey 2023 and Singh 2024 come from the motif-family table; Zhu 2023's openness is the mean of the same
# three motifs over days 32 and 39 from the Zhu 2023 chromatin table (its GSE199636 column must reproduce the
# motif-family value).
# RNA: CDX2 log2(1 + CPM) in the study's own EC population (motif-family table), and for Zhu 2023 the mean over its
# SC-SIC samples (GSE202497, samples with at least 20 cells). McCauley 2023 (HIO graft EC, no chromatin data) is drawn on
# the RNA axis only, from its CDX2 row of the regulator RNA table for the HIO samples (not computed here).
#
# Inputs:  family=<motif_family_vs_rna.tsv>  families=<motif_families.tsv>  zhu=<zhu_candidate_tf_directions.tsv>
#          rna=<tf_expression_by_sample.tsv>
# Output:  out=<f5_s4_cdx_sites_vs_rna.tsv>  (cohort, kind, z, rna)
args <- commandArgs(trailingOnly = TRUE); a <- setNames(sub("^[^=]*=", "", args), sub("=.*", "", args))
rtsv <- function(p) read.delim(p, stringsAsFactors = FALSE)

FAMN <- "CDX1/CDX2/HOXA10"
FM <- rtsv(a[["families"]]); mem <- FM$tf[FM$family == FAMN]
stopifnot(setequal(mem, c("CDX1", "CDX2", "HOXA10")), !anyDuplicated(mem))     # one motif per member, so a TF mean = its motif
V <- rtsv(a[["family"]]); V <- V[V$family == FAMN, ]
ZT <- rtsv(a[["zhu"]]); ZT <- ZT[ZT$tf %in% mem, ]
# Zhu's GSE199636 column must reproduce the motif-family value, or the two averages are not the same metric
stopifnot(abs(mean(ZT$gse199636) - V$z[V$cohort == "gse199636"]) < 0.005)
TR <- rtsv(a[["rna"]]); TR <- TR[TR$tf == "CDX2" & as.logical(TR$min20), ]
sic <- TR[TR$population == "SC_SIC" & TR$cohort == "sc", ]
per_study <- tapply(sic$mean_log2cpm, sic$family, mean)                         # each study: mean over its samples
zhu_rna <- per_study[["GSE202497"]]
D <- data.frame(cohort = c(V$cohort, "zhu"), z = c(V$z, mean(ZT$zhu_D32_D39)), rna = c(V$rna, zhu_rna))
D$kind <- ifelse(D$cohort %in% c("gse199636", "zhu"), "SC-SIC (SC-islet)", ifelse(D$cohort == "hickey", "adult gut EC", "induced EC"))
write.table(D[, c("cohort", "kind", "z", "rna")], a[["out"]], sep = "\t", quote = TRUE, row.names = FALSE, na = "NA")
