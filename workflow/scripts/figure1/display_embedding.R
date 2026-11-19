# Figure 1C: the endocrine cells of the six SC-islet studies, each study on its own display embedding, colored by
# published label, and the SC-SIC share of each study's labeled endocrine cells.
# The per-study UMAP coordinates are an input (made per study from the counts; see inputs/figure1.tsv). This script
# keeps the SC-islet cells and gives every cell the color group the panel draws: SC-SIC (group sc_sic), SC-beta (beta),
# SC-alpha (alpha), SC-delta (delta), every other endocrine label (other). The SC-SIC share is the SC-SIC cells over all
# cells of the study's UMAP, other endocrine cells included.
# Inputs:  cells   per-study UMAP coordinates and published labels, one row per cell
# Outputs: out     study, cell_id, umap_1, umap_2, group (gzip-compressed TSV), cells in input order
#          share   study, sc_sic_cells, tile_cells, share_pct (percent, 2 decimals), printed (as drawn, 1 decimal);
#                  studies in the panel's order
#   Rscript display_embedding.R cells=<...tsv.gz> out=<...tsv.gz> share=<sc_sic_share.tsv>
args <- commandArgs(trailingOnly = TRUE)
arg <- setNames(sub("^[^=]*=", "", args), sub("=.*$", "", args))
stopifnot(all(c("cells", "out", "share") %in% names(arg)))

STUDIES <- c("Augsornworawat 2020", "Balboa 2022", "Du 2022", "Augsornworawat 2023", "Zhu 2023", "Siehler 2024")
GRP <- c("SC-SIC" = "sc_sic", "SC-beta" = "beta", "SC-alpha" = "alpha", "SC-delta" = "delta", "Other endocrine" = "other")

cl <- read.delim(gzfile(arg[["cells"]]), check.names = FALSE, stringsAsFactors = FALSE)
cl <- cl[cl$context == "SC-islet", ]
stopifnot(setequal(unique(cl$study), STUDIES), all(cl$class %in% names(GRP)), nrow(cl) == 55586)
cl$group <- unname(GRP[cl$class])

out <- data.frame(study = cl$study, cell_id = cl$cell_id, umap_1 = cl$umap_1, umap_2 = cl$umap_2, group = cl$group)
con <- gzfile(arg[["out"]], "w")
write.table(out, con, sep = "\t", quote = FALSE, row.names = FALSE)
close(con)

share <- vapply(STUDIES, function(s) 100 * mean(out$group[out$study == s] == "sc_sic"), 0)
tab <- data.frame(study = STUDIES, sc_sic_cells = vapply(STUDIES, function(s) sum(out$study == s & out$group == "sc_sic"), 0L),
                  tile_cells = vapply(STUDIES, function(s) sum(out$study == s), 0L), share_pct = round(share, 2),
                  printed = sprintf("%.1f%%", share))
write.table(tab, arg[["share"]], sep = "\t", quote = FALSE, row.names = FALSE)
