# Figure 1C: each study's endocrine cells on its own display embedding, colored by cell group.
# The per-study UMAP coordinates are an input (made per study from the counts; see inputs/figure1.tsv). This script
# keeps the studies of the Figure 1B roster and assigns every cell the color group the panel draws:
#   * primary study labels: EC -> intestinal EC (group gut_ec); SC-SIC; (SC-)beta, (SC-)alpha, (SC-)delta; every
#     other endocrine label -> other;
#   * Singh 2024 EC -> induced EC; McCauley 2023 "Endocrine_1" -> HIO graft EC (group hio_ec), its other labels -> other;
#   * Wang 2020, Burclaff 2022 and Yu 2021 have no EC in their primary study labels: cells the Figure 2 four-marker
#     classifier classifies as high-confidence EC (HIGH_CONFIDENCE_EC) are drawn as intestinal EC.
# Inputs: cells    per-study UMAP coordinates and primary study labels, one row per cell
#         calls    Figure 2 four-marker classifier results per cell (cell_key, presentation_class)
#         roles    Figure 1B role table as drawn (study roster)
# Output: out      study, cell_id, umap_1, umap_2, group (gzip-compressed TSV)
#   Rscript display_embedding.R cells=<...tsv.gz> calls=<...tsv.gz> roles=<study_roles_drawn.tsv> out=<...tsv.gz>
args <- commandArgs(trailingOnly = TRUE)
arg <- setNames(sub("^[^=]*=", "", args), sub("=.*$", "", args))
stopifnot(all(c("cells", "calls", "roles", "out") %in% names(arg)))

roster <- read.delim(arg[["roles"]], colClasses = "character", na.strings = character(0), encoding = "UTF-8")$study
cl <- read.delim(gzfile(arg[["cells"]]), check.names = FALSE, stringsAsFactors = FALSE)
cl$study[cl$study == "Elmentaite 2020 + Elmentaite 2021"] <- "Elmentaite 2020 + 2021"
cl$study[cl$study == "Yu Q 2021"] <- "Yu 2021"                      # paper-wide study name
cl <- cl[cl$study %in% roster, ]
stopifnot(nrow(cl) == 140805)

GRP <- c("EC" = "gut_ec", "SC-SIC" = "sc_sic", "SC-beta" = "beta", "Beta" = "beta", "SC-alpha" = "alpha", "Alpha" = "alpha",
         "SC-delta" = "delta", "Delta" = "delta", "Other EEC" = "other", "EEC (no subtype label)" = "other",
         "EEC progenitor" = "other", "Other endocrine" = "other", "PP" = "other", "Epsilon" = "other", "EC (marker-rule call)" = "other")
stopifnot(all(cl$class %in% names(GRP)))
cl$grp <- unname(GRP[cl$class]); cl$grp[cl$study == "Singh 2024" & cl$class == "EC"] <- "induced_ec"
mc <- cl$study == "McCauley 2023"; cl$grp[mc] <- ifelse(cl$author_label[mc] == "Endocrine_1", "hio_ec", "other")
f2 <- read.delim(gzfile(arg[["calls"]]), stringsAsFactors = FALSE)
NO_EC_LABEL <- c("Wang 2020", "Burclaff 2022", "Yu 2021")      # studies whose primary study labels have no EC
cl$grp[cl$study %in% NO_EC_LABEL & cl$cell_id %in% f2$cell_key[f2$presentation_class == "HIGH_CONFIDENCE_EC"]] <- "gut_ec"
stopifnot(identical(as.integer(table(factor(cl$study[cl$study %in% NO_EC_LABEL & cl$grp == "gut_ec"], NO_EC_LABEL))), c(40L, 36L, 173L)))

out <- data.frame(study = cl$study, cell_id = cl$cell_id, umap_1 = cl$umap_1, umap_2 = cl$umap_2, group = cl$grp)
con <- gzfile(arg[["out"]], "w")
write.table(out, con, sep = "\t", quote = FALSE, row.names = FALSE)
close(con)
