# Figure S10 (A, B) and Figure S11: the genes of the three EC gene groups, row by row, with the values the heatmaps show.
#   S10  G1 (EC genes low in SC-SIC) and G2 (EC genes SC-SIC have; split by its Fig. 3D result)
#   S11  G3 (SC-SIC's own program), in two side-by-side blocks of 43 and 42 genes
# Each row: five fold changes (EC vs other EEC in the three intestinal studies; SC-SIC vs other SC endocrine cells, median
# of six studies; SC-SIC minus intestinal EC from the Fig. 3D donor-block fit) and nine per-study expression values, log2(1 + CPM) of
# the study's pooled counts (intestinal EC pooled over donors, Elmentaite 2020 + 2021 together and its Crohn's donors left
# out; SC-SIC pooled over the study's samples), plus how many fold changes are beyond 6 (drawn clipped at 6).
# Row order: G1 and G3 by the size of the Fig. 3D t statistic (largest first); G2 by the Fig. 3D difference.
# Inputs:  groups = gene_groups_definitions.tsv.gz (main Fig. 3F gene groups)
#          gut_sums, gut_groups, gut_roster, sc_sums, sc_groups (per-sample pseudobulk sums and group tables;
#          gut_roster lists the excluded Crohn's donors)
# Outputs: s6 = table for S10, s7 = table for S11
a <- commandArgs(TRUE); A <- setNames(sub("^[^=]*=", "", a), sub("=.*", "", a))
rd <- function(p) read.delim(if (grepl("\\.gz$", p)) gzfile(p) else p, check.names = FALSE, stringsAsFactors = FALSE)
G <- rd(A[["groups"]])
# per-study log2(1 + CPM) of pooled counts: counts summed over the study's units, divided by the summed full-cell totals
rds <- function(f) { m <- read.delim(gzfile(f), check.names = FALSE); g <- m$gene; m <- as.matrix(m[, -1]); rownames(m) <- g; m }
gs <- rds(A[["gut_sums"]]); gg <- rd(A[["gut_groups"]]); ss <- rds(A[["sc_sums"]]); sg <- rd(A[["sc_groups"]])
crohns <- jsonlite::fromJSON(A[["gut_roster"]])$excluded_units
gg <- gg[gg$class == "POSITIVE" & !(gg$family_id == "SF_ELMENTAITE2020_DEVGUT" & gg$unit %in% crohns), ]
sg <- sg[sg$class == "SC_SIC", ]
gg$key <- c(SF_ELMENTAITE2020_DEVGUT = "ELM_PAIR", SF_ELMENTAITE2021_GUT_ATLAS = "ELM_PAIR", SF_HICKEY2023_ADULT_EEC = "HICKEY2023",
            SF_EGOZI2021_DEVELOPMENTAL_EEC = "EGOZI2021")[gg$family_id]
pool <- function(m, grp, key) vapply(split(grp$group, key), function(g) log2(1 + rowSums(m[, g, drop = FALSE]) / sum(grp$lib[match(g, grp$group)]) * 1e6),
                                     numeric(nrow(m)))
sh <- intersect(rownames(gs), rownames(ss))                     # genes in both count tables
LVM <- cbind(pool(gs[sh, ], gg, gg$key), pool(ss[sh, ], sg, sg$family_id))

GROUPS <- c(G1 = "G1_EC_GENES_SC_SIC_LACK", G2 = "G2_EC_GENES_SC_SIC_HAVE", G3 = "G3_SC_SIC_OWN_PROGRAM")
stopifnot(sum(G$group == GROUPS[["G1"]]) == 28, sum(G$group == GROUPS[["G2"]]) == 34, sum(G$group == GROUPS[["G3"]]) == 84)

# studies of the level columns: three intestinal EC studies, then the six SC-islet studies in main Fig. 3A order
LV_KEYS <- c("ELM_PAIR", "HICKEY2023", "EGOZI2021", "GSE151117", "GSE167880", "GSE202497", "GSE226346", "GSE185036", "GSE199636")
LV_LAB <- c("Elmentaite 2020 + 2021 (donors pooled)", "Hickey 2023 (nuclei)", "Egozi 2021", "Augsornworawat 2020", "Balboa 2022",
            "Zhu 2023", "Siehler 2024", "Du 2022", "Augsornworawat 2023 (nuclei)")
stopifnot(setequal(colnames(LVM), LV_KEYS), !anyNA(gg$key))

G$paneld_logFC <- G$panel5_logFC
FC_F <- c("gut_lfc_Elmentaite 2020 + 2021", "gut_lfc_Hickey 2023", "gut_lfc_Egozi 2021", "sic_in_dish_lfc", "paneld_logFC")
FC_COL <- c("fc_EC vs other EEC: Elmentaite 2020 + 2021 (mean of the 2 datasets)", "fc_EC vs other EEC: Hickey 2023",
            "fc_EC vs other EEC: Egozi 2021", "fc_SC-SIC vs other SC endocrine (median of 6 studies)", "fc_SC-SIC minus gut EC (Fig. 3D)")
CAPV <- 6

g1 <- G[G$group == GROUPS[["G1"]], ]; g1 <- g1[order(-abs(g1$paneld_t)), ]
g2 <- G[G$group == GROUPS[["G2"]], ]; g2 <- g2[order(-g2$paneld_logFC), ]
g2h <- g2[g2$paneld_call == "SC_SIC_HIGHER", ]; g2n <- g2[g2$paneld_call == "NO_CALL", ]
stopifnot(nrow(g2h) == 9, nrow(g2n) == 25)
g3 <- G[G$group == GROUPS[["G3"]], ]; g3 <- g3[order(-abs(g3$paneld_t)), ]
g3l <- g3[1:42, ]; g3r <- g3[43:84, ]

table_of <- function(blocks) do.call(rbind, lapply(names(blocks), function(b) {
  g0 <- blocks[[b]]
  stopifnot(all(g0$gene %in% rownames(LVM)))
  lv <- LVM[g0$gene, LV_KEYS, drop = FALSE]
  d <- data.frame(block = b, row = seq_len(nrow(g0)), gene = g0$gene, stringsAsFactors = FALSE)
  for (j in 1:5) d[[FC_COL[j]]] <- g0[[FC_F[j]]]
  for (j in 1:9) d[[paste0("expression_log2(1+CPM): ", LV_LAB[j])]] <- lv[, j]
  d$fc_cells_beyond_6 <- rowSums(abs(as.matrix(g0[, FC_F])) > CAPV)
  d
}))
wt <- function(d, p) write.table(d, p, sep = "\t", quote = TRUE, row.names = FALSE, fileEncoding = "UTF-8")
wt(table_of(list(G1 = g1, `G2 higher in SC-SIC, Fig. 3D FDR < 0.05` = g2h, `G2 no significant Fig. 3D difference` = g2n)), A[["s6"]])
wt(table_of(list(`G3 left` = g3l, `G3 right` = g3r)), A[["s7"]])
