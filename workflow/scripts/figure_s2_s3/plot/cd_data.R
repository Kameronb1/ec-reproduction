# Figure S2C-D data, read by S2C.R and S2D.R: the marker-free map of the 1,164 fitting cells (the embedding with the
# four score genes withheld), each cell's cluster and four-marker class, and one row per cluster.
EMB <- "Harmony, score genes withheld"
um <- read_table("umap_coordinates.tsv.gz"); um <- um[um$reduction == EMB, ]
cc <- read_table("cell_clusters.tsv.gz")
cls <- read_table("fitting_cell_classes.tsv.gz")
summ <- read_table("cluster_summary.tsv")
stopifnot(nrow(um) == 1164L, nrow(cc) == 1164L, setequal(um$cell_key, cc$cell_key), setequal(cc$cell_key, cls$cell_key))
cc$presentation_class <- cls$presentation_class[match(cc$cell_key, cls$cell_key)]
v <- merge(um[c("cell_key", "disp_1", "disp_2")], cc[c("cell_key", "arm", "presentation_class", "cluster")], by = "cell_key")
summ <- summ[order(summ$cluster), ]
NCL <- nrow(summ); stopifnot(NCL <= length(CLUSTER_PAL), identical(sort(unique(cc$cluster)), summ$cluster))
CL_COL <- setNames(CLUSTER_PAL[seq_len(NCL)], summ$cluster)
ID_SHORT <- c("EC (enterochromaffin)" = "EC", "unassigned" = "Unassigned")
summ$ident <- ifelse(summ$assigned_identity %in% names(ID_SHORT), ID_SHORT[summ$assigned_identity], summ$assigned_identity)
summ$genes <- gsub(" \\([0-9.]+\\)", "", summ$support_genes); summ$genes[is.na(summ$genes) | summ$genes == ""] <- "none"
ORD <- c("INSUFFICIENT_EC_EVIDENCE", "LOW_CONFIDENCE_EC", "HIGH_CONFIDENCE_EC")
