# Figure S2C-D: clusters of the 1,164 fitting cells and what marks each one.
#
# What it computes (score-blinded Harmony embedding from umap.R; the score genes were withheld)
#   Louvain clustering (Seurat, k = 20) at resolutions 0.1-2.0. The resolution kept is the one
#   recovering the most distinct reference enteroendocrine identities among resolutions whose
#   smallest cluster has >= 15 cells (ties: higher silhouette, then fewer clusters). One-vs-rest
#   Wilcoxon markers (positive, min.pct 0.25, logFC 0.25, adjusted p < 0.05). Each cluster gets
#   the reference identity with the highest module score, accepted only when that set's genes are
#   among the cluster's markers (otherwise the next-best supported set, or "unassigned").
#   The unblinded embedding is clustered at the same resolution for comparison.
#
#   The fitting cells are those of umap.R (without the seven pediatric Crohn's disease donors, roster).
#
# Inputs : matrix, features, cells, meta, roster, classes, embeddings (from umap.R)
# Outputs: summary (one row per cluster), assign (cluster of every cell)
#   Rscript clusters.R matrix=... features=... cells=... meta=... roster=... classes=... embeddings=... summary=... assign=...
args <- commandArgs(trailingOnly = TRUE)
a <- setNames(sub("^[^=]*=", "", args), sub("=.*$", "", args))
suppressPackageStartupMessages({library(Matrix); library(Seurat)})

R_BLIND <- "Harmony, score genes withheld"; R_FULL <- "Harmony on study family"
N_PCS <- 30L; K_PARAM <- 20L
RESOLUTIONS <- seq(0.1, 2.0, by = 0.1)
MIN_CLUSTER_CELLS <- 15L; MIN_CELLS_PER_GENE <- 10L
SEED_CLUSTER <- 2026082305L
REF_SETS <- list(
  "EC (enterochromaffin)" = c("TPH1", "SLC18A1", "FEV", "LMX1A", "DDC", "TAC1"),
  "L-cell" = c("GCG", "PYY", "INSL5"), "K-cell" = c("GIP"), "I-cell" = c("CCK"), "S-cell" = c("SCT"),
  "D-cell" = c("SST", "HHEX"), "M / ghrelin cell" = c("MLN", "GHRL"), "N-cell" = c("NTS"),
  "EEC progenitor" = c("NEUROG3", "NEUROD1", "PAX4", "SOX4", "HES6"),
  "Enterocyte" = c("FABP1", "APOA1", "APOA4", "ALPI", "RBP2"),
  "Goblet" = c("MUC2", "TFF3", "SPINK4", "CLCA1"), "Paneth" = c("LYZ", "DEFA5", "DEFA6", "REG3A"),
  "Tuft" = c("POU2F3", "TRPM5", "AVIL", "SH2D6"), "Stem / cycling" = c("LGR5", "OLFM4", "ASCL2", "MKI67", "TOP2A"))
EEC_IDENTITIES <- names(REF_SETS)[1:9]

# ---- load ------------------------------------------------------------------------------------
features <- read.delim(a[["features"]], stringsAsFactors = FALSE)
cells <- read.delim(a[["cells"]], stringsAsFactors = FALSE)
meta <- read.delim(gzfile(a[["meta"]]), stringsAsFactors = FALSE, check.names = FALSE)
meta <- meta[!duplicated(meta$capture_qualified_key), ]
idx <- match(cells$actual_cell_key, meta$capture_qualified_key)
stopifnot(!anyNA(idx))
n_all <- nrow(cells); n_feat <- nrow(features)
con <- file(a[["matrix"]], "rb")
X <- readBin(con, "double", n = as.numeric(n_all) * n_feat, size = 8, endian = "little")
close(con)
dim(X) <- c(n_feat, n_all)
keep <- meta$model_context[idx] == "primary_intestine"
ROSTER <- jsonlite::fromJSON(a[["roster"]])                         # the seven Crohn's donors are left out
crohn <- meta$study_family[idx] == ROSTER$family_id & meta$biological_unit[idx] %in% ROSTER$excluded_units
stopifnot(sum(keep & crohn) == 35L)
keep <- keep & !crohn
X <- X[, keep, drop = FALSE]; invisible(gc(verbose = FALSE))
mrow <- meta[idx[keep], ]
n_cells <- ncol(X)
stopifnot(n_cells == 1164L)
rownames(X) <- make.unique(gsub("_", "-", features$feature_symbol))
colnames(X) <- mrow$capture_qualified_key
D <- Matrix::Matrix(X, sparse = TRUE); rm(X); invisible(gc(verbose = FALSE))

cls <- read.delim(gzfile(a[["classes"]]), stringsAsFactors = FALSE)
cls <- cls[match(colnames(D), cls$cell_key), ]
stopifnot(!anyNA(cls$cell_key))
md <- data.frame(arm = mrow$arm_id, presentation_class = cls$presentation_class, candidate = cls$candidate,
                 row.names = colnames(D), stringsAsFactors = FALSE)

harm_all <- read.delim(gzfile(a[["embeddings"]]), stringsAsFactors = FALSE, check.names = FALSE)
dim_cols <- grep("^dim_", names(harm_all), value = TRUE)
harmony_of <- function(tag) {
  h <- harm_all[harm_all$embedding == tag, ]
  m <- as.matrix(h[, dim_cols]); rownames(m) <- h$cell_key
  m[colnames(D), , drop = FALSE]
}
H_main <- harmony_of(R_BLIND); H_other <- harmony_of(R_FULL)

# ---- cluster at every resolution -------------------------------------------------------------
obj <- Seurat::CreateSeuratObject(counts = D, meta.data = md, assay = "RNA", min.cells = 0, min.features = 0)
obj[["RNA"]]$data <- D
obj[["harmony"]] <- SeuratObject::CreateDimReducObject(embeddings = H_main, key = "harmony_", assay = "RNA")
obj <- Seurat::FindNeighbors(obj, reduction = "harmony", dims = 1:N_PCS, k.param = K_PARAM, verbose = FALSE)
set.seed(SEED_CLUSTER)
obj <- Seurat::FindClusters(obj, resolution = RESOLUTIONS, algorithm = 1L, random.seed = SEED_CLUSTER, verbose = FALSE)
res_cols <- sprintf("RNA_snn_res.%s", RESOLUTIONS)
d_harm <- dist(H_main)

# module score of each reference set per cluster: z-score of cluster mean expression across clusters
identity_of <- function(cl) {
  g <- unique(unlist(REF_SETS)); g <- g[g %in% rownames(D)]
  m <- as.matrix(D[g, , drop = FALSE])
  cm <- t(vapply(split(seq_len(n_cells), cl), function(i) Matrix::rowMeans(m[, i, drop = FALSE]), numeric(length(g))))
  z <- scale(cm); z[is.na(z)] <- 0
  sc <- vapply(REF_SETS, function(gs) { gs <- gs[gs %in% colnames(z)]
    if (!length(gs)) rep(NA_real_, nrow(z)) else rowMeans(z[, gs, drop = FALSE]) }, numeric(nrow(z)))
  list(identity = apply(sc, 1, function(r) names(REF_SETS)[which.max(r)]), scores = sc)
}

sweep_rows <- do.call(rbind, lapply(seq_along(RESOLUTIONS), function(i) {
  cl <- droplevels(factor(obj@meta.data[[res_cols[i]]]))
  tab <- table(cl)
  dom <- vapply(split(md$arm, cl), function(x) max(table(x)), numeric(1))
  sil <- if (nlevels(cl) > 1L) mean(cluster::silhouette(as.integer(cl), d_harm)[, "sil_width"]) else NA_real_
  idn <- identity_of(cl)
  data.frame(resolution = RESOLUTIONS[i], n_clusters = nlevels(cl), min_cluster_size = min(tab),
             mean_silhouette = sil, dominant_arm_share = sum(dom) / n_cells,
             n_eec_identities = length(intersect(unique(idn$identity), EEC_IDENTITIES)),
             n_identities = length(unique(idn$identity)), stringsAsFactors = FALSE)
}))
cand <- sweep_rows[sweep_rows$min_cluster_size >= MIN_CLUSTER_CELLS, ]
cand <- cand[order(-cand$n_eec_identities, -cand$mean_silhouette, cand$n_clusters), ]
CHOSEN <- cand$resolution[1]
cl_chosen <- droplevels(factor(obj@meta.data[[sprintf("RNA_snn_res.%s", CHOSEN)]]))
Seurat::Idents(obj) <- cl_chosen

# the unblinded embedding, clustered at the same resolution
obj_other <- Seurat::CreateSeuratObject(counts = D[1:2, , drop = FALSE], assay = "RNA", min.cells = 0, min.features = 0)
obj_other[["harmony"]] <- SeuratObject::CreateDimReducObject(embeddings = H_other, key = "harmony_", assay = "RNA")
obj_other <- Seurat::FindNeighbors(obj_other, reduction = "harmony", dims = 1:N_PCS, k.param = K_PARAM, verbose = FALSE)
set.seed(SEED_CLUSTER)
obj_other <- Seurat::FindClusters(obj_other, resolution = CHOSEN, algorithm = 1L, random.seed = SEED_CLUSTER, verbose = FALSE)
cl_other <- droplevels(factor(Seurat::Idents(obj_other)))

# ---- markers and identities ------------------------------------------------------------------
det <- Matrix::rowSums(D > 0)
mk_obj <- subset(obj, features = names(det)[det >= MIN_CELLS_PER_GENE])
markers <- Seurat::FindAllMarkers(mk_obj, only.pos = TRUE, min.pct = 0.25, logfc.threshold = 0.25,
                                  test.use = "wilcox", verbose = FALSE)
markers <- markers[markers$p_val_adj < 0.05, ]
markers <- markers[order(markers$cluster, -markers$avg_log2FC), ]

idn <- identity_of(cl_chosen)
support_for <- function(nm, mk) {
  gs <- REF_SETS[[nm]]; gs <- gs[gs %in% rownames(D)]
  hit <- mk[mk$gene %in% gs, ]; hit <- hit[order(-hit$avg_log2FC), ]
  list(n = nrow(hit), of = length(gs),
       txt = if (nrow(hit)) paste(sprintf("%s (%.1f)", hit$gene, hit$avg_log2FC), collapse = ", ") else "")
}
ver <- do.call(rbind, lapply(seq_along(levels(cl_chosen)), function(i) {
  k <- levels(cl_chosen)[i]; mk <- markers[markers$cluster == k, ]
  sc <- idn$scores[i, ]; ord <- order(sc, decreasing = TRUE)
  prop <- names(REF_SETS)[ord[1]]; ps <- support_for(prop, mk)
  alt <- NA_character_; as_ <- NULL
  for (j in ord) { cs <- support_for(names(REF_SETS)[j], mk); if (cs$n > 0) { alt <- names(REF_SETS)[j]; as_ <- cs; break } }
  if (ps$n > 0) { assigned <- prop; verdict <- "supported"; use <- ps
  } else if (!is.na(alt)) { assigned <- alt; verdict <- sprintf("FELL BACK from %s", prop); use <- as_
  } else { assigned <- "unassigned"; verdict <- sprintf("UNSUPPORTED, proposal was %s", prop); use <- list(n = 0L, of = ps$of, txt = "") }
  data.frame(proposed_identity = prop, score_margin = sc[ord[1]] - sc[ord[2]], assigned_identity = assigned,
             verdict = verdict, n_support = use$n, n_set = use$of, support_genes = use$txt, stringsAsFactors = FALSE)
}))
summ <- do.call(rbind, lapply(levels(cl_chosen), function(k) {
  i <- which(cl_chosen == k); at <- table(md$arm[i]); ct <- table(md$presentation_class[i])
  data.frame(cluster = k, n_cells = length(i),
             pct_high_confidence = 100 * as.numeric(ct["HIGH_CONFIDENCE_EC"]) / length(i),
             dominant_arm = names(at)[which.max(at)], dominant_arm_pct = 100 * max(at) / length(i),
             n_arms_present = sum(at > 0), n_markers = sum(markers$cluster == k),
             top_markers = paste(utils::head(markers$gene[markers$cluster == k], 12), collapse = ", "),
             stringsAsFactors = FALSE)
}))
summ <- cbind(summ, ver)
summ$pct_high_confidence[is.na(summ$pct_high_confidence)] <- 0

write.table(summ, a[["summary"]], sep = "\t", row.names = FALSE, quote = FALSE)
write.table(data.frame(cell_key = colnames(D), arm = md$arm, presentation_class = md$presentation_class,
                       cluster = as.character(cl_chosen), cluster_other_embedding = as.character(cl_other),
                       stringsAsFactors = FALSE),
            gzfile(a[["assign"]]), sep = "\t", row.names = FALSE, quote = FALSE)
