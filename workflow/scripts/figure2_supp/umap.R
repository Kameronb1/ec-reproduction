# Figure S2A: the 2-D map of the 1,199 fitting cells (and the 30-dimensional embeddings that
# clusters.R groups).
#
# What it computes
#   Raw counts are recovered by inverting log1p-CPM. Seurat: 3,000 variable genes (vst), scaling,
#   30 principal components, Harmony on study family (theta 2, lambda 1, sigma 0.1), and UMAP
#   (cosine, 30 neighbors, 500 epochs, min.dist 0.3). Three embeddings: unintegrated PCA,
#   Harmony on all genes, and Harmony with the four score genes (TPH1, SLC18A1, FEV, LMX1A)
#   withheld before variable-gene selection. Also: the share of each cell's 25 nearest neighbors
#   on the map that share its four-marker class, the maximum and detection of the 15 genes drawn on the map, and
#   each cell's log1p-CPM of those 15 genes (the maps are colored by value / gene maximum).
#
# Inputs : matrix, features, cells, meta, classes (per-cell classes from marker_fit.R)
# Outputs: coords, purity, scaling, embeddings, gene_expression
#   Rscript umap.R matrix=... features=... cells=... meta=... classes=... coords=... purity=... scaling=... embeddings=...
#     gene_expression=...
args <- commandArgs(trailingOnly = TRUE)
a <- setNames(sub("^[^=]*=", "", args), sub("=.*$", "", args))
suppressPackageStartupMessages({library(Matrix); library(Seurat)})

PANEL_GENES <- c("ARX","CBLN1","CCK","CHGA","CHGB","DDC","FEV","GAST","GCG","GHRL","GIP","HHEX",
                 "IAPP","IL20RA","INS","LMX1A","MLN","NKX6-1","NTS","PCSK1","PCSK2","PDX1","PPY",
                 "PYY","SCG5","SCT","SLC18A1","SST","TAC1","TPH1")
EC_CORE <- c("TPH1", "SLC18A1", "FEV", "LMX1A")
CLASSES <- c("HIGH_CONFIDENCE_EC", "LOW_CONFIDENCE_EC", "INSUFFICIENT_EC_EVIDENCE")
N_HVG <- 3000L; N_PCS <- 30L
SEED_PCA <- 2026082301L; SEED_HARMONY <- 2026082302L; SEED_UMAP <- 2026082304L
R_UNINT <- "Unintegrated"; R_FULL <- "Harmony on study family"; R_BLIND <- "Harmony, score genes withheld"
FEATURES <- c("TPH1", "SLC18A1", "FEV", "LMX1A", "CHGA", "NEUROG3", "DDC", "TAC1", "CCK", "GCG",
              "SST", "GHRL", "NTS", "OLFM4", "MKI67")

# ---- load the fitting cells ------------------------------------------------------------------
features <- read.delim(a[["features"]], stringsAsFactors = FALSE)
cells <- read.delim(a[["cells"]], stringsAsFactors = FALSE)
n_all <- nrow(cells); n_feat <- nrow(features)
meta <- read.delim(gzfile(a[["meta"]]), stringsAsFactors = FALSE, check.names = FALSE)
meta <- meta[!duplicated(meta$capture_qualified_key), ]
idx <- match(cells$actual_cell_key, meta$capture_qualified_key)
stopifnot(!anyNA(idx))
con <- file(a[["matrix"]], "rb")
X <- readBin(con, "double", n = as.numeric(n_all) * n_feat, size = 8, endian = "little")
close(con)
dim(X) <- c(n_feat, n_all)
keep <- meta$model_context[idx] == "primary_intestine"
X <- X[, keep, drop = FALSE]; invisible(gc(verbose = FALSE))
mrow <- meta[idx[keep], ]
n_cells <- ncol(X)
stopifnot(n_cells == 1199L)

# Seurat-safe unique feature names (21 symbols repeat; underscores become dashes)
feat <- make.unique(gsub("_", "-", features$feature_symbol))
stopifnot(!anyDuplicated(feat), identical(feat[match(PANEL_GENES, feat)], PANEL_GENES))
rownames(X) <- feat
colnames(X) <- mrow$capture_qualified_key

# ---- counts, by inverting log1p(raw * 1e6 / total); must equal the stored raw counts ----------
total <- as.numeric(mrow$cell_total_counts)
C <- expm1(X)
for (j in seq_len(n_cells)) C[, j] <- C[, j] * (total[j] / 1e6)
C <- round(C)
ref_counts <- vapply(PANEL_GENES, function(g) as.numeric(mrow[[paste0("raw_count__", g)]]), numeric(n_cells))
stopifnot(max(abs(t(C[match(PANEL_GENES, rownames(C)), , drop = FALSE]) - ref_counts)) == 0)
C <- Matrix::Matrix(C, sparse = TRUE)
D <- Matrix::Matrix(X, sparse = TRUE)
rm(X); invisible(gc(verbose = FALSE))

cls <- read.delim(gzfile(a[["classes"]]), stringsAsFactors = FALSE)
j <- match(colnames(C), cls$cell_key)
stopifnot(!anyNA(j), !anyDuplicated(j), nrow(cls) == n_cells)
cls <- cls[j, ]
md <- data.frame(arm = mrow$arm_id, study_family = mrow$study_family, biological_unit = mrow$biological_unit,
                 ec_score = cls$ec_core_score, markers_detected = cls$markers_detected,
                 presentation_class = cls$presentation_class, candidate = cls$candidate,
                 row.names = colnames(C), stringsAsFactors = FALSE)
obj <- Seurat::CreateSeuratObject(counts = C, meta.data = md, assay = "RNA", min.cells = 0, min.features = 0)
obj[["RNA"]]$data <- D                       # the stored log1p-CPM, not a recomputation
FAMILIES <- sort(unique(md$study_family))

build_rep <- function(universe, blind) {
  work <- subset(obj, features = universe)
  work <- Seurat::FindVariableFeatures(work, assay = "RNA", selection.method = "vst", loess.span = 0.3,
                                       clip.max = sqrt(n_cells), nfeatures = N_HVG, verbose = FALSE)
  hvg <- Seurat::VariableFeatures(work, assay = "RNA")
  stopifnot(length(hvg) == N_HVG, !blind || !any(hvg %in% EC_CORE))
  work <- Seurat::ScaleData(work, assay = "RNA", features = hvg, do.scale = TRUE, do.center = TRUE,
                            scale.max = 10, verbose = FALSE)
  set.seed(SEED_PCA)
  work <- Seurat::RunPCA(work, assay = "RNA", features = hvg, npcs = N_PCS, reduction.name = "pca",
                         seed.use = SEED_PCA, approx = TRUE, verbose = FALSE)
  work$fam <- factor(work$study_family, levels = FAMILIES)
  work[["RNA"]] <- split(work[["RNA"]], f = work$fam, layers = c("counts", "data", "scale.data"))
  set.seed(SEED_HARMONY)
  work <- Seurat::IntegrateLayers(object = work, method = Seurat::HarmonyIntegration, orig.reduction = "pca",
                                  assay = "RNA", layers = "data", scale.layer = "scale.data",
                                  new.reduction = "harmony", theta = 2, lambda = 1, sigma = 0.1, verbose = FALSE)
  work
}
umap_of <- function(work, reduction) {
  set.seed(SEED_UMAP)
  u <- Seurat::RunUMAP(work, reduction = reduction, dims = 1:N_PCS, umap.method = "uwot", n.neighbors = 30L,
                       n.components = 2L, metric = "cosine", n.epochs = 500L, min.dist = 0.3, spread = 1,
                       seed.use = SEED_UMAP, reduction.name = paste0("umap_", reduction), verbose = FALSE)
  Seurat::Embeddings(u, paste0("umap_", reduction))
}
full <- build_rep(rownames(obj), FALSE)
blind <- build_rep(setdiff(rownames(obj), EC_CORE), TRUE)
E <- list()
E[[R_UNINT]] <- umap_of(full, "pca")
E[[R_FULL]] <- umap_of(full, "harmony")
E[[R_BLIND]] <- umap_of(blind, "harmony")

# ---- same-class share among each cell's 25 nearest neighbors on the map ------------------------
knn_purity <- function(e, lab, k = 25L) {
  d <- as.matrix(dist(e)); diag(d) <- Inf
  nn <- t(apply(d, 1, function(r) order(r)[seq_len(k)]))
  vapply(seq_along(lab), function(i) mean(lab[nn[i, ]] == lab[i]), numeric(1))
}
purity_for <- function(e, tag) {
  lab <- md[rownames(e), "presentation_class"]; cand <- md[rownames(e), "candidate"]
  pc <- knn_purity(e, cand); pl <- knn_purity(e, lab)
  out <- rbind(data.frame(embedding = tag, class = "CANDIDATE_BEFORE_VETO", n = sum(cand), observed = mean(pc[cand]),
                          chance = mean(cand), stringsAsFactors = FALSE),
               do.call(rbind, lapply(CLASSES, function(cl) data.frame(embedding = tag, class = cl, n = sum(lab == cl),
                          observed = mean(pl[lab == cl]), chance = mean(lab == cl), stringsAsFactors = FALSE))))
  out$enrichment <- out$observed / out$chance
  out
}
pur_tab <- rbind(purity_for(E[[R_BLIND]], R_BLIND), purity_for(E[[R_FULL]], R_FULL))

# ---- coordinates, with a per-map centering and single-scalar rescaling used for drawing ---------
dat <- do.call(rbind, lapply(names(E), function(r) data.frame(reduction = r, cell_key = rownames(E[[r]]),
  umap_1 = E[[r]][, 1], umap_2 = E[[r]][, 2], md[rownames(E[[r]]), ], row.names = NULL, stringsAsFactors = FALSE)))
dat$disp_1 <- dat$umap_1; dat$disp_2 <- dat$umap_2
for (r in names(E)) {
  i <- dat$reduction == r
  cx <- mean(range(dat$umap_1[i])); cy <- mean(range(dat$umap_2[i]))
  sc <- max(diff(range(dat$umap_1[i])), diff(range(dat$umap_2[i]))) / 2
  dat$disp_1[i] <- (dat$umap_1[i] - cx) / sc; dat$disp_2[i] <- (dat$umap_2[i] - cy) / sc
}
class_lab <- factor(c(HIGH_CONFIDENCE_EC = "High-confidence EC", LOW_CONFIDENCE_EC = "Low-confidence EC",
                      INSUFFICIENT_EC_EVIDENCE = "Insufficient EC evidence")[dat$presentation_class],
                    levels = c("High-confidence EC", "Low-confidence EC", "Insufficient EC evidence"))
dat <- dat[order(class_lab, decreasing = TRUE), ]
write.table(dat[, c("reduction", "cell_key", "umap_1", "umap_2", "disp_1", "disp_2", "arm", "study_family",
                    "presentation_class", "ec_score", "markers_detected")],
            gzfile(a[["coords"]]), sep = "\t", row.names = FALSE, quote = FALSE)
write.table(pur_tab, a[["purity"]], sep = "\t", row.names = FALSE, quote = FALSE)

GMAX <- apply(as.matrix(D[FEATURES, , drop = FALSE]), 1, max)
ndet <- as.integer(Matrix::rowSums(D[FEATURES, , drop = FALSE] > 0))
write.table(data.frame(gene = FEATURES, panel_row = rep(1:3, each = 5), max_log1p_cpm = as.numeric(GMAX[FEATURES]),
                       strip_label = sprintf("%.1f", GMAX[FEATURES]), n_cells_detecting = ndet,
                       pct_cells_detecting = round(100 * ndet / n_cells, 2), stringsAsFactors = FALSE),
            a[["scaling"]], sep = "\t", row.names = FALSE, quote = FALSE)

# ---- per-cell values of the 15 genes, in the coordinates' cell order, 17 significant digits ------
stopifnot(identical(match(FEATURES, feat), match(FEATURES, features$feature_symbol)))   # first matrix column of each symbol
keys <- unique(dat$cell_key); stopifnot(length(keys) == n_cells)
GE <- t(as.matrix(D[FEATURES, keys, drop = FALSE]))
write.table(data.frame(cell_key = keys, apply(GE, 2, sprintf, fmt = "%.17g"), check.names = FALSE),
            a[["gene_expression"]], sep = "\t", row.names = FALSE, quote = FALSE)

# ---- the 30-dimensional reductions, read by clusters.R ----------------------------------------
SRC <- list(); SRC[[R_UNINT]] <- list(full, "pca"); SRC[[R_FULL]] <- list(full, "harmony"); SRC[[R_BLIND]] <- list(blind, "harmony")
harm <- do.call(rbind, lapply(names(SRC), function(tag) {
  h <- Seurat::Embeddings(SRC[[tag]][[1]], SRC[[tag]][[2]])
  colnames(h) <- sprintf("dim_%02d", seq_len(N_PCS))
  cbind(data.frame(embedding = tag, cell_key = rownames(h), stringsAsFactors = FALSE), as.data.frame(h))
}))
write.table(harm, gzfile(a[["embeddings"]]), sep = "\t", row.names = FALSE, quote = FALSE)
