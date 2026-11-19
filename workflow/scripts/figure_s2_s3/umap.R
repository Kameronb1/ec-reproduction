# Figure S2C: the 2-D map of the 1,164 fitting cells (and the 30-dimensional embeddings that
# clusters.R groups, S2C-D).
#
# What it computes
#   The fitting cells without the seven pediatric Crohn's disease donors of Elmentaite 2020 (roster).
#   Raw counts are recovered by inverting log1p-CPM. Seurat: 3,000 variable genes (vst), scaling,
#   30 principal components, Harmony on study family (theta 2, lambda 1, sigma 0.1), and UMAP
#   (cosine, 30 neighbors, 500 epochs, min.dist 0.3). Three embeddings: unintegrated PCA,
#   Harmony on all genes, and Harmony with the four score genes (TPH1, SLC18A1, FEV, LMX1A)
#   withheld before variable-gene selection (the one drawn).
#
# Inputs : matrix, features, cells, meta, roster, classes (per-cell classes from marker_fit.R)
# Outputs: coords, embeddings
#   Rscript umap.R matrix=... features=... cells=... meta=... roster=... classes=... coords=... embeddings=...
args <- commandArgs(trailingOnly = TRUE)
a <- setNames(sub("^[^=]*=", "", args), sub("=.*$", "", args))
suppressPackageStartupMessages({library(Matrix); library(Seurat)})

PANEL_GENES <- c("ARX","CBLN1","CCK","CHGA","CHGB","DDC","FEV","GAST","GCG","GHRL","GIP","HHEX",
                 "IAPP","IL20RA","INS","LMX1A","MLN","NKX6-1","NTS","PCSK1","PCSK2","PDX1","PPY",
                 "PYY","SCG5","SCT","SLC18A1","SST","TAC1","TPH1")
EC_CORE <- c("TPH1", "SLC18A1", "FEV", "LMX1A")
N_HVG <- 3000L; N_PCS <- 30L
SEED_PCA <- 2026082301L; SEED_HARMONY <- 2026082302L; SEED_UMAP <- 2026082304L
R_UNINT <- "Unintegrated"; R_FULL <- "Harmony on study family"; R_BLIND <- "Harmony, score genes withheld"

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
ROSTER <- jsonlite::fromJSON(a[["roster"]])                         # the seven Crohn's donors are left out
crohn <- meta$study_family[idx] == ROSTER$family_id & meta$biological_unit[idx] %in% ROSTER$excluded_units
stopifnot(sum(keep & crohn) == 35L)
keep <- keep & !crohn
X <- X[, keep, drop = FALSE]; invisible(gc(verbose = FALSE))
mrow <- meta[idx[keep], ]
n_cells <- ncol(X)
stopifnot(n_cells == 1164L)

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

# ---- the 30-dimensional reductions, read by clusters.R ----------------------------------------
SRC <- list(); SRC[[R_UNINT]] <- list(full, "pca"); SRC[[R_FULL]] <- list(full, "harmony"); SRC[[R_BLIND]] <- list(blind, "harmony")
harm <- do.call(rbind, lapply(names(SRC), function(tag) {
  h <- Seurat::Embeddings(SRC[[tag]][[1]], SRC[[tag]][[2]])
  colnames(h) <- sprintf("dim_%02d", seq_len(N_PCS))
  cbind(data.frame(embedding = tag, cell_key = rownames(h), stringsAsFactors = FALSE), as.data.frame(h))
}))
write.table(harm, gzfile(a[["embeddings"]]), sep = "\t", row.names = FALSE, quote = FALSE)
