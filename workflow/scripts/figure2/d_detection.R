# Figure 2D dot sizes: for each dataset of panel D and each of its 20 genes, the share of EC cells and of
# other EEC cells with at least one read of the gene (raw count > 0), next to the drawn effect.
#
# What it computes
#   Seven gene-ranking datasets (training object): EC and other = the training labels (primary study
#     labels where the study has them, else the four-marker classifier's results); the shares use every
#     labeled cell of the dataset; effect = mean over donors/samples (>= 5 cells per donor-class) of EC
#     minus other mean log-normalized expression, as in gene_panel_effects.R.
#   Singh 2024: EC and other by its primary study labels; effect = difference of the mean log-normalized
#     expression over cells.
#   McCauley 2023 HIO grafts: Endocrine_1 = EC, Endocrine_2 = other; expression log1p(counts / library
#     total x 1e6) with each cell's whole-library total. Shares and effects per HIO sample group (wt1 =
#     libraries wt1.1 + wt1.2, wt2), then the mean of the two groups; cell counts are summed.
#   Controls: the effects equal panel D's (D_effect_heatmap.tsv) and, per HIO sample group,
#   hio_D_hio_effects.tsv.
#
# Inputs : analysis (training RDS), secondary (Singh 2024 RDS), sv_labels (Singh 2024 primary study
#          labels as compared), counts, genes, cells (the GSE214852 endocrine export), hio (four-marker
#          classifier results from four_marker_fit.R, for library totals), heatmap (D_effect_heatmap.tsv),
#          hio_d (hio_D_hio_effects.tsv)
# Outputs: out (D_detection.tsv)
args <- commandArgs(trailingOnly = TRUE)
a <- setNames(sub("^[^=]*=", "", args), sub("=.*$", "", args))
source(file.path(dirname(sub("^--file=", "", grep("^--file=", commandArgs(), value = TRUE))), "classifier_functions.R"))
SINGH <- "SF_INTESTINAL_GSE238274"; HIO <- "SF_GSE214852_HIO_GRAFT"

D <- read.delim(a[["heatmap"]], stringsAsFactors = FALSE)
genes <- unique(D$gene_symbol[order(D$panel_rank)]); stopifnot(length(genes) == 20L)

# ---- the seven gene-ranking datasets -------------------------------------------------------------
x <- readRDS(a[["analysis"]]); m <- x$metadata; labels <- training_labels(m)
u <- unit_class_means(x$normalized, m, labels, rep(TRUE, nrow(m)), genes, MIN_CELLS_PER_UNIT_CLASS)
um <- u$metadata
gi <- match(genes, x$genes); stopifnot(!anyNA(gi))
det <- function(mat, cells) as.numeric(Matrix::rowMeans(mat[gi, cells, drop = FALSE] > 0))
res <- do.call(rbind, lapply(sort(unique(um$family_id)), function(f) {
  p <- um$family_id == f & um$label == "POSITIVE"; n <- um$family_id == f & um$label == "NEGATIVE"
  cp <- m$family_id == f & labels$label == "POSITIVE"; cn <- m$family_id == f & labels$label == "NEGATIVE"
  data.frame(family_id = f, gene_symbol = genes,
             effect = rowMeans(u$values[, p, drop = FALSE]) - rowMeans(u$values[, n, drop = FALSE]),
             det_ec = det(x$counts, which(cp)), det_other = det(x$counts, which(cn)), ec_cells = sum(cp),
             other_cells = sum(cn), stringsAsFactors = FALSE)
}))
rm(x); invisible(gc(verbose = FALSE))

# ---- Singh 2024 ----------------------------------------------------------------------------------
s <- readRDS(a[["secondary"]])
lab <- read.delim(gzfile(a[["sv_labels"]]), stringsAsFactors = FALSE)
cls <- lab$comparison_label[match(s$metadata$project_cell_id, lab$project_cell_id)]
gs <- match(genes, s$genes); stopifnot(!anyNA(gs))
pos <- cls %in% "POSITIVE"; neg <- cls %in% "NEGATIVE"
res <- rbind(res, data.frame(family_id = SINGH, gene_symbol = genes,
  effect = as.numeric(Matrix::rowMeans(s$normalized[gs, pos, drop = FALSE]) - Matrix::rowMeans(s$normalized[gs, neg, drop = FALSE])),
  det_ec = as.numeric(Matrix::rowMeans(s$counts[gs, pos, drop = FALSE] > 0)),
  det_other = as.numeric(Matrix::rowMeans(s$counts[gs, neg, drop = FALSE] > 0)),
  ec_cells = sum(pos), other_cells = sum(neg), stringsAsFactors = FALSE))
rm(s); invisible(gc(verbose = FALSE))

# control: every effect equals panel D's
k <- match(paste(D$family_id, D$gene_symbol), paste(res$family_id, res$gene_symbol)); stopifnot(!anyNA(k))
stopifnot(nrow(D) == nrow(res), max(abs(D$effect - res$effect[k])) < 1e-8)

# ---- McCauley 2023 HIO grafts --------------------------------------------------------------------
cnt <- as(Matrix::readMM(a[["counts"]]), "CsparseMatrix")
ga <- readLines(a[["genes"]]); cells <- read.delim(a[["cells"]], stringsAsFactors = FALSE)
if (nrow(cnt) != length(ga)) cnt <- Matrix::t(cnt)
stopifnot(nrow(cnt) == length(ga), ncol(cnt) == nrow(cells), !anyDuplicated(ga), all(genes %in% ga))
dimnames(cnt) <- list(ga, cells$figure2_cell_id)
cells$sample <- sub("[.][0-9]+$", "", cells$sample)                 # wt1.1 + wt1.2 -> wt1
stopifnot(setequal(cells$sample, c("wt1", "wt2")), all(cells$author_label %in% c("Endocrine_1", "Endocrine_2")))
hio <- read.delim(a[["hio"]], stringsAsFactors = FALSE)
lib <- hio$gex_total_counts[match(cells$figure2_cell_id, hio$cell_key)]; stopifnot(!anyNA(lib))
z <- cnt[genes, , drop = FALSE] %*% Matrix::Diagonal(x = 1e6 / lib); z@x <- log1p(z@x)
hd <- read.delim(a[["hio_d"]], stringsAsFactors = FALSE)
per <- lapply(c("wt1", "wt2"), function(w) {
  p <- cells$sample == w & cells$author_label == "Endocrine_1"; n <- cells$sample == w & cells$author_label == "Endocrine_2"
  e <- as.numeric(Matrix::rowMeans(z[, p, drop = FALSE]) - Matrix::rowMeans(z[, n, drop = FALSE]))
  h <- hd[hd$sample == w, ]; stopifnot(max(abs(h$effect[match(genes, h$gene_symbol)] - e)) < 1e-8)   # control
  list(e = e, de = as.numeric(Matrix::rowMeans(cnt[genes, p, drop = FALSE] > 0)),
       do = as.numeric(Matrix::rowMeans(cnt[genes, n, drop = FALSE] > 0)), np = sum(p), nn = sum(n))
})
res <- rbind(res, data.frame(family_id = HIO, gene_symbol = genes, effect = (per[[1]]$e + per[[2]]$e) / 2,
  det_ec = (per[[1]]$de + per[[2]]$de) / 2, det_other = (per[[1]]$do + per[[2]]$do) / 2,
  ec_cells = per[[1]]$np + per[[2]]$np, other_cells = per[[1]]$nn + per[[2]]$nn, stringsAsFactors = FALSE))
rownames(res) <- NULL
write.table(res, a[["out"]], sep = "\t", row.names = FALSE, quote = FALSE)
