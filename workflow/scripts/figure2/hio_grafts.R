# Figure 2B-E, McCauley 2023 HIO grafts (GSE214852): the HIO rows of every panel.
#
# Primary study labels: the authors' Endocrine_1 cluster is taken as EC and Endocrine_2 as other
# terminal EEC. Units are the two project-defined HIO sample groups wt1 (libraries wt1.1 and wt1.2
# combined) and wt2; they are groups of cells from one cell line, not independent organoids.
#
# What it computes
#   B/C  cells the four-marker classifier classifies as EC (from four_marker_fit.R) within each primary
#        study label, per HIO sample group.
#   D    for the 20 panel-D genes: mean log-normalized expression in Endocrine_1 minus Endocrine_2,
#        over all cells and per HIO sample group. Expression is log1p(counts / library total x 1e6)
#        with the whole-library total of each cell, the scale of the training object.
#   E    the label-trained classifier with the top 10, 20 and 30 ranked genes, weights and threshold
#        fitted on the three primary-label training studies (as for Singh 2024), >= 25 % of genes
#        detected, and a neighborhood gate on a graph built from these cells (excluding the four
#        markers and the top 30 genes); cells classified as EC by this classifier and by the four-marker
#        classifier within the primary study labels. The 20-gene classifier's per-cell score, threshold
#        and result are also written (scores), for e_cells.R (panel E).
#
# Inputs : counts (Matrix Market, genes x cells), genes, cells (the GSE214852 endocrine export),
#          hio (four-marker classifier results from four_marker_fit.R), analysis (training RDS), ranking (gene_ranking.tsv)
# Outputs: bc (BC_rule_calls_by_label.tsv), d (D_hio_effects.tsv), e (E_hio_summary.tsv),
#          scores (per-cell 20-gene classifier results)
args <- commandArgs(trailingOnly = TRUE)
a <- setNames(sub("^[^=]*=", "", args), sub("=.*$", "", args))
source(file.path(dirname(sub("^--file=", "", grep("^--file=", commandArgs(), value = TRUE))), "classifier_functions.R"))
wt <- function(x, path) write.table(x, path, sep = "\t", quote = TRUE, row.names = FALSE)
PRIM <- c("SF_EGOZI2021_DEVELOPMENTAL_EEC", "SF_ELMENTAITE2021_GUT_ATLAS", "SF_HICKEY2023_ADULT_EEC")

x <- readRDS(a[["analysis"]]); m <- x$metadata; labels <- training_labels(m)
ranking <- read.delim(a[["ranking"]], stringsAsFactors = FALSE)
TOP <- head(ranking$gene_symbol, 30L); D_GENES <- head(TOP, 20L)

cnt <- as(Matrix::readMM(a[["counts"]]), "CsparseMatrix")
genes_all <- readLines(a[["genes"]]); cells <- read.delim(a[["cells"]], stringsAsFactors = FALSE)
if (nrow(cnt) != length(genes_all)) cnt <- Matrix::t(cnt)
stopifnot(nrow(cnt) == length(genes_all), ncol(cnt) == nrow(cells), ncol(cnt) == 598L, !anyDuplicated(genes_all))
dimnames(cnt) <- list(genes_all, cells$figure2_cell_id)
cells$sample <- sub("[.][0-9]+$", "", cells$sample)            # wt1.1 + wt1.2 -> wt1
stopifnot(setequal(cells$sample, c("wt1", "wt2")))
hio <- read.delim(a[["hio"]], stringsAsFactors = FALSE)
hio <- hio[match(cells$figure2_cell_id, hio$cell_key), ]
stopifnot(!anyNA(hio$cell_key), all(hio$arm_id == "GSE214852_WT"),
          all(hio$presentation_class %in% c("HIGH_CONFIDENCE_EC", "INSUFFICIENT_EC_EVIDENCE")))
rule_label <- ifelse(hio$presentation_class == "HIGH_CONFIDENCE_EC", "POSITIVE", "NEGATIVE")
stopifnot(all(cells$author_label %in% c("Endocrine_1", "Endocrine_2")))
label <- ifelse(cells$author_label == "Endocrine_1", "POSITIVE", "NEGATIVE")

# ---- B/C ----------------------------------------------------------------------------------------
bc <- aggregate(cbind(n = 1, rule_called = rule_label == "POSITIVE") ~ sample + author_label,
                data = data.frame(sample = cells$sample, author_label = cells$author_label, rule_label), FUN = sum)
wt(bc, a[["bc"]])

# ---- expression on the training scale ----------------------------------------------------------
lib <- hio$gex_total_counts                                    # whole-library total per cell
stopifnot(all(lib >= as.numeric(Matrix::colSums(cnt))))
keep_genes <- intersect(x$genes, genes_all); stopifnot(all(TOP %in% keep_genes))
z <- cnt[keep_genes, , drop = FALSE] %*% Matrix::Diagonal(x = 1e6 / lib); z@x <- log1p(z@x)
dimnames(z) <- list(keep_genes, cells$figure2_cell_id)

# ---- D ----------------------------------------------------------------------------------------
pos <- label == "POSITIVE"; neg <- label == "NEGATIVE"
d_row <- function(s, k) data.frame(family_id = "SF_GSE214852_HIO_GRAFT", sample = s, gene_symbol = D_GENES,
  effect = as.numeric(Matrix::rowMeans(z[D_GENES, k & pos]) - Matrix::rowMeans(z[D_GENES, k & neg])),
  ec_cells = sum(k & pos), other_cells = sum(k & neg), label_source = "AUTHOR", stringsAsFactors = FALSE)
D <- rbind(d_row("all", rep(TRUE, ncol(z))), do.call(rbind, lapply(sort(unique(cells$sample)), function(s) d_row(s, cells$sample == s))))
D$panel_rank <- match(D$gene_symbol, D_GENES)
wt(D, a[["d"]])

# ---- E ----------------------------------------------------------------------------------------
graph <- family_graph(z, unique(c(SELECTION_EXCLUSIONS, TOP)))
scores20 <- NULL
summ <- do.call(rbind, lapply(c(10L, 20L, 30L), function(k) {
  model <- fit_panel_model(x$normalized, m, labels, m$family_id %in% PRIM, head(TOP, k))
  sc <- score_panel_model(z, model)
  call <- gated_call(sc$score, sc$detected_genes, model$detection_requirement, model$threshold, graph)$call
  if (k == 20L) scores20 <<- data.frame(cell_key = cells$figure2_cell_id, score = sc$score, threshold = model$threshold,
                                        call = call, stringsAsFactors = FALSE)
  mr <- rule_label == "POSITIVE"
  data.frame(panel_size = k, ec_n = sum(pos), ec_called = sum(call[pos]), ec_rule_called = sum(mr[pos]),
             other_n = sum(!pos), other_called = sum(call[!pos]), other_rule_called = sum(mr[!pos]),
             roc_auc = binary_auc(as.integer(pos), sc$score), threshold = model$threshold,
             below_detection = sum(sc$detected_genes < model$detection_requirement), stringsAsFactors = FALSE)
}))
rownames(summ) <- NULL
wt(summ, a[["e"]])
write.table(scores20, a[["scores"]], sep = "\t", quote = FALSE, row.names = FALSE)
print(summ)
