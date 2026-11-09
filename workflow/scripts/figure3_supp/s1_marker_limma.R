# Figure S15A: the 20 genes that best separate SC-SIC from the other SC endocrine cells (beta, alpha, delta, epsilon,
# gamma), ranked exactly as Figure 2's classifier genes: genes detected (>= 5% of cells) in at least 3 of the 6 SC-islet
# studies, the four EC marker genes TPH1, SLC18A1, FEV and LMX1A not ranked; sample-by-class mean expression
# (ln(1 + CPM), >= 5 cells per sample and class); limma (~ 0 + sample + SC-SIC, each sample weighted 1 / samples in its
# study, eBayes trend); ranked by moderated t. The heatmap value is, per study, the mean over its samples of SC-SIC minus
# other endocrine mean log expression; sic_detection is, per study, the mean over its samples of the share of SC-SICs
# with >= 1 read for the gene (the S15A dot size). Primary study labels.
# Inputs:  lib (s1_lib.R); the SC endocrine count package: meta, genes, indptr, indices, data
# Outputs: top20 (top20.tsv), effects (heatmap_effects.tsv)
a <- commandArgs(TRUE); A <- setNames(sub("^[^=]*=", "", a), sub("=.*", "", a))
source(A[["lib"]])
TOP_N <- 20
sc <- read_sc_endocrine(A[["meta"]], A[["genes"]], A[["indptr"]], A[["indices"]], A[["data"]])
expr <- sc$expr; genes <- sc$genes; meta <- sc$meta; nC <- nrow(meta)
stopifnot(length(unique(meta$family_id)) == 6)

pos <- meta$target %in% "SC_SIC"; neg <- meta$target %in% c("SC_BETA", "SC_ALPHA", "SC_DELTA", "SC_EPSILON", "SC_GAMMA")
stopifnot(all(pos | neg))
labels <- data.frame(project_cell_id = meta$project_cell_id, label = ifelse(pos, "POSITIVE", "NEGATIVE"), stringsAsFactors = FALSE)

s <- rank_genes(expr, genes, meta, labels, rep(TRUE, nC))
rk <- s$ranking; rk$rank <- seq_len(nrow(rk))
top <- head(rk, TOP_N)
stopifnot(all(top$limma_t > 0))
um <- s$unit_class$metadata
stopifnot(length(unique(paste(um$family_id, um$biological_unit_id))) == 10)   # ten samples, both classes in each

studies <- sort(unique(um$family_id))
eff <- do.call(rbind, lapply(studies, function(f) {
  pp <- um$family_id == f & um$label == "POSITIVE"; nn <- um$family_id == f & um$label == "NEGATIVE"
  data.frame(family_id = f, gene_symbol = top$gene_symbol,
             effect = rowMeans(s$unit_class$values[top$gene_symbol, pp, drop = FALSE]) -
                      rowMeans(s$unit_class$values[top$gene_symbol, nn, drop = FALSE]))
}))
# dot size: SC-SIC detection, the share of a sample's SC-SICs with >= 1 read, averaged over the study's samples as above
dt <- cp_unit_class_means((expr[top$gene_symbol, , drop = FALSE] > 0) * 1, meta, labels, rep(TRUE, nC), top$gene_symbol,
                          PARAMETERS$minimum_cells_per_unit_class)
stopifnot(identical(dt$metadata, um))
eff$sic_detection <- unlist(lapply(studies, function(f)
  rowMeans(dt$values[top$gene_symbol, dt$metadata$family_id == f & dt$metadata$label == "POSITIVE", drop = FALSE])))
eff$panel_rank <- match(eff$gene_symbol, top$gene_symbol)
studies_higher <- tapply(eff$effect > 0, eff$gene_symbol, sum)[top$gene_symbol]

wt <- function(d, f) write.table(d, f, sep = "\t", quote = FALSE, row.names = FALSE, na = "")
wt(data.frame(rank = top$rank, gene = top$gene_symbol, limma_logFC = top$limma_logfc, limma_t = top$limma_t,
              limma_fdr = top$limma_fdr, studies_higher = as.integer(studies_higher),
              median_study_effect = top$median_family_effect, median_detection = top$median_detection_fraction), A[["top20"]])
wt(eff, A[["effects"]])
