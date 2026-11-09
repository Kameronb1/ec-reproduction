# Figure S3C (four-marker classifier half): the four-marker EC classifier under read thinning,
# like-for-like with the 20-gene classifier's depth test in classifier.R.
#
# Every count of the scored study is kept with probability f (binomial thinning), f = 0.5 and 0.25,
# three draws each. The whole four-marker classifier is re-run on the thinned counts: log1p-CPM, the
# fixed fitting standardization and cutoff, the >= 2 of 4 markers floor, and the within-group graphs rebuilt
# before the local-consistency check. Reported: ROC AUC of the marker score, sensitivity (cells with
# the EC primary study label that are classified as EC) and specificity (other terminal EEC not
# classified as EC).
#
# Inputs : h5ad_dir, fit (JSON from marker_fit.R), scoring (test-study classes from heldout.R)
# Output : out
args <- commandArgs(trailingOnly = TRUE)
a <- setNames(sub("^[^=]*=", "", args), sub("=.*$", "", args))
HERE <- dirname(sub("^--file=", "", grep("^--file=", commandArgs(), value = TRUE)))
source(file.path(HERE, "ec_classifier_lib.R")); source(file.path(HERE, "read_h5ad.R"))
suppressPackageStartupMessages({library(Matrix); library(jsonlite)})

FAMILIES <- c(gut_cell_atlas = "SF_ELMENTAITE2021_GUT_ATLAS", hickey = "SF_HICKEY2023_ADULT_EEC",
              gse238274 = "SF_INTESTINAL_GSE238274", egozi = "SF_EGOZI2021_DEVELOPMENTAL_EEC")
EC_CORE <- c("TPH1", "SLC18A1", "FEV", "LMX1A")
cells <- read.delim(gzfile(a[["scoring"]]), stringsAsFactors = FALSE); stopifnot(nrow(cells) == 6318L)
live <- fromJSON(a[["fit"]], simplifyVector = FALSE)
CUT <- as.numeric(live$ec_core_cutoff); MINM <- as.integer(live$min_markers_detected); VQ <- as.numeric(live$neighbourhood$q_threshold)
std <- do.call(rbind, lapply(live$standardisation, as.data.frame)); rownames(std) <- std$gene
MU <- std[EC_CORE, "mean"]; SD <- std[EC_CORE, "sd"]
binary_auc <- function(y, s) { r <- rank(s); n1 <- sum(y == 1L); n0 <- sum(y == 0L); (sum(r[y == 1L]) - n1 * (n1 + 1) / 2) / (n1 * n0) }
thin <- function(X, f, seed) { set.seed(as.integer(seed)); Y <- X; Y@x <- as.numeric(rbinom(length(X@x), size = X@x, prob = f)); drop0(Y) }

rows <- list()
for (fam in names(FAMILIES)) {
  sce <- read_h5ad_as_sce(file.path(a[["h5ad_dir"]], paste0(fam, ".h5ad")))
  X <- SummarizedExperiment::assay(sce, "X"); rd <- as.data.frame(SummarizedExperiment::rowData(sce), stringsAsFactors = FALSE)
  ft <- intersect(c("feature_type", "original_feature_type"), colnames(rd))
  keep <- if (length(ft)) as.character(rd[[ft[1]]]) == "Gene Expression" else rep(TRUE, nrow(X))
  X <- X[keep, , drop = FALSE]; rd <- rd[keep, , drop = FALSE]
  sym <- if ("gene_symbol" %in% colnames(rd)) as.character(rd$gene_symbol) else rownames(X)
  gi <- match(EC_CORE, sym); stopifnot(!anyNA(gi))
  ff <- cells[cells$family == fam, ]; ff <- ff[match(colnames(X), ff$cell_key), ]
  stopifnot(identical(ff$cell_key, colnames(X)))
  arm <- ff$arm_id
  lab <- ff$author_label_class; ev <- lab %in% c("EXPLICIT_EC", "TERMINAL_NON_EC_EEC"); resp <- as.integer(lab[ev] == "EXPLICIT_EC")
  runs <- rbind(data.frame(fraction = 1, replicate = 0L), expand.grid(fraction = c(.5, .25), replicate = 1:3))
  for (i in seq_len(nrow(runs))) {
    fr <- runs$fraction[i]; rp <- runs$replicate[i]
    Xt <- if (fr == 1) X else thin(X, fr, 20261008L + rp + 100L * match(fam, names(FAMILIES)) + round(1000 * fr))
    tot <- colSums(Xt); stopifnot(all(tot > 0))
    L <- log1p(sweep(Xt, 2, tot, "/") * 1e6)
    score <- colMeans((as.matrix(L[gi, , drop = FALSE]) - MU) / SD)
    det <- as.integer(colSums(Xt[gi, , drop = FALSE] > 0))
    cand <- score >= CUT & det >= MINM
    g <- ec_graphs_build(L, arm, excluded_rows = integer(0))
    stopifnot(!any(vapply(g, is.null, logical(1))))
    call <- if (any(cand)) ec_graph_veto(cand, g, VQ) else cand
    if (fr == 1) stopifnot(identical(unname(call), as.logical(ff$locally_consistent)))   # unthinned = the classes from heldout.R
    cl <- call[ev]
    rows[[length(rows) + 1L]] <- data.frame(study = FAMILIES[[fam]], fraction = fr, replicate = rp,
      roc_auc = binary_auc(resp, score[ev]), sensitivity = mean(cl[resp == 1L]), specificity = mean(!cl[resp == 0L]),
      ec_called = sum(cl[resp == 1L]), ec_n = sum(resp == 1L), other_called = sum(cl[resp == 0L]), other_n = sum(resp == 0L),
      median_total_counts = median(tot))
  }
  rm(sce, X); invisible(gc(verbose = FALSE))
}
write.table(do.call(rbind, rows), a[["out"]], sep = "\t", quote = TRUE, row.names = FALSE)
