# Figure S2B: the label-free four-marker EC classifier applied unchanged to the three intestinal test studies,
# which have published labels and were not used to fit it (Elmentaite 2021, Hickey 2023, Egozi 2021).
#
# What it computes, per study
#   1. From raw counts (Gene Expression features only): log1p counts-per-million, the four-marker
#      score with the fitted standardization constants, and markers detected.
#   2. Candidate = score >= fitted cutoff and >= 2 of 4 markers detected; classified as EC =
#      candidate that passes the within-group (per deposited dataset) 25-nearest-neighbor
#      local-consistency test on all genes (BH q <= 0.05).
#   3. Cells, and cells classified as EC, within every published label (labels).
#   4. Every cell's study, published label, library total, markers detected, score and result
#      (cells), read by c_marker_detection.R (the marker dots of S2B).
#
# Inputs : fit (fitted_constants.json from marker_fit.R), h5ad_<study> for gut_cell_atlas, hickey, egozi
# Outputs: labels (C_primary_labels.tsv), cells (per-cell table)
#   Rscript heldout_marker_calls.R fit=... h5ad_gut_cell_atlas=... h5ad_hickey=... h5ad_egozi=... labels=... cells=...
args <- commandArgs(trailingOnly = TRUE)
a <- setNames(sub("^[^=]*=", "", args), sub("=.*$", "", args))
suppressPackageStartupMessages({library(Matrix); library(rhdf5); library(jsonlite)})
source(file.path(dirname(sub("^--file=", "", grep("^--file=", commandArgs(), value = TRUE))), "four_marker_lib.R"))

FAMILIES <- c("gut_cell_atlas", "hickey", "egozi")
STUDY <- c(gut_cell_atlas = "Elmentaite 2021 · fetal + adult intestine",
           hickey = "Hickey 2023 · adult intestine",
           egozi = "Egozi 2021 · fetal + neonatal intestine")
EXPECTED_CELLS <- c(gut_cell_atlas = 897L, hickey = 3533L, egozi = 784L)
CLASS_ORDER <- c("EXPLICIT_EC", "TERMINAL_NON_EC_EEC", "BROAD_EEC_PARENT_UNRESOLVED", "ENDOCRINE_PROGENITOR")

fit <- fromJSON(a[["fit"]])
EC_CORE <- fit$ec_core_genes                                    # the fit's marker genes (TPH1, SLC18A1, FEV, LMX1A)
STANDARDIZED_GENES <- c(EC_CORE, "DDC", "TAC1")
mu <- setNames(fit$standardisation$mean, fit$standardisation$gene)[STANDARDIZED_GENES]
sdv <- setNames(fit$standardisation$sd, fit$standardisation$gene)[STANDARDIZED_GENES]

# Minimal h5ad reader: a column of /obs or /var (plain or categorical), and X (CSR, cells x genes).
h5_col <- function(f, group, k) {
  nm <- paste0("/", group, "/", k); ls <- h5ls(f)
  if (any(ls$group == nm & ls$name == "categories")) {
    cats <- as.character(h5read(f, paste0(nm, "/categories"))); cods <- as.integer(h5read(f, paste0(nm, "/codes")))
    stopifnot(all(cods >= 0L)); cats[cods + 1L]
  } else as.character(h5read(f, nm))
}
h5_counts <- function(f, n_obs, n_var) {
  new("dgCMatrix", i = as.integer(h5read(f, "/X/indices")), p = as.integer(h5read(f, "/X/indptr")),
      x = as.numeric(h5read(f, "/X/data")), Dim = c(as.integer(n_var), as.integer(n_obs)))
}

lab_rows <- list(); cell_rows <- list()
for (fam in FAMILIES) {
  f <- a[[paste0("h5ad_", fam)]]
  cell <- h5_col(f, "obs", "project_cell_id")
  # gene symbol: the gene_symbol column, else the var index; feature type: feature_type or
  # original_feature_type, else every feature counts as Gene Expression
  vcols <- h5ls(f)$name[h5ls(f)$group == "/var"]
  sym <- h5_col(f, "var", if ("gene_symbol" %in% vcols) "gene_symbol" else h5readAttributes(f, "var")[["_index"]])
  tcol <- intersect(c("feature_type", "original_feature_type"), vcols)
  ftype <- if (length(tcol)) h5_col(f, "var", tcol[1]) else rep("Gene Expression", length(sym))
  X <- h5_counts(f, length(cell), length(sym))
  obs <- data.frame(arm = h5_col(f, "obs", "object_dataset_id"),
                    orig = h5_col(f, "obs", "author_label_original"), cls = h5_col(f, "obs", "author_label_class"),
                    stringsAsFactors = FALSE)
  h5closeAll()
  stopifnot(ncol(X) == EXPECTED_CELLS[[fam]], all(X@x == round(X@x)))         # raw integer counts
  keep <- ftype == "Gene Expression"
  X <- X[keep, , drop = FALSE]; sym <- sym[keep]
  total <- Matrix::colSums(X); stopifnot(all(total > 0))
  gr <- match(STANDARDIZED_GENES, sym); stopifnot(!anyNA(gr), !any(duplicated(sym[sym %in% STANDARDIZED_GENES])))
  counts6 <- as.matrix(X[gr, , drop = FALSE]); rownames(counts6) <- STANDARDIZED_GENES
  detected <- as.integer(Matrix::colSums(X[gr[seq_along(EC_CORE)], , drop = FALSE] > 0))
  L6 <- log1p(sweep(counts6, 2, total, "/") * 1e6)
  score <- colMeans(((L6 - mu[rownames(L6)]) / sdv[rownames(L6)])[EC_CORE, , drop = FALSE])
  cand <- (score >= fit$ec_core_cutoff) & (detected >= fit$min_markers_detected)
  full_lcpm <- log1p(sweep(X, 2, total, "/") * 1e6)        # the graph reads every gene
  high <- locally_consistent(cand, full_lcpm, obs$arm)
  rm(full_lcpm, X); invisible(gc(verbose = FALSE))
  cat(sprintf("%-15s %4d cells, %4d candidates, %4d classified as EC\n", fam, length(cand), sum(cand), sum(high)))
  cell_rows[[fam]] <- data.frame(family = fam, project_cell_id = cell, author_label_original = obs$orig,
                                 author_label_class = obs$cls, gex_total_counts = as.numeric(total),
                                 markers_detected = detected, ec_core_score = unname(score), called = unname(high),
                                 stringsAsFactors = FALSE)
  # every published label
  agg <- aggregate(list(n = rep(1L, length(high)), called = as.integer(high)),
                   by = list(author_label_original = obs$orig, author_label_class = obs$cls), FUN = sum)
  agg$family <- fam
  lab_rows[[fam]] <- agg
}

fmt <- function(v) format(v, trim = TRUE, big.mark = ",")
C <- do.call(rbind, lab_rows); C$rate <- 100 * C$called / C$n
C <- C[order(match(C$family, FAMILIES), match(C$author_label_class, CLASS_ORDER), -C$rate), ]
C$study <- unname(STUDY[C$family])
C$label <- sprintf("%s/%s", fmt(C$called), fmt(C$n))
C <- C[c("family", "study", "author_label_original", "author_label_class", "n", "called", "rate", "label")]
write.table(C, a[["labels"]], sep = "\t", row.names = FALSE, quote = TRUE, fileEncoding = "UTF-8")
write.table(do.call(rbind, cell_rows), a[["cells"]], sep = "\t", row.names = FALSE, quote = TRUE, fileEncoding = "UTF-8")
