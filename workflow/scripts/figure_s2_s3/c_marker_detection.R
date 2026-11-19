# Figure S2B marker dots: for each published label of S2B and each of the four classifier markers
# (TPH1, SLC18A1, FEV, LMX1A), the cells with at least one raw read of the marker, their share of the
# label's cells (dot size) and the marker's mean log1p counts per million among those cells (dot color).
#
# What it computes
#   Labels in the order of C_primary_labels.tsv, then the McCauley 2023 HIO grafts' Endocrine_1 and
#   Endocrine_2 (HIO sample groups pooled). Raw counts of the markers: Elmentaite 2021, Hickey 2023 and
#   Egozi 2021 from their deposited objects, McCauley 2023 from the GSE214852 export. Log expression is
#   log1p(count / library total x 1e6) with each cell's Gene Expression library total (test studies) or
#   whole-library total (HIO grafts). mean_log_detecting is NA where no cell of the label has a read.
#   Controls: per label, the cells and the cells classified as EC equal those of the labels table; the
#   published labels equal those stored in the objects; per cell, the markers with a read equal the
#   four-marker classifier's markers detected.
#
# Inputs : cells (per-cell table from heldout_marker_calls.R), h5ad_<study> for gut_cell_atlas, hickey, egozi,
#          counts, genes, hio_cells (the GSE214852 endocrine export), hio (HIO classes from four_marker_fit.R),
#          labels (C_primary_labels.tsv), bc (hio_BC_classifier_calls_by_label.tsv)
# Outputs: out (C_marker_detection.tsv)
args <- commandArgs(trailingOnly = TRUE)
a <- setNames(sub("^[^=]*=", "", args), sub("=.*$", "", args))
suppressPackageStartupMessages({library(Matrix); library(rhdf5)})
rd <- function(p) read.delim(p, stringsAsFactors = FALSE, encoding = "UTF-8")
GENES <- c("TPH1", "SLC18A1", "FEV", "LMX1A")
FAMILIES <- c("gut_cell_atlas", "hickey", "egozi")

# ---- the labels, with the HIO grafts' two clusters ----------------------------------------------
labels_C <- rd(a[["labels"]])[c("family", "author_label_original", "n", "called")]
bc <- aggregate(cbind(n, rule_called) ~ author_label, rd(a[["bc"]]), sum)
bc <- bc[order(bc$author_label), ]
labels_C <- rbind(labels_C, data.frame(family = "hio", author_label_original = bc$author_label, n = bc$n,
                                       called = bc$rule_called, stringsAsFactors = FALSE))

# ---- cells: label, markers detected, call, library total ---------------------------------------
cl <- rd(a[["cells"]])
hc <- rd(a[["hio_cells"]]); hm <- rd(a[["hio"]])
k <- match(hc$figure2_cell_id, hm$cell_key); stopifnot(!anyNA(k))
cells <- rbind(data.frame(key = cl$project_cell_id, fam = cl$family, label = cl$author_label_original,
                          markers = cl$markers_detected, called = as.logical(cl$called), lib = cl$gex_total_counts),
               data.frame(key = hc$figure2_cell_id, fam = "hio", label = hc$author_label, markers = hm$markers_detected[k],
                          called = hm$presentation_class[k] == "HIGH_CONFIDENCE_EC", lib = hm$gex_total_counts[k]))
for (i in seq_len(nrow(labels_C))) {                        # control: cells and calls equal the labels table's
  r <- cells$fam == labels_C$family[i] & cells$label == labels_C$author_label_original[i]
  stopifnot(sum(r) == labels_C$n[i], sum(cells$called[r]) == labels_C$called[i])
}
stopifnot(nrow(cells) == sum(labels_C$n))

# ---- the four markers' raw counts for every cell -----------------------------------------------
h5_col <- function(f, group, k) {
  nm <- paste0("/", group, "/", k); ls <- h5ls(f)
  if (any(ls$group == nm & ls$name == "categories")) {
    cats <- as.character(h5read(f, paste0(nm, "/categories"))); cods <- as.integer(h5read(f, paste0(nm, "/codes")))
    stopifnot(all(cods >= 0L)); cats[cods + 1L]
  } else as.character(h5read(f, nm))
}
CNT <- matrix(NA_real_, nrow(cells), length(GENES), dimnames = list(NULL, GENES))
for (fam in FAMILIES) {
  f <- a[[paste0("h5ad_", fam)]]
  id <- h5_col(f, "obs", "project_cell_id"); orig <- h5_col(f, "obs", "author_label_original")
  vcols <- h5ls(f)$name[h5ls(f)$group == "/var"]
  sym <- h5_col(f, "var", if ("gene_symbol" %in% vcols) "gene_symbol" else h5readAttributes(f, "var")[["_index"]])
  tcol <- intersect(c("feature_type", "original_feature_type"), vcols)
  ftype <- if (length(tcol)) h5_col(f, "var", tcol[1]) else rep("Gene Expression", length(sym))
  X <- new("dgCMatrix", i = as.integer(h5read(f, "/X/indices")), p = as.integer(h5read(f, "/X/indptr")),
           x = as.numeric(h5read(f, "/X/data")), Dim = c(length(sym), length(id)))
  h5closeAll()
  gi <- which(ftype == "Gene Expression")[match(GENES, sym[ftype == "Gene Expression"])]; stopifnot(!anyNA(gi))
  j <- match(cells$key, id); r <- which(!is.na(j) & cells$fam == fam)
  stopifnot(length(r) == sum(cells$fam == fam), identical(orig[j[r]], cells$label[r]))   # same published labels
  CNT[r, ] <- t(as.matrix(X[gi, j[r], drop = FALSE]))
  rm(X); invisible(gc(verbose = FALSE))
}
cnt <- as(Matrix::readMM(a[["counts"]]), "CsparseMatrix")
ga <- readLines(a[["genes"]])
if (nrow(cnt) != length(ga)) cnt <- Matrix::t(cnt)
stopifnot(nrow(cnt) == length(ga), ncol(cnt) == nrow(hc), all(GENES %in% ga))
r <- which(cells$fam == "hio")
CNT[r, ] <- t(as.matrix(cnt[match(GENES, ga), , drop = FALSE]))
stopifnot(!anyNA(CNT), !anyNA(cells$lib))
EXP <- log1p(CNT / cells$lib * 1e6)
EXP[r, ] <- log1p(1e6 * CNT[r, ] / cells$lib[r])
stopifnot(all(rowSums(CNT > 0) == cells$markers))                # control: markers with a read = markers detected

# ---- per label and marker ------------------------------------------------------------------------
mm <- do.call(rbind, lapply(seq_len(nrow(labels_C)), function(i) {
  r <- cells$fam == labels_C$family[i] & cells$label == labels_C$author_label_original[i]
  data.frame(family = labels_C$family[i], author_label_original = labels_C$author_label_original[i], n = sum(r),
             gene = GENES, cells_detecting = colSums(CNT[r, , drop = FALSE] > 0),
             share_detecting = colMeans(CNT[r, , drop = FALSE] > 0),
             mean_log_detecting = vapply(GENES, function(g) { k <- CNT[r, g] > 0; if (any(k)) mean(EXP[r, g][k]) else NA_real_ }, 0),
             stringsAsFactors = FALSE)
}))
rownames(mm) <- NULL
write.table(mm, a[["out"]], sep = "\t", row.names = FALSE, quote = FALSE)
