# Figure 2C marker dots: for each primary study label of panel C and each of the four classifier markers
# (TPH1, SLC18A1, FEV, LMX1A), the cells with at least one raw read of the marker, their share of the
# label's cells (dot size) and the marker's mean log1p counts per million among those cells (dot color).
#
# What it computes
#   Labels in panel C's order (C_primary_labels.tsv), then the McCauley 2023 HIO grafts' Endocrine_1 and
#   Endocrine_2 (HIO sample groups pooled). Raw counts and log-normalized expression of the markers:
#   Elmentaite 2021, Hickey 2023 and Egozi 2021 from the training object, Singh 2024 from its object
#   (both store log1p(count / library total x 1e6)); McCauley 2023 from the GSE214852 export as
#   log1p(1e6 x count / library total) with each cell's whole-library total. mean_log_detecting is NA
#   where no cell of the label has a read.
#   Controls: per label, the cells and the cells classified as EC equal panel C's; the primary study labels equal those stored
#   in the objects; the stored log expression equals log1p(1e6 x count / library total); per cell, the
#   markers with a read equal the four-marker classifier's markers detected.
#
# Inputs : cells (per-cell table from heldout_marker_calls.R), analysis (training RDS), secondary
#          (Singh 2024 RDS), counts, genes, hio_cells (the GSE214852 endocrine export), hio (four-marker
#          classifier results from four_marker_fit.R), labels (C_primary_labels.tsv), bc (hio_BC_classifier_calls_by_label.tsv)
# Outputs: out (C_marker_detection.tsv)
args <- commandArgs(trailingOnly = TRUE)
a <- setNames(sub("^[^=]*=", "", args), sub("=.*$", "", args))
suppressPackageStartupMessages(library(Matrix))
rd <- function(p) read.delim(p, stringsAsFactors = FALSE, encoding = "UTF-8")
GENES <- c("TPH1", "SLC18A1", "FEV", "LMX1A")

# ---- panel C's labels, with the HIO grafts' two clusters ---------------------------------------
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
for (i in seq_len(nrow(labels_C))) {                        # control: cells and calls equal panel C's
  r <- cells$fam == labels_C$family[i] & cells$label == labels_C$author_label_original[i]
  stopifnot(sum(r) == labels_C$n[i], sum(cells$called[r]) == labels_C$called[i])
}
stopifnot(nrow(cells) == sum(labels_C$n))

# ---- the four markers' raw counts and log expression for every cell ------------------------------
CNT <- EXP <- matrix(NA_real_, nrow(cells), length(GENES), dimnames = list(NULL, GENES))
for (f in c("analysis", "secondary")) {
  x <- readRDS(a[[f]]); gi <- match(GENES, x$genes); stopifnot(!anyNA(gi))
  j <- match(cells$key, x$metadata$project_cell_id); r <- which(!is.na(j) & cells$fam != "hio")
  stopifnot(identical(x$metadata$author_label_original[j[r]], cells$label[r]))      # same primary study labels
  CNT[r, ] <- t(as.matrix(x$counts[gi, j[r], drop = FALSE])); EXP[r, ] <- t(as.matrix(x$normalized[gi, j[r], drop = FALSE]))
  rm(x); invisible(gc(verbose = FALSE))
}
cnt <- as(Matrix::readMM(a[["counts"]]), "CsparseMatrix")
ga <- readLines(a[["genes"]])
if (nrow(cnt) != length(ga)) cnt <- Matrix::t(cnt)
stopifnot(nrow(cnt) == length(ga), ncol(cnt) == nrow(hc), all(GENES %in% ga))
r <- which(cells$fam == "hio")
CNT[r, ] <- t(as.matrix(cnt[match(GENES, ga), , drop = FALSE]))
EXP[r, ] <- log1p(1e6 * CNT[r, ] / cells$lib[r])
stopifnot(!anyNA(CNT), !anyNA(EXP), !anyNA(cells$lib))
stopifnot(max(abs(EXP - log1p(1e6 * CNT / cells$lib))) < 1e-8)   # control: one log scale for every study
stopifnot(all(rowSums(CNT > 0) == cells$markers))                # control: markers with a read = markers detected

# ---- per label and marker ------------------------------------------------------------------------
mm <- do.call(rbind, lapply(seq_len(nrow(labels_C)), function(i) {
  r <- cells$fam == labels_C$family[i] & cells$label == labels_C$author_label_original[i]
  data.frame(family = labels_C$family[i], author_label_original = labels_C$author_label_original[i], n = sum(r),
             called = sum(cells$called[r]), gene = GENES, cells_detecting = colSums(CNT[r, , drop = FALSE] > 0),
             share_detecting = colMeans(CNT[r, , drop = FALSE] > 0),
             mean_log_detecting = vapply(GENES, function(g) { k <- CNT[r, g] > 0; if (any(k)) mean(EXP[r, g][k]) else NA_real_ }, 0),
             stringsAsFactors = FALSE)
}))
rownames(mm) <- NULL
write.table(mm, a[["out"]], sep = "\t", row.names = FALSE, quote = FALSE)
