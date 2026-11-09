# Figure S1C: leave-one-group-out and leave-one-study-out refits of the four-marker EC classifier.
#
# Each row drops one fitting group (or one fitting study), refits the per-gene standardization, the
# Gaussian-mixture cutoff, the candidate set and the neighborhood check on what is left, and
# reclassifies all 1,199 fitting cells. The within-group graphs do not depend on the fit set and are
# built once. Reported per row: the refitted cutoff, candidate and high-confidence counts, Jaccard
# overlap with the full fit, and how many cells moved in each direction.
#
# Inputs : matrix, features, cells (fitting-cell matrix, shared/fitting_matrix/), meta (per-cell raw marker counts and totals)
# Output : out
#   Rscript cutoff_stability.R matrix=... features=... cells=... meta=... out=...
args <- commandArgs(trailingOnly = TRUE)
a <- setNames(sub("^[^=]*=", "", args), sub("=.*$", "", args))
source(file.path(dirname(sub("^--file=", "", grep("^--file=", commandArgs(), value = TRUE))), "ec_classifier_lib.R"))

EC_CORE <- c("TPH1", "SLC18A1", "FEV", "LMX1A")
MIN_DETECTED <- 2L
raw_cols <- paste0("raw_count__", EC_CORE)

cg <- read.delim(gzfile(a[["meta"]]), stringsAsFactors = FALSE, colClasses = "character")[
  , c("capture_qualified_key", "arm_id", "study_family", "model_context", "cell_total_counts", raw_cols)]
cg <- cg[cg$arm_id != "GSE178342", ]
features <- read.delim(a[["features"]], stringsAsFactors = FALSE)
cells <- read.delim(a[["cells"]], stringsAsFactors = FALSE)
n_cells <- nrow(cells); n_feat <- nrow(features)
stopifnot(nrow(cg) == n_cells)
ord <- match(cells$actual_cell_key, cg$capture_qualified_key)
stopifnot(!anyNA(ord), !any(duplicated(ord)))
cg <- cg[ord, ]
con <- file(a[["matrix"]], "rb")
X <- readBin(con, "double", n = as.numeric(n_cells) * n_feat, size = 8, endian = "little")
close(con)
dim(X) <- c(n_feat, n_cells)

# log1p counts per million of the four markers, from the raw counts and each cell's total
total <- as.numeric(cg$cell_total_counts)
counts <- vapply(raw_cols, function(cc) as.numeric(cg[[cc]]), numeric(n_cells))
colnames(counts) <- EC_CORE
L <- log1p(counts * (1e6 / total))
detected <- as.integer(rowSums(counts > 0))

keep <- cg$model_context == "primary_intestine"
cg <- cg[keep, ]; L <- L[keep, , drop = FALSE]; detected <- detected[keep]
X <- X[, keep, drop = FALSE]; invisible(gc(verbose = FALSE))
arm <- cg$arm_id
n_cells <- sum(keep)
stopifnot(n_cells == 1199L)
graphs <- ec_graphs_build(X, arm)

run_fit <- function(fit) {
  mu_g <- colMeans(L[fit, , drop = FALSE])
  sd_g <- sqrt(colMeans(sweep(L[fit, , drop = FALSE], 2, mu_g, "-")^2))
  Z <- sweep(sweep(L, 2, mu_g, "-"), 2, sd_g, "/")
  score <- rowMeans(Z)
  g <- ec_gmm_fit(score[fit])
  cand <- score >= g$threshold & detected >= MIN_DETECTED
  final <- ec_graph_veto(cand, graphs)
  list(threshold = g$threshold, cand = cand, final = final, calls = sum(cand), high = sum(final), gmm = g)
}
jaccard <- function(x, y) { u <- sum(x | y); if (u == 0L) 1 else sum(x & y) / u }

r0 <- run_fit(rep(TRUE, n_cells))
stopifnot(r0$high == 289L)
units <- rbind(data.frame(design = "Leave one arm out", held_out = sort(unique(cg$arm_id)), stringsAsFactors = FALSE),
               data.frame(design = "Leave one family out", held_out = sort(unique(cg$study_family)), stringsAsFactors = FALSE))
rows <- do.call(rbind, lapply(seq_len(nrow(units)), function(i) {
  u <- units[i, ]
  key <- if (u$design == "Leave one arm out") cg$arm_id else cg$study_family
  fit <- key != u$held_out
  r <- run_fit(fit)
  data.frame(design = u$design, held_out = u$held_out, fit_cells = sum(fit), held_cells = sum(!fit),
             threshold = r$threshold, calls = r$calls, jaccard = jaccard(r$cand, r0$cand),
             cand_gained = sum(r$cand & !r0$cand), cand_lost = sum(!r$cand & r0$cand),
             high = r$high, jaccard_high = jaccard(r$final, r0$final),
             jaccard_nonhigh = jaccard(!r$final, !r0$final),
             gained = sum(r$final & !r0$final), lost = sum(!r$final & r0$final),
             veto_only_lost = sum(!r$final & r0$final & r$cand & r0$cand),
             veto_only_gained = sum(r$final & !r0$final & r$cand & r0$cand),
             init_index = r$gmm$init_index, iterations = r$gmm$iterations, stringsAsFactors = FALSE)
}))

# A study with more than one fitting group (Elmentaite 2020: child and fetal) is tagged "A"; a study
# with a single group repeats that group's row and is not drawn.
af <- unique(cg[, c("arm_id", "study_family")])
family_arms <- split(af$arm_id, af$study_family)
multi_arm <- names(family_arms)[lengths(family_arms) > 1L]
rows$drawn <- !(rows$design == "Leave one family out" & !(rows$held_out %in% multi_arm))
TAGS <- setNames(LETTERS[seq_along(multi_arm)], multi_arm)
rows$tag <- mapply(function(design, held_out) {
  if (design == "Leave one family out") return(unname(TAGS[held_out]))
  fam <- unique(cg$study_family[cg$arm_id == held_out])
  if (fam %in% multi_arm) unname(TAGS[fam]) else NA_character_
}, rows$design, rows$held_out)

write.table(rows[order(rows$design, rows$held_out),
                 c("design", "held_out", "held_cells", "fit_cells", "threshold", "calls", "jaccard", "cand_lost",
                   "cand_gained", "high", "jaccard_high", "jaccard_nonhigh", "lost", "gained",
                   "veto_only_lost", "veto_only_gained", "init_index", "iterations", "tag", "drawn")],
            a[["out"]], sep = "\t", row.names = FALSE, quote = FALSE)
