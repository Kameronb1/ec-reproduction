# Figure 2D (gene-ranking datasets: the 20 genes and their per-study effects) and the table behind
# Figure 2E (label-trained 20-gene classifier and four-marker classifier calls within primary study
# labels; the drawn classifier rows come from primary_label_classifier.R's calls).
#
# What it computes
#   D  Ranks every eligible gene on all seven training datasets (limma, see classifier_functions.R);
#      the top 20 must equal the 20-gene panel of the deposited leave-one-study-out computation. Effect per study = mean
#      over donors/samples (>= 5 cells per donor-class) of EC minus other mean log-normalized
#      expression under the training labels; Singh 2024 (not in training) by its primary study
#      labels, over cells. limma_fdr is the gene's FDR in the all-training ranking.
#   E  For the four studies with primary study labels: the cells the deposited computation's 20-gene
#      classifier classifies as EC, among EC-labeled and other-terminal-EEC-labeled cells (strict =
#      classified as EC at the unrelaxed threshold; extra = added by a relaxed threshold, which the
#      deposited computation does not use, so it must be zero), next to the cells the four-marker
#      classifier classifies as EC (from panel C).
#
# Inputs : analysis (training object, RDS), secondary (Singh 2024 object, RDS), and from the
#          deposited computation: run_metadata, nested_predictions, sv_predictions, sv_labels,
#          sv_panels; labels (C_primary_labels.tsv from heldout_marker_calls.R)
# Outputs: heatmap (D_effect_heatmap.tsv), concordance (E_concordance.tsv), ranking (gene_ranking.tsv)
args <- commandArgs(trailingOnly = TRUE)
a <- setNames(sub("^[^=]*=", "", args), sub("=.*$", "", args))
source(file.path(dirname(sub("^--file=", "", grep("^--file=", commandArgs(), value = TRUE))), "classifier_functions.R"))
gz <- function(p) read.delim(gzfile(p), stringsAsFactors = FALSE, check.names = FALSE)
lg <- function(v) as.logical(v)
SINGH <- "SF_INTESTINAL_GSE238274"
AUTHOR_STUDIES <- c("SF_EGOZI2021_DEVELOPMENTAL_EEC", "SF_ELMENTAITE2021_GUT_ATLAS", "SF_HICKEY2023_ADULT_EEC")

# ---- D ----------------------------------------------------------------------------------------
x <- readRDS(a[["analysis"]])
stopifnot(identical(colnames(x$normalized), x$metadata$project_cell_id), !any(x$metadata$family_id == SINGH))
m <- x$metadata; labels <- training_labels(m)
ranking <- rank_training_genes(x$normalized, x$genes, m, labels, rep(TRUE, nrow(m)))
genes <- head(ranking$gene_symbol, 20L)
span <- read.delim(a[["sv_panels"]], stringsAsFactors = FALSE)
run_genes <- span$gene_symbol[span$panel_size == 20][order(span$selected_rank[span$panel_size == 20])]
if (!identical(genes, run_genes)) stop("this ranking does not reproduce the deposited 20 genes: ",paste(setdiff(genes, run_genes), collapse = ","))
write.table(data.frame(rank = seq_len(nrow(ranking)), ranking), a[["ranking"]], sep = "\t", row.names = FALSE, quote = FALSE)

u <- unit_class_means(x$normalized, m, labels, rep(TRUE, nrow(m)), genes, MIN_CELLS_PER_UNIT_CLASS)
um <- u$metadata
eff <- do.call(rbind, lapply(sort(unique(um$family_id)), function(f) {
  p <- um$family_id == f & um$label == "POSITIVE"; n <- um$family_id == f & um$label == "NEGATIVE"; stopifnot(any(p), any(n))
  data.frame(family_id = f, gene_symbol = genes,
             effect = rowMeans(u$values[, p, drop = FALSE]) - rowMeans(u$values[, n, drop = FALSE]), stringsAsFactors = FALSE)
}))
author_families <- sort(unique(m$family_id[labels$reason %in% c("ORIGINAL_AUTHOR_EC", "ORIGINAL_AUTHOR_TERMINAL_NON_EC_EEC")]))
stopifnot(identical(author_families, AUTHOR_STUDIES))
eff$label_source <- ifelse(eff$family_id %in% author_families, "AUTHOR", "MARKER_RULE")
s <- readRDS(a[["secondary"]])
sv_lab <- gz(a[["sv_labels"]])
cls <- sv_lab$comparison_label[match(s$metadata$project_cell_id, sv_lab$project_cell_id)]
gi <- match(genes, s$genes); stopifnot(!anyNA(gi))
held <- data.frame(family_id = SINGH, gene_symbol = genes,
  effect = as.numeric(Matrix::rowMeans(s$normalized[gi, cls %in% "POSITIVE", drop = FALSE]) -
                      Matrix::rowMeans(s$normalized[gi, cls %in% "NEGATIVE", drop = FALSE])),
  label_source = "AUTHOR", stringsAsFactors = FALSE)
D <- rbind(eff, held); D$panel_rank <- match(D$gene_symbol, genes)
D$limma_fdr <- ranking$limma_fdr[match(D$gene_symbol, ranking$gene_symbol)]
rownames(D) <- NULL
write.table(D, a[["heatmap"]], sep = "\t", row.names = FALSE, quote = FALSE)

# ---- E ----------------------------------------------------------------------------------------
# The deposited 20-gene classifier results on cells with primary study labels. two = classified as EC
# (score and neighborhood check); unrelaxed = at or above the unrelaxed threshold (column v05_call).
meta <- gz(a[["run_metadata"]])
meta <- meta[meta$author_label_evidence_layer == "ORIGINAL_AUTHOR_OBJECT" & meta$author_group %in% c("EXPLICIT_EC", "TERMINAL_NON_EC_EEC"), ]
meta$class <- ifelse(meta$author_group == "EXPLICIT_EC", "POSITIVE", "NEGATIVE")
np <- gz(a[["nested_predictions"]])
np <- np[np$panel_id == "NESTED_SIZE_20" & np$outer_family %in% AUTHOR_STUDIES, ]
np$two <- lg(np$stage_one_call) & lg(np$gate_pass); np$unrelaxed <- lg(np$v05_call)
np <- merge(np, meta[c("project_cell_id", "family_id", "class")], by = "project_cell_id")
stopifnot(all(np$family_id == np$outer_family)); np$family <- np$outer_family
sl <- sv_lab[sv_lab$comparison_label %in% c("POSITIVE", "NEGATIVE"), ]
sp <- gz(a[["sv_predictions"]]); sp <- sp[sp$panel_size == 20, ]
sp$two <- lg(sp$two_stage_call); sp$unrelaxed <- lg(sp$v05_call)
sp <- merge(sp, data.frame(project_cell_id = sl$project_cell_id, class = sl$comparison_label), by = "project_cell_id")
sp$family <- SINGH
calls <- rbind(np[c("family", "class", "two", "unrelaxed")], sp[c("family", "class", "two", "unrelaxed")])
r <- do.call(rbind, lapply(split(calls, list(calls$family, calls$class), drop = TRUE), function(d)
  data.frame(family = d$family[1], class = d$class[1], n = nrow(d), strict = sum(d$unrelaxed & d$two),
             extra = sum(d$two & !d$unrelaxed), total = sum(d$two), gate_dropped = sum(d$unrelaxed & !d$two))))
stopifnot(all(r$extra == 0L))                 # no relaxed threshold: every EC cell meets the unrelaxed threshold

# The cells the four-marker classifier classifies as EC (panel C's EC and other-terminal-EEC rows).
C <- read.delim(a[["labels"]], stringsAsFactors = FALSE, encoding = "UTF-8")
fam <- c(gut_cell_atlas = "SF_ELMENTAITE2021_GUT_ATLAS", hickey = "SF_HICKEY2023_ADULT_EEC",
         egozi = "SF_EGOZI2021_DEVELOPMENTAL_EEC", gse238274 = SINGH)
EY <- 0.230
E <- do.call(rbind, lapply(seq_len(nrow(r)), function(i) {
  q <- r[i, ]; short <- names(fam)[fam == q$family]
  ref <- C[C$family == short & C$author_label_class == if (q$class == "POSITIVE") "EXPLICIT_EC" else "TERMINAL_NON_EC_EEC", ]
  stopifnot(sum(ref$n) == q$n)
  mk <- sum(ref$called)
  data.frame(panel_size = 20, family = q$family, label_source = "ORIGINAL_AUTHOR_OBJECT", class = q$class,
             rule = c("panel", "marker", "panel"), seg = c("extra", "marker", "strict"),
             called = c(q$extra, mk, q$strict), rule_total = c(q$total, mk, q$total), n = q$n,
             gate_dropped = q$gate_dropped, marker_rule_fitted_here = FALSE, fam = short,
             study_y = 5 - match(short, names(fam)), stringsAsFactors = FALSE)
}))
E$pct <- 100 * E$called / E$n
E$y <- E$study_y + ifelse(E$rule == "marker", EY, -EY)
E$x0 <- ifelse(E$seg == "extra", 100 * (E$rule_total - E$called) / E$n, 0); E$x1 <- E$x0 + E$pct
E$facet <- ifelse(E$class == "POSITIVE", "Primary-labeled EC", "Other terminal EEC")
E <- E[order(E$class, E$family, E$seg), c("panel_size", "family", "label_source", "class", "rule", "seg", "called",
          "rule_total", "n", "pct", "gate_dropped", "marker_rule_fitted_here", "fam", "study_y", "y", "x0", "x1", "facet")]
write.table(E, a[["concordance"]], sep = "\t", row.names = FALSE, quote = TRUE)
