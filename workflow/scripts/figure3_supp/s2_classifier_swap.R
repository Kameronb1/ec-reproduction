# Figure S6: "Does the 20-gene classifier find SC-SIC?" (gene list x training set swap grid).
#
# What it computes. Figure 2's label-trained 20-gene classifier is scored on each of the six SC-islet studies, with the
# gene list and the training set swapped:
#   GUTG_GUTT  Figure 2 panel D's 20 intestinal genes, trained on the three intestinal studies with primary study labels
#              (Egozi 2021, Elmentaite 2021, Hickey 2023; EC vs other enteroendocrine) = Figure 2E's classifier
#   GUTG_SCT   the same 20 intestinal genes, trained on the other five SC-islet studies (SC-SIC vs other SC endocrine)
#   SCG_GUTT   Figure 3's 20 SC-SIC genes, trained on the three intestinal studies (CFAP92 read from KIAA1257 in
#              Hickey 2023 and Egozi 2021)
#   SCG_SCT    the 20 SC-SIC genes, trained on the other five SC-islet studies
# (GUT/SC in the codes = intestinal/SC-islet gene list (G) or training set (T).)
# A fit sets gene weights (median over training studies of the per-sample EC-minus-other mean, negatives set to 0,
# summing to 1), a per-gene standardizer and a two-component mixture threshold; a cell is classified as EC when it meets
# the threshold, has at least 25% of the panel genes detected, and passes the neighborhood check on its own study's
# graph (k 25, BH q <= 0.05). The four-marker classifier (Figure 2A's fixed cutoff), with and without its
# neighborhood check, is reported beside them, and the TPH1-only classifier (TPH1 alone, fitted on the same cells by
# the Figure 2 supplement; at least one TPH1 read) with the four-marker classifier's neighborhood check.
# Output: cells classified as EC / n per test study and class, and for Figure S6B every SC endocrine cell's four-marker
# score and classification and its score and classification from Figure 2E's classifier (GUTG_GUTT).
#
# Inputs (name=path):
#   lib          s2_lib.R (the classifier functions)
#   gut_cells    intestinal training cells: project_cell_id, family_id, biological_unit_id, training label
#   gut_expr     log-normalized expression of the panel genes in those cells (long format, non-zero entries)
#   sc_meta      SC endocrine cell table (metadata.tsv: study, sample, primary-study-label class, total UMIs)
#   sc_counts    raw UMI counts of the panel and marker genes in those cells (long format, non-zero entries)
#   sc_graphs    neighborhood graphs of each SC study (classifier check and four-marker check)
#   markers      the four-marker classifier's fitted constants (Figure 2A)
#   tph1_fit     the TPH1-only classifier's fitted constants (Figure 2 supplement, tph1_only/four_marker_fit.json)
#   gut_genes    Figure 2 panel D gene table (D_effect_heatmap.tsv; 20 genes by panel_rank)
#   sc_genes     Figure 3's 20 SC-SIC genes (top20.tsv)
# Outputs:
#   out          calls_by_study.tsv
#   cells        sc_cells.tsv.gz: one row per SC endocrine cell (studies in turn, cells in metadata order): study, class
#                (SC_SIC or OTHER), four-marker score, cutoff and call, 20-gene score, threshold and call
a <- commandArgs(TRUE); A <- setNames(sub("^[^=]*=", "", a), sub("=.*", "", a))
suppressPackageStartupMessages({ library(Matrix); library(jsonlite) })
source(A[["lib"]])

# Classifier settings, as in Figure 2's 20-gene classifier: mixture fit and neighborhood check.
gmm <- list(deterministic_initial_quantile_pairs = list(c(0.1, 0.9), c(0.25, 0.75), c(0.4, 0.8)), max_iterations = 1000L,
            tolerance = 1e-08, variance_floor = 0.0025, minimum_component_weight = 0.02, minimum_component_sd = 0.05,
            minimum_mean_separation_pooled_sd = 0.5)
gate <- list(k = 25L, q_threshold = 0.05)

# ---- gene lists
dt <- read.delim(A[["gut_genes"]], stringsAsFactors = FALSE)
D_GENES <- unique(dt$gene_symbol[order(dt$panel_rank)]); stopifnot(length(D_GENES) == 20L)
top20 <- read.delim(A[["sc_genes"]], stringsAsFactors = FALSE)
SC_GENES <- top20$gene[order(top20$rank)]; stopifnot(length(SC_GENES) == 20L)

# ---- intestinal training data (three intestinal studies with primary study labels; CFAP92 included)
m <- read.delim(A[["gut_cells"]], stringsAsFactors = FALSE, colClasses = "character")
labels <- data.frame(project_cell_id = m$project_cell_id, label = m$label, stringsAsFactors = FALSE)
gv <- read.delim(A[["gut_expr"]], stringsAsFactors = FALSE, colClasses = c("character", "character", "numeric"))
ggenes <- sort(unique(gv$gene), method = "radix")
gut <- sparseMatrix(i = match(gv$gene, ggenes), j = match(gv$project_cell_id, m$project_cell_id), x = gv$value,
                    dims = c(length(ggenes), nrow(m)), dimnames = list(ggenes, m$project_cell_id))
PRIM <- c("SF_EGOZI2021_DEVELOPMENTAL_EEC", "SF_ELMENTAITE2021_GUT_ATLAS", "SF_HICKEY2023_ADULT_EEC")
stopifnot(nrow(m) == 5214L, all(m$family_id %in% PRIM), all(c(D_GENES, SC_GENES) %in% ggenes), !anyNA(gv$value))

# ---- SC cells: ln(1 + count / cell's total UMIs x 1e6), as the SC input package is normalized
meta <- read.delim(A[["sc_meta"]], stringsAsFactors = FALSE, colClasses = "character")
meta$total_umis <- as.numeric(meta$total_umis)
# The neighborhood graphs list each study's cells in their own order; reorder the rows within each (contiguous) study
# block to match, which leaves every block where it is
local({ sg0 <- readRDS(A[["sc_graphs"]])$classifier_gate; for (f in names(sg0)) { idx <- which(meta$family_id == f)
  stopifnot(all(diff(idx) == 1L)); o <- match(sg0[[f]]$cells, meta$project_cell_id[idx]); stopifnot(!anyNA(o), length(o) == length(idx))
  meta[idx, ] <<- meta[idx[o], ] } })
sv <- read.delim(A[["sc_counts"]], stringsAsFactors = FALSE, colClasses = c("character", "character", "numeric"))
sgenes <- sort(unique(sv$gene), method = "radix")
j <- match(sv$project_cell_id, meta$project_cell_id); stopifnot(!anyNA(j))
expr <- sparseMatrix(i = match(sv$gene, sgenes), j = j, x = log1p(sv$count / meta$total_umis[j] * 1e6),
                     dims = c(length(sgenes), nrow(meta)), dimnames = list(sgenes, meta$project_cell_id))
stopifnot(nrow(meta) == 55120L, all(is.finite(expr@x)))
STUDIES <- sort(unique(meta$family_id)); stopifnot(length(STUDIES) == 6L)
POS <- "SC_SIC"; NEG <- c("SC_BETA", "SC_ALPHA", "SC_DELTA", "SC_EPSILON", "SC_GAMMA")
stopifnot(all(meta$target %in% c(POS, NEG)))
sc_labels <- data.frame(project_cell_id = meta$project_cell_id, label = ifelse(meta$target == POS, "POSITIVE", "NEGATIVE"), stringsAsFactors = FALSE)

sg <- readRDS(A[["sc_graphs"]])
graphs <- sg$classifier_gate; mgraphs <- sg$marker_rule_veto
stopifnot(identical(names(graphs), STUDIES), identical(names(mgraphs), STUDIES))
for (f in STUDIES) stopifnot(identical(graphs[[f]]$cells, meta$project_cell_id[meta$family_id == f]),
                             identical(mgraphs[[f]]$cols, which(meta$family_id == f)))

# ---- fit and classify (as fit_and_call in workflow/scripts/figure2/primary_label_classifier.R, with the training data
# as an argument)
fit_model <- function(expr, genes_all, meta, lab, train, genes) {
  fit <- tryCatch(cp_fit_panel_model(expr, genes_all, meta, lab, train, genes, gmm), error = function(e) e)
  if (inherits(fit, "error")) list(error = conditionMessage(fit)) else list(model = fit$model)
}
call_with <- function(fm, expr, graph) {
  if (!is.null(fm$error)) return(fm)
  sc <- cp_score_panel_model(expr, fm$model)
  cl <- cp_two_stage_call(sc$score, sc$detected_genes, fm$model$detection_requirement, fm$model$gmm$threshold, 0, graph, gate)
  list(model = fm$model, score = sc$score, detected = sc$detected_genes, stage_one = cl$stage_one, call = cl$call,
       threshold = fm$model$gmm$threshold)
}
GUT_FIT_D <- fit_model(gut, rownames(gut), m, labels, m$family_id %in% PRIM, D_GENES)
GUT_FIT_SC20 <- fit_model(gut, rownames(gut), m, labels, m$family_id %in% PRIM, SC_GENES)
stopifnot(is.null(GUT_FIT_D$error), is.null(GUT_FIT_SC20$error))
fold_fit <- function(combo, f) switch(combo,
  GUTG_GUTT = GUT_FIT_D,
  SCG_GUTT = GUT_FIT_SC20,
  GUTG_SCT = fit_model(expr, rownames(expr), meta, sc_labels, meta$family_id != f, D_GENES),
  SCG_SCT = fit_model(expr, rownames(expr), meta, sc_labels, meta$family_id != f, SC_GENES))
run_combo <- function(combo) lapply(setNames(STUDIES, STUDIES), function(f) {
  test <- meta$family_id == f
  call_with(fold_fit(combo, f), expr[, test, drop = FALSE], graphs[[f]])
})
COMBOS <- data.frame(
  combination = c("GUTG_GUTT", "GUTG_SCT", "SCG_GUTT", "SCG_SCT"),
  label = c("Gut genes + gut training", "Gut genes + SC training", "SC genes + gut training", "SC genes + SC training (= A)"),
  stringsAsFactors = FALSE)
RES <- lapply(setNames(COMBOS$combination, COMBOS$combination), run_combo)

# ---- the four-marker classifier: score of the four marker z-scores >= cut-off and >= 2 of 4 markers detected,
# then its own neighborhood check within each SC study (MARKER_RULE_GATED); without the check: MARKER_RULE_UNGATED
mc <- fromJSON(A[["markers"]])
CORE4 <- mc$ec_core_genes; stopifnot(identical(sort(CORE4), c("FEV", "LMX1A", "SLC18A1", "TPH1")))
st <- mc$standardisation; rownames(st) <- st$gene
Z <- vapply(CORE4, function(g) (as.numeric(expr[g, ]) - st[g, "mean"]) / st[g, "sd"], numeric(ncol(expr)))
mk_score <- rowMeans(Z)
mk_detected <- as.numeric(Matrix::colSums(expr[CORE4, , drop = FALSE] > 0))
mk_candidate <- mk_score >= mc$ec_core_cutoff & mk_detected >= mc$min_markers_detected
mk_high <- ec_graph_veto(mk_candidate, mgraphs, mc$neighbourhood$q_threshold)
# ---- the TPH1-only classifier: TPH1 z-score >= its cutoff and >= 1 TPH1 read, then the same check (TPH1_GATED)
tc <- fromJSON(A[["tph1_fit"]]); stopifnot(identical(tc$ec_core_genes, "TPH1"))
ts <- tc$standardisation; rownames(ts) <- ts$gene
t1_candidate <- (as.numeric(expr["TPH1", ]) - ts["TPH1", "mean"]) / ts["TPH1", "sd"] >= tc$ec_core_cutoff &
  as.numeric(expr["TPH1", ] > 0) >= tc$min_markers_detected
t1_high <- ec_graph_veto(t1_candidate, mgraphs, mc$neighbourhood$q_threshold)

# ---- table: one row per combination x test study x class
NAME <- c(GSE151117 = "Augsornworawat 2020 HUES8", GSE167880 = "Balboa 2022", GSE185036 = "Du 2022", GSE199636 = "Augsornworawat 2023",
          GSE202497 = "Zhu 2023", GSE226346 = "Siehler 2024")   # study names by first author and year
CLS <- c(POSITIVE = "SC_SIC", NEGATIVE = "OTHER_MATURE_SC_ENDOCRINE")
calls <- list()
for (i in seq_len(nrow(COMBOS))) {
  cmb <- COMBOS$combination[i]
  for (f in STUDIES) {
    r <- RES[[cmb]][[f]]; test <- meta$family_id == f; lab <- sc_labels$label[test]; resp <- as.integer(lab == "POSITIVE")
    if (!is.null(r$error)) {
      calls[[paste(cmb, f)]] <- data.frame(combination = cmb, combination_label = COMBOS$label[i], held_out_study = f, study_name = NAME[[f]],
        class = CLS, called = NA, n = c(sum(resp == 1L), sum(resp == 0L)), pct = NA, stage_one_called = NA, threshold = NA,
        detection_requirement = NA, n_genes = NA, positive_weight_genes = NA, trained_on = NA, status = "FIT_FAILED", note = r$error)
      next
    }
    md <- r$model
    calls[[paste(cmb, f)]] <- data.frame(combination = cmb, combination_label = COMBOS$label[i], held_out_study = f, study_name = NAME[[f]],
      class = CLS, called = c(sum(r$call[resp == 1L]), sum(r$call[resp == 0L])), n = c(sum(resp == 1L), sum(resp == 0L)),
      pct = NA, stage_one_called = c(sum(r$stage_one[resp == 1L]), sum(r$stage_one[resp == 0L])), threshold = r$threshold,
      detection_requirement = md$detection_requirement, n_genes = length(md$panel_genes), positive_weight_genes = sum(md$gene_weights > 0),
      trained_on = paste(md$training_families, collapse = ";"), status = "PASS", note = "")
  }
}
for (f in STUDIES) {
  test <- meta$family_id == f; lab <- sc_labels$label[test]; resp <- as.integer(lab == "POSITIVE")
  for (rule in c("MARKER_RULE_GATED", "MARKER_RULE_UNGATED")) {
    cl <- if (rule == "MARKER_RULE_GATED") mk_high[test] else mk_candidate[test]
    calls[[paste(rule, f)]] <- data.frame(combination = rule,
      combination_label = if (rule == "MARKER_RULE_GATED") "Four-marker classifier with its local-consistency check" else "Four-marker classifier without its local-consistency check (as Figure 4)",
      held_out_study = f, study_name = NAME[[f]], class = CLS, called = c(sum(cl[resp == 1L]), sum(cl[resp == 0L])), n = c(sum(resp == 1L), sum(resp == 0L)),
      pct = NA, stage_one_called = c(sum(mk_candidate[test][resp == 1L]), sum(mk_candidate[test][resp == 0L])), threshold = mc$ec_core_cutoff,
      detection_requirement = mc$min_markers_detected, n_genes = 4L, positive_weight_genes = 4L, trained_on = "fixed: fitted on the 1,199 primary intestinal cells of the fitting studies (Figure 2A)",
      status = "PASS", note = "")
  }
  calls[[paste("TPH1_GATED", f)]] <- data.frame(combination = "TPH1_GATED",
    combination_label = "TPH1-only classifier with the four-marker classifier's local-consistency check",
    held_out_study = f, study_name = NAME[[f]], class = CLS, called = c(sum(t1_high[test][resp == 1L]), sum(t1_high[test][resp == 0L])),
    n = c(sum(resp == 1L), sum(resp == 0L)), pct = NA,
    stage_one_called = c(sum(t1_candidate[test][resp == 1L]), sum(t1_candidate[test][resp == 0L])), threshold = tc$ec_core_cutoff,
    detection_requirement = tc$min_markers_detected, n_genes = 1L, positive_weight_genes = 1L,
    trained_on = "fixed: fitted on the 1,199 primary intestinal cells of the fitting studies (Figure 2A), TPH1 alone",
    status = "PASS", note = "")
}
calls <- do.call(rbind, calls); rownames(calls) <- NULL; calls$pct <- round(100 * calls$called / calls$n, 2)
stopifnot(nrow(calls) == 84L)
con <- file(A[["out"]], "wb")
write.table(calls, con, sep = "\t", quote = FALSE, row.names = FALSE, na = "NA", eol = "\n")
close(con)

# ---- per-cell table (Figure S6B): Figure 2E's classifier (intestinal genes and training) and the four-marker classifier
cells <- do.call(rbind, lapply(STUDIES, function(f) {
  r <- RES$GUTG_GUTT[[f]]; test <- meta$family_id == f
  data.frame(study = f, class = ifelse(sc_labels$label[test] == "POSITIVE", "SC_SIC", "OTHER"), score4 = mk_score[test],
             cut4 = mc$ec_core_cutoff, call4 = mk_high[test], score20 = r$score, thr20 = r$threshold, call20 = r$call)
}))
stopifnot(nrow(cells) == nrow(meta))
con <- gzfile(A[["cells"]], "wb")
write.table(cells, con, sep = "\t", quote = FALSE, row.names = FALSE, eol = "\n")
close(con)
