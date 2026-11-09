# Figures S3C, S4C-D, S5B-C: the 20-gene classifier scored on the four test studies
# with primary study labels and on McCauley 2023 HIO grafts.
#
# The classifier: the top 20 genes of the gene ranking on all training studies (study-weighted
# limma, as in Figure 2's leave-one-study-out run; sizes 10 and 30 are also scored), with gene
# weights and a Gaussian-mixture score threshold fitted only on cells with primary study labels from
# the OTHER primary-label studies (Egozi 2021, Elmentaite 2021, Hickey 2023; all three for Singh 2024
# and McCauley 2023), at least 25% of the genes detected, and the within-study neighborhood check.
# McCauley is trained on the same cells as Singh, so every McCauley number uses Singh's fitted models.
#
# Tables
#   predictions / hio_predictions   per cell, sizes 10, 20, 30                       (S5C)
#   weights / hio_weights           genes, weights, training effects                  (S5B)
#   matched_random / hio_matched    200 expression-, detection- and variance-matched random 20-gene sets (S4C);
#                                   written only when these two arguments are given (about 55 min on one core)
#   generic / hio_generic           CHGA, CHGB, SCG2, SCG3, PCSK1, PCSK2 fitted the same way (S4D)
#   depth / hio_depth               scored study's counts thinned to 50% / 25%, 3 draws, graph rebuilt, model fixed (S3C)
#   hio_by_organoid, hio_ranking    McCauley per organoid and AUC / sensitivity / specificity
#
# The classifier's functions are in classifier_lib/. Figure 2's leave-one-study-out run used R 4.2.0,
# Matrix 1.6.1.1 and single-threaded BLAS; this script attaches Matrix and runs on the installed versions.
#
# Inputs : core (core result of Figure 2's leave-one-study-out run, RDS: the test studies' within-study
#          graphs), analysis (training corpus, RDS), secondary (Singh 2024, RDS), singh_labels,
#          hio_counts, hio_genes, hio_cells (McCauley export), hio_classes (from marker_fit.R)
args <- commandArgs(trailingOnly = TRUE)
a <- setNames(sub("^[^=]*=", "", args), sub("=.*$", "", args))
HERE <- dirname(sub("^--file=", "", grep("^--file=", commandArgs(), value = TRUE)))
pkg <- file.path(HERE, "classifier_lib")
suppressPackageStartupMessages({library(Matrix); library(jsonlite)})
for (f in c("figure2_postvalidation_core.R", "figure2_postvalidation_analysis.R"))
  source(file.path(pkg, "vendor/base/R", f))
source(file.path(pkg, "R/candidate_panel_core.R")); source(file.path(pkg, "R/population_policy.R"))
gz <- function(p) read.delim(gzfile(p), stringsAsFactors = FALSE, check.names = FALSE)
wt <- function(x, path) {
  if (grepl("[.]gz$", path)) { con <- gzfile(path, "w"); write.table(x, con, sep = "\t", quote = TRUE, row.names = FALSE); close(con) }
  else write.table(x, path, sep = "\t", quote = TRUE, row.names = FALSE)
}
log_step <- function(...) cat(format(Sys.time(), "%H:%M:%S"), ..., "\n")

core <- readRDS(a[["core"]])
x <- readRDS(a[["analysis"]]); classifier_assert_training_data(x)
y <- readRDS(a[["secondary"]]); classifier_assert_secondary_data(y)
m <- x$metadata; labels <- classifier_preflight_labels(m)
contract <- cp_read_json(file.path(pkg, "config/candidate_panel_contract.json"))
parameters <- read_numerical_contract(file.path(pkg, "vendor/base/config/numerical_contract.json"))$parameters
g <- contract$gmm
gmm <- list(deterministic_initial_quantile_pairs = g$initial_quantile_pairs, max_iterations = as.integer(g$max_iterations),
            tolerance = as.numeric(g$tolerance), variance_floor = as.numeric(g$variance_floor),
            minimum_component_weight = as.numeric(g$minimum_component_weight), minimum_component_sd = as.numeric(g$minimum_component_sd),
            minimum_mean_separation_pooled_sd = as.numeric(g$minimum_mean_separation_pooled_sd))
gate <- cp_gate_spec(contract)
PRIM <- c("SF_EGOZI2021_DEVELOPMENTAL_EEC", "SF_ELMENTAITE2021_GUT_ATLAS", "SF_HICKEY2023_ADULT_EEC")
SINGH <- "SF_INTESTINAL_GSE238274"; HIO <- "SF_GSE214852_HIO_GRAFT"
STUDIES <- c(PRIM, SINGH)
EXCL <- sort(unique(classifier_selection_exclusions))

# ---- the ranking, the graphs and the scored populations ----------------------------------------
all_cells <- rep(TRUE, nrow(m))
ranked <- cp_rank_training_features(x$normalized, x$genes, m, labels, all_cells, core$excluded_genes, parameters, 30L, NULL,
                                    as.numeric(contract$gene_ranking$redundancy_ceiling))
TOP <- ranked$genes; stopifnot(length(TOP) >= 30L)
D_GENES <- head(TOP, 20L)
GRAPH_EXCL <- unique(c(core$excluded_genes, TOP))
singh_lab <- gz(a[["singh_labels"]])
sgraph <- cp_family_graph(y$normalized, rep(TRUE, nrow(y$metadata)), GRAPH_EXCL, gate)
pop <- lapply(setNames(STUDIES, STUDIES), function(s) {
  if (s == SINGH) {
    list(expr = y$normalized, counts = y$counts, graph = sgraph, cells = y$metadata$project_cell_id,
         unit = y$metadata$biological_unit_id, rule = y$metadata$presentation_class,
         label = singh_lab$comparison_label[match(y$metadata$project_cell_id, singh_lab$project_cell_id)],
         graph_exclusions = GRAPH_EXCL, train = m$family_id %in% PRIM)
  } else {
    test <- m$family_id == s
    gr <- core$nested$graphs[[s]]; stopifnot(identical(gr$cells, m$project_cell_id[test]))
    list(expr = x$normalized[, test, drop = FALSE], counts = x$counts[, test, drop = FALSE], graph = gr,
         cells = m$project_cell_id[test], unit = m$biological_unit_id[test], rule = m$presentation_class[test],
         label = labels$label[test], graph_exclusions = gr$excluded_genes, train = m$family_id %in% setdiff(PRIM, s))
  }
})

# McCauley 2023: log1p counts per million of the whole library (library total from the HIO table of marker_fit.R)
cnt <- as(Matrix::readMM(a[["hio_counts"]]), "CsparseMatrix")
genes_all <- readLines(a[["hio_genes"]]); hcells <- read.delim(a[["hio_cells"]], stringsAsFactors = FALSE)
if (nrow(cnt) != length(genes_all)) cnt <- Matrix::t(cnt)
stopifnot(ncol(cnt) == 598L, !anyDuplicated(genes_all))
dimnames(cnt) <- list(genes_all, hcells$figure2_cell_id)
hcells$sample <- sub("[.][0-9]+$", "", hcells$sample)               # organoids wt1 (wt1.1 + wt1.2) and wt2
hclass <- gz(a[["hio_classes"]]); hclass <- hclass[match(hcells$figure2_cell_id, hclass$cell_key), ]
stopifnot(!anyNA(hclass$cell_key), all(hclass$arm_id == "GSE214852_WT"))
lib <- hclass$gex_total_counts
hio_label <- ifelse(hcells$author_label == "Endocrine_1", "POSITIVE", "NEGATIVE")
stopifnot(all(hcells$author_label %in% c("Endocrine_1", "Endocrine_2")), sum(hio_label == "POSITIVE") == 322L)
keep_genes <- intersect(x$genes, genes_all); stopifnot(all(TOP %in% keep_genes))
hcounts <- cnt[keep_genes, , drop = FALSE]
z <- hcounts %*% Matrix::Diagonal(x = 1e6 / lib); z@x <- log1p(z@x); dimnames(z) <- list(keep_genes, hcells$figure2_cell_id)
pop[[HIO]] <- list(expr = z, counts = hcounts, graph = cp_family_graph(z, rep(TRUE, ncol(z)), GRAPH_EXCL, gate),
                   cells = hcells$figure2_cell_id, unit = hcells$sample, rule = hclass$presentation_class, label = hio_label,
                   graph_exclusions = GRAPH_EXCL, train = pop[[SINGH]]$train)

# ---- fitting and scoring ----------------------------------------------------------------------
cp_score_panel_model_free <- function(expr, model) {
  standardized <- cp_apply_standardizer(expr, model$standardizer)
  index <- match(model$panel_genes, rownames(expr)); stopifnot(!anyNA(index))
  list(score = as.numeric(Matrix::crossprod(model$gene_weights, standardized)),
       detected_genes = as.numeric(Matrix::colSums(expr[index, , drop = FALSE] > 0)))
}
fit_only <- function(genes, s, detection_fraction = 0.25, null_policy = FALSE) {
  fit <- tryCatch(cp_fit_panel_model(x$normalized, x$genes, m, labels, pop[[s]]$train, genes, gmm,
                                     detection_fraction = detection_fraction, null_panel_policy = null_policy),
                  error = function(e) e)
  if (inherits(fit, "error")) list(error = conditionMessage(fit)) else fit
}
score_fit <- function(fit, s) {
  if (!is.null(fit$error)) return(list(error = fit$error))
  P <- pop[[s]]
  sc <- cp_score_panel_model_free(P$expr, fit$model); thr <- fit$model$gmm$threshold
  call <- if (is.finite(thr)) cp_two_stage_call(sc$score, sc$detected_genes, fit$model$detection_requirement, thr, 0, P$graph, gate) else NULL
  list(model = fit$model, score = sc$score, detected = sc$detected_genes, threshold = thr,
       stage_one = if (is.null(call)) rep(NA, length(sc$score)) else call$stage_one,
       call = if (is.null(call)) rep(NA, length(sc$score)) else call$call)
}
metrics <- function(s, r) {
  P <- pop[[s]]; keep <- P$label %in% c("POSITIVE", "NEGATIVE"); resp <- as.integer(P$label[keep] == "POSITIVE")
  if (!is.null(r$error)) return(data.frame(status = "FIT_FAILED", note = r$error, roc_auc = NA, sensitivity = NA, specificity = NA,
                                           ec_called = NA, ec_n = sum(resp == 1L), other_called = NA, other_n = sum(resp == 0L), threshold = NA))
  called <- r$call[keep]
  data.frame(status = if (is.finite(r$threshold)) "PASS" else "NO_THRESHOLD", note = "", roc_auc = binary_auc(resp, r$score[keep]),
             sensitivity = if (anyNA(called)) NA else mean(called[resp == 1L]),
             specificity = if (anyNA(called)) NA else mean(!called[resp == 0L]),
             ec_called = if (anyNA(called)) NA else sum(called[resp == 1L]), ec_n = sum(resp == 1L),
             other_called = if (anyNA(called)) NA else sum(called[resp == 0L]), other_n = sum(resp == 0L), threshold = r$threshold)
}
pred_rows <- function(s, k, r) { P <- pop[[s]]
  data.frame(study = s, panel_size = k, project_cell_id = P$cells, biological_unit_id = P$unit, primary_label = P$label,
             rule_grade = P$rule, score = r$score, detected_genes = r$detected, detection_requirement = r$model$detection_requirement,
             threshold = r$threshold, stage_one_call = r$stage_one, call = r$call, stringsAsFactors = FALSE) }
weight_rows <- function(s, k, r) data.frame(study = s, panel_size = k, rank = seq_len(k), gene_symbol = r$model$panel_genes,
  weight = r$model$gene_weights, training_effect = as.numeric(r$model$training_family_effects),
  trained_on = paste(r$model$training_families, collapse = ";"), stringsAsFactors = FALSE)

# ---- 1. predictions and weights at sizes 10 / 20 / 30 ------------------------------------------
log_step("predictions")
pred <- list(); wts <- list(); base <- list(); hpred <- list(); hwts <- list()
for (s in STUDIES) for (k in c(10L, 20L, 30L)) {
  fit <- fit_only(head(TOP, k), s); r <- score_fit(fit, s); stopifnot(is.null(r$error), is.finite(r$threshold))
  pred[[paste(s, k)]] <- pred_rows(s, k, r); wts[[paste(s, k)]] <- weight_rows(s, k, r)
  if (k == 20L) base[[s]] <- r
  if (s == SINGH) {
    rh <- score_fit(fit, HIO)
    hpred[[as.character(k)]] <- pred_rows(HIO, k, rh); hwts[[as.character(k)]] <- weight_rows(HIO, k, rh)
    if (k == 20L) base[[HIO]] <- rh
  }
}
pred <- do.call(rbind, pred); hpred <- do.call(rbind, hpred)
wt(pred, a[["predictions"]]); wt(do.call(rbind, wts), a[["weights"]])
wt(hpred, a[["hio_predictions"]]); wt(do.call(rbind, hwts), a[["hio_weights"]])

# ---- 2. matched random 20-gene sets (McCauley: Singh's 200 sets and fits) ------------------------
eligible <- x$genes[training_feature_mask(x$normalized, x$genes, m, all_cells, as.numeric(parameters$training_detection_floor),
                                          as.integer(parameters$minimum_training_families), EXCL)]
eligible <- intersect(eligible, y$genes); stopifnot(all(D_GENES %in% eligible))
if ("matched_random" %in% names(a)) {   # about 55 min single-threaded; the workflow uses the deposited tables
log_step("matched random")
REPS <- 200L
mr <- list(); mr_hio <- NULL
for (s in STUDIES) {
  draws <- cp_draw_matched_random_panels(x$normalized, x$genes, m, pop[[s]]$train, D_GENES, eligible, REPS,
                                         seed = 20261004L + match(s, STUDIES) * 10000L, forbidden_genes = EXCL)
  rows <- lapply(seq_len(REPS), function(i) {
    fit <- fit_only(draws$panels$gene_symbol[draws$panels$replicate == i], s, null_policy = TRUE)
    out <- cbind(data.frame(study = s, replicate = i, panel = "MATCHED_RANDOM"), metrics(s, score_fit(fit, s)))
    if (s == SINGH) attr(out, "hio") <- cbind(data.frame(study = HIO, replicate = i, panel = "MATCHED_RANDOM"), metrics(HIO, score_fit(fit, HIO)))
    out
  })
  mr[[s]] <- rbind(cbind(data.frame(study = s, replicate = 0L, panel = "REAL"), metrics(s, base[[s]])), do.call(rbind, rows))
  if (s == SINGH) mr_hio <- rbind(cbind(data.frame(study = HIO, replicate = 0L, panel = "REAL"), metrics(HIO, base[[HIO]])),
                                  do.call(rbind, lapply(rows, attr, "hio")))
  log_step(" ", s)
}
wt(do.call(rbind, mr), a[["matched_random"]]); wt(mr_hio, a[["hio_matched"]])
}

# ---- 3. generic neuroendocrine genes ------------------------------------------------------------
GENERIC <- c("CHGA", "CHGB", "SCG2", "SCG3", "PCSK1", "PCSK2"); stopifnot(all(GENERIC %in% eligible))
gen_rows <- function(s, gfit) rbind(
  cbind(data.frame(study = s, panel = "20_GENE_PANEL", genes = paste(D_GENES, collapse = ";")), metrics(s, base[[s]])),
  cbind(data.frame(study = s, panel = "GENERIC_NEUROENDOCRINE", genes = paste(GENERIC, collapse = ";")), metrics(s, score_fit(gfit, s))))
gfits <- lapply(setNames(STUDIES, STUDIES), function(s) fit_only(GENERIC, s))
wt(do.call(rbind, lapply(STUDIES, function(s) gen_rows(s, gfits[[s]]))), a[["generic"]])
wt(gen_rows(HIO, gfits[[SINGH]]), a[["hio_generic"]])

# ---- 4. sequencing depth: thin the scored study only, rebuild its graph, model fixed --------------
log_step("depth")
depth_rows <- function(s, idx) {
  P <- pop[[s]]; model <- base[[s]]$model; out <- list()
  # stored expression = log1p(counts x per-cell factor); recover the factor and rescale it after thinning
  csum <- as.numeric(Matrix::colSums(P$counts))
  sf <- as.numeric(Matrix::colSums(expm1(P$expr))) / csum
  renorm <- function(cn, f) { zz <- cn %*% Matrix::Diagonal(x = f); zz@x <- log1p(zz@x); dimnames(zz) <- dimnames(cn); zz }
  stopifnot(max(abs(renorm(P$counts, sf) - P$expr)) < 1e-6)
  for (fr in c(.5, .25)) for (rep in 1:3) {
    tc <- cp_downsample_counts(P$counts, fr, seed = 20261004L + rep + 100L * idx + round(1000 * fr))
    thin <- renorm(tc, sf * csum / pmax(as.numeric(Matrix::colSums(tc)), 1))
    gr <- cp_family_graph(thin, rep(TRUE, ncol(thin)), P$graph_exclusions, gate)
    sc <- cp_score_panel_model_free(thin, model)
    cl <- cp_two_stage_call(sc$score, sc$detected_genes, model$detection_requirement, model$gmm$threshold, 0, gr, gate)
    out[[paste(fr, rep)]] <- cbind(data.frame(study = s, fraction = fr, replicate = rep),
      metrics(s, list(score = sc$score, call = cl$call, threshold = model$gmm$threshold)))
  }
  out[["1"]] <- cbind(data.frame(study = s, fraction = 1, replicate = 0L), metrics(s, base[[s]]))
  do.call(rbind, out)
}
wt(do.call(rbind, lapply(STUDIES, function(s) depth_rows(s, match(s, STUDIES)))), a[["depth"]])
wt(depth_rows(HIO, length(STUDIES) + 1L), a[["hio_depth"]])          # McCauley = 5th scored study

# ---- 5. McCauley per organoid, and ranking versus calling -----------------------------------------
p20 <- hpred[hpred$panel_size == 20L, ]; ec <- p20[p20$primary_label == "POSITIVE", ]
wt(do.call(rbind, lapply(c("wt1", "wt2"), function(u) { k <- ec$biological_unit_id == u
  data.frame(study = HIO, unit = u, n_rule = sum(k), called_rule = sum(ec$rule_grade[k] == "HIGH_CONFIDENCE_EC"),
             called_cls = sum(ec$call[k]), stringsAsFactors = FALSE) })), a[["hio_by_organoid"]])
resp <- p20$primary_label == "POSITIVE"
wt(data.frame(study = HIO, metric = c("Ranking: ROC AUC", "Calling: sensitivity", "Calling: specificity"),
              value = c(binary_auc(as.integer(resp), p20$score), mean(p20$call[resp]), mean(!p20$call[!resp]))), a[["hio_ranking"]])
log_step("classifier done")
