# Figure 2E: the label-trained 20-gene classifier as drawn, i.e. trained on primary study labels only.
#
# What it computes, for each study with primary study labels (Egozi 2021, Elmentaite 2021,
# Hickey 2023, Singh 2024)
#   D_GENES_PRIMARY_ONLY     panel D's 20 genes; weights and threshold refitted on the OTHER
#                            primary-label studies only (all three for Singh 2024). Drawn in E.
#                            D's genes were ranked with the labels of Egozi 2021, Elmentaite 2021
#                            and Hickey 2023, so for those three it is not an independent test.
#   FOLD_GENES_PRIMARY_ONLY  the 20-gene list the deposited leave-one-study-out computation ranked
#                            without this study, refitted the same way (an independent test for
#                            every study). Not drawn.
#   CONTROL                  the deposited computation's own classifier (genes and training studies),
#                            refitted here; it must reproduce the deposited results exactly. Not drawn.
# Results use the within-study neighborhood graphs of the deposited computation (Singh 2024's rebuilt
# from its cells, excluding the four markers and the top 30 ranked genes). Counts are cells classified
# as EC and all cells among the EC-labeled (POSITIVE) and other-terminal-EEC-labeled (NEGATIVE) cells.
# The drawn variant's per-cell score, threshold and result are also written (scores), for e_cells.R
# (panel E).
#
# Inputs : analysis, secondary (RDS), core (the deposited core result RDS, for its graphs), ranking
#          (gene_ranking.tsv), and the deposited nested_panels, nested_predictions, sv_predictions, sv_labels
# Outputs: calls (primary_only_training_calls.tsv), scores (per-cell D_GENES_PRIMARY_ONLY results)
args <- commandArgs(trailingOnly = TRUE)
a <- setNames(sub("^[^=]*=", "", args), sub("=.*$", "", args))
source(file.path(dirname(sub("^--file=", "", grep("^--file=", commandArgs(), value = TRUE))), "classifier_functions.R"))
gz <- function(p) read.delim(gzfile(p), stringsAsFactors = FALSE, check.names = FALSE)
PRIM <- c("SF_EGOZI2021_DEVELOPMENTAL_EEC", "SF_ELMENTAITE2021_GUT_ATLAS", "SF_HICKEY2023_ADULT_EEC")
SINGH <- "SF_INTESTINAL_GSE238274"

x <- readRDS(a[["analysis"]]); y <- readRDS(a[["secondary"]]); core <- readRDS(a[["core"]])
m <- x$metadata; labels <- training_labels(m)
ranking <- read.delim(a[["ranking"]], stringsAsFactors = FALSE)
TOP30 <- head(ranking$gene_symbol, 30L); D_GENES <- head(TOP30, 20L)
nested_panels <- gz(a[["nested_panels"]]); nested_pred <- gz(a[["nested_predictions"]])
sv_pred <- gz(a[["sv_predictions"]]); singh_lab <- gz(a[["sv_labels"]])

fit_and_call <- function(genes, train, expr, graph) {
  model <- tryCatch(fit_panel_model(x$normalized, m, labels, train, genes), error = function(e) e)
  if (inherits(model, "error")) return(list(error = conditionMessage(model)))
  sc <- score_panel_model(expr, model)
  list(call = gated_call(sc$score, sc$detected_genes, model$detection_requirement, model$threshold, graph)$call,
       threshold = model$threshold, score = sc$score)
}
cell_rows <- function(study, cells, res)
  data.frame(study = study, project_cell_id = cells, score = res$score, threshold = res$threshold, call = res$call,
             stringsAsFactors = FALSE)
count_rows <- function(study, variant, genes, res, cls, trained_on) {
  if (!is.null(res$error))
    return(data.frame(study = study, variant = variant, class = c("POSITIVE", "NEGATIVE"), called = NA, n = NA, threshold = NA,
                      genes_shared_with_D = sum(genes %in% D_GENES), trained_on = trained_on, note = res$error))
  do.call(rbind, lapply(c("POSITIVE", "NEGATIVE"), function(k)
    data.frame(study = study, variant = variant, class = k, called = sum(res$call[cls %in% k]), n = sum(cls %in% k),
               threshold = res$threshold, genes_shared_with_D = sum(genes %in% D_GENES), trained_on = trained_on, note = "")))
}

rows <- list(); scores <- list()
for (h in PRIM) {
  test <- m$family_id == h
  graph <- core$nested$graphs[[h]]; stopifnot(identical(graph$cells, m$project_cell_id[test]))
  fold <- nested_panels[nested_panels$outer_family == h & nested_panels$panel_id == "NESTED_SIZE_20", ]
  fold_genes <- fold$gene_symbol[order(fold$selected_rank)]; stopifnot(length(fold_genes) == 20L)
  cls <- labels$label[test]; expr <- x$normalized[, test, drop = FALSE]
  ctl <- fit_and_call(fold_genes, m$family_id != h, expr, graph)
  run <- nested_pred[nested_pred$outer_family == h & nested_pred$panel_id == "NESTED_SIZE_20", ]
  run <- run[match(m$project_cell_id[test], run$project_cell_id), ]
  stopifnot(identical(ctl$call, as.logical(run$stage_one_call) & as.logical(run$gate_pass)))   # reproduces the deposited results
  prim_train <- m$family_id %in% setdiff(PRIM, h)
  other <- paste(sort(setdiff(PRIM, h)), collapse = " + ")
  drawn <- fit_and_call(D_GENES, prim_train, expr, graph); stopifnot(is.null(drawn$error))
  scores[[h]] <- cell_rows(h, m$project_cell_id[test], drawn)
  rows[[h]] <- rbind(count_rows(h, "CONTROL", fold_genes, ctl, cls, "six other training studies (shipped)"),
                     count_rows(h, "D_GENES_PRIMARY_ONLY", D_GENES, drawn, cls, other),
                     count_rows(h, "FOLD_GENES_PRIMARY_ONLY", fold_genes, fit_and_call(fold_genes, prim_train, expr, graph), cls, other))
}

# Singh 2024 (not in training): graph from its own cells.
sgraph <- family_graph(y$normalized, unique(c(SELECTION_EXCLUSIONS, TOP30)))
scls <- singh_lab$comparison_label[match(y$metadata$project_cell_id, singh_lab$project_cell_id)]
all_cells <- rep(TRUE, nrow(m))
sctl <- fit_and_call(D_GENES, all_cells, y$normalized, sgraph)
sp <- sv_pred[sv_pred$panel_size == 20, ]; sp <- sp[match(y$metadata$project_cell_id, sp$project_cell_id), ]
stopifnot(identical(sctl$call, as.logical(sp$two_stage_call)))                                   # reproduces the deposited results
sdrawn <- fit_and_call(D_GENES, m$family_id %in% PRIM, y$normalized, sgraph); stopifnot(is.null(sdrawn$error))
scores[[SINGH]] <- cell_rows(SINGH, y$metadata$project_cell_id, sdrawn)
rows[[SINGH]] <- rbind(count_rows(SINGH, "CONTROL", D_GENES, sctl, scls, "all seven training studies (shipped)"),
                       count_rows(SINGH, "D_GENES_PRIMARY_ONLY", D_GENES, sdrawn, scls, paste(PRIM, collapse = " + ")))

res <- do.call(rbind, rows); rownames(res) <- NULL
res$pct <- round(100 * res$called / res$n, 1)
write.table(res, a[["calls"]], sep = "\t", quote = FALSE, row.names = FALSE)
write.table(do.call(rbind, scores), a[["scores"]], sep = "\t", quote = FALSE, row.names = FALSE)
print(res[c("study", "variant", "class", "called", "n", "pct")], row.names = FALSE)
