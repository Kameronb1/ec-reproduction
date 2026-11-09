# Figure 2E per-study scatters: every EC and other terminal EEC cell (primary study labels) of the four
# test studies and the McCauley 2023 HIO grafts, with its four-marker classifier score and result and
# its label-trained 20-gene classifier score, result and threshold. The corner counts of the panel
# count the cells each classifier classifies as EC.
#
# What it computes
#   A join of per-cell results of the other Figure 2 scripts; nothing is recalculated.
#   Test studies: four-marker score and result from heldout_marker_calls.R; 20-gene score, result and
#     threshold (panel D's 20 genes, weights and threshold fitted on the primary study labels of the
#     other studies) from primary_label_classifier.R.
#   HIO grafts: Endocrine_1 = EC, Endocrine_2 = other terminal EEC; four-marker score and call from
#     four_marker_fit.R; 20-gene score, result and threshold from hio_grafts.R.
#   Row order sets which dots are drawn on top: Elmentaite 2021, Hickey 2023, Singh 2024, Egozi 2021,
#   each in deposited cell order, then the HIO cells in deposited order.
#
# Inputs : cells (per-cell table from heldout_marker_calls.R), scores (per-cell table from
#          primary_label_classifier.R), hio (four-marker classifier results from four_marker_fit.R), hio_scores (per-cell
#          table from hio_grafts.R), hio_cells (cells of the GSE214852 endocrine export)
# Outputs: out (E_cells.tsv)
args <- commandArgs(trailingOnly = TRUE)
a <- setNames(sub("^[^=]*=", "", args), sub("=.*$", "", args))
rd <- function(p) read.delim(p, stringsAsFactors = FALSE, encoding = "UTF-8")
ROW_ORDER <- c("gut_cell_atlas", "hickey", "gse238274", "egozi")
CLASSES <- c("EXPLICIT_EC", "TERMINAL_NON_EC_EEC")

# ---- test studies --------------------------------------------------------------------------
cl <- rd(a[["cells"]])
stopifnot(setequal(cl$family, ROW_ORDER))
cl <- cl[order(match(cl$family, ROW_ORDER), seq_len(nrow(cl))), ]
cl <- cl[cl$author_label_class %in% CLASSES, ]
sc <- rd(a[["scores"]])
k <- match(cl$project_cell_id, sc$project_cell_id); stopifnot(!anyNA(k), !anyDuplicated(sc$project_cell_id))
tissue <- data.frame(fam = cl$family, class = cl$author_label_class, score4 = cl$ec_core_score, call4 = as.logical(cl$called),
                     score20 = sc$score[k], call20 = as.logical(sc$call[k]), thr20 = sc$threshold[k], stringsAsFactors = FALSE)

# ---- McCauley 2023 HIO grafts ------------------------------------------------------------------
hc <- rd(a[["hio_cells"]]); stopifnot(all(hc$author_label %in% c("Endocrine_1", "Endocrine_2")))
hs <- rd(a[["hio_scores"]]); stopifnot(identical(hs$cell_key, hc$figure2_cell_id))
hm <- rd(a[["hio"]]); k <- match(hc$figure2_cell_id, hm$cell_key); stopifnot(!anyNA(k))
hio <- data.frame(fam = "hio", class = ifelse(hc$author_label == "Endocrine_1", "EXPLICIT_EC", "TERMINAL_NON_EC_EEC"),
                  score4 = hm$ec_core_score[k], call4 = hm$presentation_class[k] == "HIGH_CONFIDENCE_EC",
                  score20 = hs$score, call20 = as.logical(hs$call), thr20 = hs$threshold, stringsAsFactors = FALSE)

out <- rbind(tissue, hio)
stopifnot(!anyNA(out))
write.table(out, a[["out"]], sep = "\t", row.names = FALSE, quote = FALSE)
