# Figure S4B: the 84 G3 genes (rows, three blocks of 28) in each study (columns), log2(1 + CPM) of counts pooled per study,
# in decreasing size of the study-level t statistic.
draw_panel <- function() {
  reset_metric_cache()
  letter_at("B", MARGIN + .02, B_HDR - .04); group_head("G3: separate SC-SIC program (84)", GROUP_COL[["G3"]], B_HDR - .01)
  for (k in 1:3) { bars(k, B_BAR); labels_down(k, B_ROWS_END + .04); cells(k, B_ROW1, G3O[(k - 1) * 28 + 1:28]) }
}
