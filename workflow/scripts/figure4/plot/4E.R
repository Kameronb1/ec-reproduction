# Figure 4E: the median expression difference from native intestinal EC (log2; 0, dashed = native intestinal EC) of SC-SIC, SC-beta
# and SC-alpha on Figure 3's three gene groups (one facet each), as an arrow from each row's start (open circle) to its
# end; behind each arrow a pale bar spans the middle half (25th to 75th percentile) of the group's genes at the end. Rows
# as in A, with grey bands on alternate rows. The plot shares its row heights with F (ef_layout() in helpers.R).

EF <- ef_layout()

draw_panel <- function() shifted_up(ROW_UP[["ef"]], {
  band_rects(EF$bxE, label_x = .06)
  draw_at(EF$gE, E_BOX[["x"]], EF$y, E_BOX[["w"]], EF$h)
  section_labels(EF$gE, EF$bxE[[1]], EF$strip_base)
})
