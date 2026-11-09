# Figure 4F: Figure 3H's fetal–adult score of SC-SIC and SC-beta, as an arrow from the mean over each
# row's start samples (open circle) to the mean over its end samples. Ticks on the axis: Figure 3H's intestinal rows on the
# same genes (fetal orange, adult teal; long: the intestinal arms that chose the genes, short: the intestinal test sources). Rows as in
# A, with grey bands on alternate rows. The plot shares its row heights with E (ef_layout() in helpers.R).

EF <- ef_layout()

draw_panel <- function() shifted_up(ROW_UP[["ef"]], {
  band_rects(EF$bxF)
  draw_at(EF$gF, F_BOX[["x"]], EF$y, F_BOX[["w"]], EF$h)
})
