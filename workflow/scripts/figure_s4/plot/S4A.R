# Figure S4A: the 34 G2 genes (rows) in each study (columns), log2(1 + CPM) of counts pooled per study: the nine genes
# higher in SC-SICs, then the 25 without a significant difference; with the color key and the column-group key.
draw_panel <- function() {
  reset_metric_cache()
  letter_at("A", MARGIN + .02, A_HDR - .04); group_head("G2: shared EC genes (34)", GROUP_COL[["G2"]], A_HDR - .03)
  for (k in 1:2) { labels_up(k, A_BAR - .04); bars(k, A_BAR) }
  grid.text("Higher in SC-SICs, FDR < 0.05 (9)", XU(SX[1]), YU(A_DATA + SECT_H / 2), just = c("left", "centre"), gp = gp(PT, "bold"))
  cells(1, A_ROW1, G2A)
  y2 <- A_ROW1 + 9 * R + .06
  grid.text("No significant difference (25)", XU(SX[1]), YU(y2 + SECT_H / 2), just = c("left", "centre"), gp = gp(PT, "bold"))
  cells(1, y2 + SECT_H, G2B[1:8]); cells(2, A_ROW1, G2B[9:25])
  grid.text("No significant difference (continued)", XU(SX[2]), YU(A_DATA + SECT_H / 2), just = c("left", "centre"), gp = gp(PT, "bold"))
  # the key, in the third slot beside the G2 blocks
  kx <- SX[3] + .10; ky <- A_ROW1
  grid.text(expression("Expression, log"[2]*"(1 + CPM),"), XU(kx), YU(ky), just = c("left", "top"), gp = gp())
  grid.text("counts pooled per study", XU(kx), YU(ky + .15), just = c("left", "top"), gp = gp())
  bw <- 1.4; yb <- ky + .36
  ramp_bar(lv_col(seq(LV_RANGE[1], LV_RANGE[2], length.out = 256)), kx, yb, bw)
  br <- seq(0, LV_RANGE[2], 2); bx <- kx + (br - LV_RANGE[1]) / diff(LV_RANGE) * bw
  grid.segments(XU(bx), YU(yb + .085), XU(bx), YU(yb + .085 - .025), gp = gpar(col = "white", lwd = .8))
  grid.text(br, XU(bx), YU(yb + .085 + .035), just = c("centre", "top"), gp = gp())
  its <- lapply(unique(COLS$pop), function(p) sq_item(p, POP_COL[[p]]))
  for (q in seq_along(its)) { yq <- yb + .40 + (q - 1) * .16; its[[q]]$draw(kx, yq)
    grid.text(its[[q]]$label, XU(kx + its[[q]]$w_sym), YU(yq + .0425), just = c("left", "centre"), gp = gp()) }
}
