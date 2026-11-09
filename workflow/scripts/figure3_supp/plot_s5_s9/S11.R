# Figure S11: gene group G3, SC-SIC's own program, in two side-by-side blocks (strongest Figure 3D evidence on top,
# reading down the left block, then the right block), with the columns and keys of S10. "+" marks a gene 2-fold up in
# EC in that intestinal study. Column labels stand upright: the columns are 0.2 in apart here.

draw_panel <- function() {
  t7 <- T7(); rng <- LV_RANGE()
  g3l <- t7[t7$block == "G3 left", ]; g3r <- t7[t7$block == "G3 right", ]
  gl_rel <- max(str_in(c(g3l$gene, g3r$gene), PT_GENE, "italic")); GB <- .16   # widest gene name; gap between the blocks
  bw <- (PW - .05 - .04 - GB) / 2                                               # block width
  p <- (bw - gl_rel - .05 - 2 * .07 - .12) / 14                                 # column pitch
  lft <- c(.04 + gl_rel + .05, .04 + bw + GB + gl_rel + .05)                    # left edge of each block's first column
  xs <- list(col_x(lft[1], p, .07, .12), col_x(lft[2], p, .07, .12))
  blocks <- list(g3l, g3r)
  R <- .105
  y <- 0
  Y <- list(h = y); y <- y + HEADER_H
  Y$strip <- y + .15 + LINE_IN(PT_STRIP); y <- Y$strip + .04
  Y$b <- y; y <- y + nrow(g3l) * R
  Y$bot <- y; y <- y + label_depth(90) + .08
  Y$key <- y; y <- y + KEY_H
  H7 <- ceiling((y + .04) * 600) / 600                                          # whole pixels at 600 dpi

  heatmap_canvas(60, 48, H7, function() {
    block_heading(Y$h, GROUP_COL[["G3"]], sprintf("G3: SC-SIC's own program (%d genes)", nrow(g3l) + nrow(g3r)))
    for (s in 1:2) {
      draw_strips(Y$strip, xs[[s]], p, "EC vs\nother EEC")
      draw_block(cells(blocks[[s]], TRUE, rng), blocks[[s]]$gene, Y$b, xs[[s]], p, R, lft[s] - .05)
      draw_col_labels(Y$bot, xs[[s]], 90)
    }
    draw_keys(Y$key, lft[1], 2.55, 5.25, TRUE, rng)
  })
}
