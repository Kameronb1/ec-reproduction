# Figure 5A: RNA of 34 transcription factors (columns, in their display groups) by cell type (rows: native intestinal, induced and HIO graft EC,
# SC-islet SC-SIC, SC-beta and SC-alpha, adult islet beta and alpha), as mean log2(1 + CPM). Below the heatmap, SC-SIC
# minus SC-beta in each SC-islet study (dots, colored by which cell type has more) and their mean (bar). Color bar and
# dot key at the lower left. The gene names under the strip are cut where the B-D row begins; their ends below that line
# belong to 5B and 5C (helpers.R).

A <- plot_A()

# Horizontal color bar (5A's magma ramp, A_PAL), w wide with its top-left at (x, y), ticks at round values and the title above it.
color_bar <- function(x, y, w, lims, title) {
  grid.raster(matrix(A_PAL, nrow = 1), x = unit(x, "in"), y = unit(H - y, "in"), width = unit(w, "in"),
              height = unit(.075, "in"), just = c("left", "top"), interpolate = TRUE)
  br <- pretty(lims, 4); br <- br[br >= lims[1] & br <= lims[2]]
  for (b in br) {
    xx <- x + (b - lims[1]) / diff(lims) * w
    grid.lines(unit(rep(xx, 2), "in"), unit(H - y - c(.075, .105), "in"), gp = gpar(col = "black", lwd = .5))
    text_at(format(b), xx, y + .115, just = c("centre", "top"))
  }
  text_at(title, x, y - .03, just = c("left", "bottom"))
}

draw_panel <- function() {
  draw_A_plot(A$plot)
  # the strip's row label, beside it just above its 0 line (0.947 in above A's bottom)
  grid.text("SC-SIC minus SC-beta\n(log2), per study", x = unit(1.50, "in"), y = unit(H - (A_Y + A_H - 1.01), "in"),
            just = c("right", "centre"), gp = gpar(fontfamily = "Arial", fontsize = 7, fontface = "bold", lineheight = .95))
  color_bar(.30, A_Y + A_H - .19 - .18, .9, A$lims, "Mean log2(1 + CPM)")
  row_keys(list(list(type = "pt", label = "SC-SIC higher", col = POP_COL[["SC-SIC"]]),
                list(type = "dash", label = "Mean", col = "black", lwd = .2 * ggplot2::.pt * 2.5)),
           .10, 1.7, A_Y + A_H - .19 - .66, center = FALSE, gap = .12)
  row_keys(list(list(type = "pt", label = "SC-beta higher", col = POP_COL[["SC-beta"]])),
           .10, 1.7, A_Y + A_H - .19 - .54, center = FALSE, gap = .12)
  row_keys(list(list(type = "pt", label = "No RNA at this depth", col = "grey60")),
           .10, 1.7, A_Y + A_H - .19 - .42, center = FALSE, gap = .12)
  white_outside(0, 0, W, ROW_BCD_TOP)               # the gene names' ends below the B-D row's top belong to 5B and 5C
}
