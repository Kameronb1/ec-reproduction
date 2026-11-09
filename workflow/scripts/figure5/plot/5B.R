# Figure 5B: promoter accessibility over Zhu 2023's differentiation days 21, 32 and 39 within SC-SIC, SC-beta and
# SC-alpha: each gene group against the typical gene of the same cell type, as its median gene (line and dots) and the
# middle half of its genes (light bar), on 5C's fold axis. B, C and D are built together so that their plot areas line
# up (helpers.R). This panel's image also carries the ends of panel A's gene names above it.

A <- plot_A()
g <- row_BCD()

draw_panel <- function() {
  draw_A_plot(A$plot)
  draw_BCD(g, "B")
  white_outside(0, ROW_BCD_TOP, g$cut[["BC"]], H)   # ink above the B-D row's top belongs to 5A; right of the cut to 5C
}
