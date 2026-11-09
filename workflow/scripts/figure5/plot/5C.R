# Figure 5C: promoter accessibility across stage 6 culture (weeks 2-4, months 6 and 12) and in grafts (median of three;
# Augsornworawat 2023) within SC-SIC, SC-beta and SC-alpha: each gene group against the typical gene of the same cell type,
# as its median gene (line and dots; the graft median a hollow diamond) and the middle half of its genes (light bar), on a
# fold axis as deep as 5B's. The plot reaches up beside the letter C. As for 5B, this panel's image also carries the ends
# of panel A's gene names above it (helpers.R).

A <- plot_A()
g <- row_BCD()

draw_panel <- function() {
  draw_A_plot(A$plot)
  draw_BCD(g, "C")
  white_outside(g$cut[["BC"]], ROW_BCD_TOP, g$cut[["CD"]], H)   # ink above the B-D row's top belongs to 5A; beside the cuts to 5B and 5D
}
