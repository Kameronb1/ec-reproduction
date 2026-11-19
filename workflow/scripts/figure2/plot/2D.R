# Figure 2D: the three gene groups at their 2C position (x), by the number of intestinal studies in which each gene is
# 2-fold up in EC (rows); the mixed genes (2 of 3 studies) in the middle row. The plot is built in volcano_and_groups.R.

source(file.path(dirname(PLOT_SCRIPT), "volcano_and_groups.R"))
D_YTITLE <- "Intestinal studies with the gene\n2-fold up in EC"

draw_panel <- function() { text_at(D_YTITLE, at("D_ytitle")[["x"]], at("D_ytitle")[["y"]]); plot_at(plot_D, "D_plot") }
