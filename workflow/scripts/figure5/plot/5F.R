# Figure 5F: the CDX2 knockout's effect within each cell type (Zhu 2023, no SC-beta subtraction), knockout minus
# wild type, on G2's four EC markers, the other G2 genes and G3 in SC-SIC, SC-beta and SC-alpha, pseudobulk and per-cell
# detection, for the nuclei of the study's GEO deposit; grey bars: middle 95% of 500 matched random gene sets; dots in each
# cell type's color, filled outside that range and a ring inside it. Full page
# width; key below.

g <- panel_F()

draw_panel <- function() {
  draw_at(g, .04, FG_TOP, W - .08, FG_H)
  row_keys(list(list(type = "pt", pch = 21, label = "Outside random range", col = "grey35", fill = "grey35", size = 1.9, plwd = FG_PT_LWD),
                list(type = "pt", pch = 21, label = "Inside", col = "grey35", fill = "white", size = 1.57, plwd = .8 * ggplot2::.stroke / 2),
                list(type = "ln", label = "Middle 95% of 500 random sets", col = FG_BAR, lwd = 1.6 * ggplot2::.pt)),
           .04, W - .08, FG_TOP + FG_H + .06, gap = .12)
}
