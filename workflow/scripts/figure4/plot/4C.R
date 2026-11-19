# Figure 4C: median expression difference from primary intestinal EC (log2; dashed line) of SC-SIC, SC-beta and SC-alpha on
# the gene groups G1, G2 and G3, as an arrow from each row's start (open circle) to its end; behind each arrow a pale bar
# spans the middle half of the group's genes at the end.

source(file.path(dirname(PLOT_SCRIPT), "page.R"))
draw_panel <- function() { draw_letter("C"); draw_gg() }
