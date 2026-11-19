# Figure 2C: SC-SIC against intestinal EC for every tested gene (donor-level limma-voom fit, study as a block). The plot is
# built in volcano_and_groups.R.

source(file.path(dirname(PLOT_SCRIPT), "volcano_and_groups.R"))

draw_panel <- function() plot_at(plot_C, "C_plot")
