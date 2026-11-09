# Figure S16A: cholesterol synthesis genes in grafts and over in vitro culture time alone, and the NFE2L1 regulon target
# score. One block per comparison (SC-SIC, SC-beta and SC-alpha lines): a heatmap of the 26 genes' change (log2; colors
# end at +-3, larger changes drawn at the ends; a dot where every earlier x later sample pair changes the same way by at
# least 0.5 on average), each line's median change as a bar with its value, and the change of the NFE2L1 regulon target
# score (AUCell, test minus reference) as a second bar column. Rows: eight grafts, then five longer in vitro cultures.

source(file.path(dirname(PLOT_SCRIPT), "s2_page.R"))
draw_panel <- function() s2_panel("A")
