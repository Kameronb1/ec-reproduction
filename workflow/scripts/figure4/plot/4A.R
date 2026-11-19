# Figure 4A: share of mature endocrine cells (published labels) in each graft, as an arrow from the mean of its in vitro
# comparison samples (open circle; grey bar = their range) to the graft (arrowhead; open when the graft has fewer than 50
# mature endocrine cells), per cell type; sample counts start -> end beside each study.

source(file.path(dirname(PLOT_SCRIPT), "page.R"))
draw_panel <- function() { draw_letter("A"); draw_prop() }
