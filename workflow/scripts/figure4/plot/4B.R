# Figure 4B: for each of the eleven rows, the share of SC-SICs and SC-beta cells with a four-marker score at or above the
# cutoff, from the start (open circle) to the graft or extended culture (arrowhead), samples weighted equally (left), and
# the four EC markers' mean expression change (right); label table of study, line and end point and samples start -> end.

source(file.path(dirname(PLOT_SCRIPT), "page.R"))
draw_panel <- function() { draw_letter("B"); draw_fm() }
