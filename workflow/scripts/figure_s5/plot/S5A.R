# Figure S5A: pathway terms (GO biological process, Reactome) of gene groups G1, G2 and G3, each tested alone and with its
# mixed genes. The panels are built in enrichment_panels.R.

source(file.path(dirname(PLOT_SCRIPT), "enrichment_panels.R"))

draw_panel <- function() draw_at(PA$g, X0, yA, PW, PA$h)
