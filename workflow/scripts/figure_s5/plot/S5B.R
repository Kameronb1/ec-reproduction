# Figure S5B: MSigDB C8 cell-type signature terms of gene groups G1, G2 and G3, each tested alone and with its mixed genes;
# laid out as S5A, so the two plots line up. The panels are built in enrichment_panels.R.

source(file.path(dirname(PLOT_SCRIPT), "enrichment_panels.R"))

draw_panel <- function() draw_at(PB$g, X0, yB, PW, PB$h)
