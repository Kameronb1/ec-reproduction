# Figure S13A: pathway terms (GO biological process, Reactome) of gene groups G3, G2 and G1, each tested alone (open dot)
# and with its mixed genes (filled dot; mixed genes: 2-fold up in EC in 2 of the 3 intestinal studies, Figure 3F), as -log10
# adjusted p, with the mixed genes in each term at the right.

draw_panel <- function() {
  all <- read_table("s8_mixed/pathways/enrichment_all_terms.tsv.gz")
  P <- term_panel(s8_terms("A", all, paste0(all$collection, ": ", all$term)), 6)
  plot_panel(P$plot, 120, 48, PW, P$h)
}
