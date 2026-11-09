# Figure S13B: MSigDB C8 cell-type signature terms of gene groups G3, G2 and G1, each tested alone (open dot) and with its
# mixed genes (filled dot), as -log10 adjusted p, with the mixed genes in each term at the right; laid out as S13A, so
# the two plots line up.

draw_panel <- function() {
  all <- read_table("s8_mixed/celltypes/celltype_all_terms.tsv")
  B <- s8_terms("B", all, all$term)
  G3G <- c("midbrain neuroblasts" = "embryonic midbrain neuroblasts", "colon enteroendocrine" = "adult colon enteroendocrine",
           "fetal gut chromaffin" = "fetal intestine chromaffin")   # Descartes fetal intestine term, as in Figure 3G
  stopifnot(all(names(G3G) %in% B$term)); B$term <- ifelse(B$term %in% names(G3G), G3G[B$term], B$term)   # Figure 3G's names
  P <- term_panel(B, 3)
  plot_panel(P$plot, 120, 3240, PW, P$h)        # below S13A's canvas (3132 pixels from page pixel 48) and a 60-pixel gap
}
