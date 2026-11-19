# Figure 2D key, below the panel: the color of each gene group and of the mixed genes.

GROUP_KEY <- c(G1 = "G1: EC genes low in SC-SICs", G2 = "G2: shared EC genes", G3 = "G3: separate SC-SIC program",
               Mixed = "Mixed: 2 of 3 intestinal studies")

draw_panel <- function() { k <- at("D_key")
  key_grid(lapply(names(GROUP_KEY), function(g) dot_item(GROUP_KEY[[g]], GROUP_COLS[[g]])), k[["x"]], k[["w"]], k[["y"]], ncol = 1) }
