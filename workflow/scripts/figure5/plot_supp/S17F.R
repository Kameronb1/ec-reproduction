# Figure S17F: share of each transferred cell type among the CDX2 wild-type and knockout nuclei (Zhu 2023; dots, the
# five types = 100%), with the shares in Zhu 2023's published Figure 5D as dashes beside them; SC-SIC's knockout /
# wild-type ratio in its label. Key below. It shares its row with S17G.

pF <- shares_plot()

draw_panel <- function() in_row("FG", 1, function(w, vh) {
  in_vp(0, .2, w, vh - .2, function() grid.draw(ggplotGrob(pF + theme(plot.background = element_blank()))))
  ky2(list(list(type = "gp", label = "Label transfer", fill = "black", size = 1.2), list(type = "dash", label = "Published shares")), w / 2, .08)
})
