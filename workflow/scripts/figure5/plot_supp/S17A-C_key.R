# Figure S17A-C key, two centered rows under A, B and C: each gene group G1-G3 as B and C draw it (its point on its light
# bar); then A's sample symbols (Augsornworawat 2023 circles, Zhu 2023 triangles), B's and C's group median gene and the
# bar for the middle half of the group's genes, in a dotted outline.

GL <- fig5_env()$GL

draw_panel <- function() in_row("ABC", 1, function(w, vh) {
  w1 <- ky2(lapply(seq_along(GCOL), function(k) list(type = "pb", label = GL[k], fill = GCOL[[k]], bar = GCOL_LIGHT[[k]])), w / 2, .25)
  w2 <- ky2(list(list(type = "gp", pch = 21, label = "Augsornworawat 2023 sample", fill = "white"), list(type = "gp", pch = 24, label = "Zhu 2023 sample", fill = "white"),
           list(type = "gp", label = "Group's median gene", fill = "grey75"), list(type = "bars3", label = "Middle half of the group's genes")), w / 2, .09)
  grid.rect(unit(w / 2, "in"), unit(.1775, "in"), width = unit(max(w1, w2) + .2, "in"), height = unit(.315, "in"),
            gp = gpar(fill = NA, col = "grey45", lwd = .6, lty = "dotted"))
})
