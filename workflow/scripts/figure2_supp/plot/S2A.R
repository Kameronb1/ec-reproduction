# Figure S2A: UMAP of the 1,199 fitting cells, built from Harmony-corrected expression with the four score genes
# withheld: colored by cluster (numbered as in S2B) and by four-marker class, then one small map per gene colored
# by its expression as a fraction of the gene's maximum. The expression key fills the last tile slot.

cc <- read_table("cell_clusters.tsv.gz")
u <- read_table("umap_coordinates.tsv.gz")
GENES <- c("TPH1", "SLC18A1", "FEV", "LMX1A", "CHGA", "NEUROG3", "DDC", "TAC1", "CCK", "GCG", "SST", "GHRL", "NTS", "OLFM4", "MKI67")
ex <- read_table("gene_expression.tsv")                         # log1p CPM of the 15 genes
E <- as.matrix(ex[GENES]); rownames(E) <- ex$cell_key
mx <- apply(E, 2, max)
RAMP_HI <- MARKS[["4"]]

v <- u[u$reduction == "Harmony, score genes withheld", ]
v$cluster <- cc$cluster[match(v$cell_key, cc$cell_key)]
cen <- do.call(rbind, lapply(split(v, v$cluster), function(z) {      # each cluster's number sits on its cell nearest the median
  m <- c(median(z$disp_1), median(z$disp_2)); i <- which.min((z$disp_1 - m[1])^2 + (z$disp_2 - m[2])^2); z[i, c("cluster", "disp_1", "disp_2")] }))
xl <- range(v$disp_1); yl <- range(v$disp_2)
rp_dat <- v[, c("cell_key", "disp_1", "disp_2")]; rp_dat$lab <- ""    # every cell repels the numbers; 13 cells carry one
rp_dat$lab[match(paste(cen$disp_1, cen$disp_2), paste(rp_dat$disp_1, rp_dat$disp_2))] <- as.character(cen$cluster)
umap_theme <- function() theme_kit("none") + theme(axis.text = element_blank(), axis.title = element_blank(), axis.line.x = element_blank(),
                                                   axis.ticks.x = element_blank(), aspect.ratio = 1, panel.grid = element_blank(),
                                                   plot.margin = margin(0, 0, 0, 0))
cm <- ggplot(v[order(match(v$presentation_class, rev(names(CALLS)))), ], aes(disp_1, disp_2, color = presentation_class)) +
  geom_point(size = .6, stroke = 0) + scale_color_manual(values = CALLS) + umap_theme()   # high-confidence EC drawn last
v$cl <- factor(as.character(v$cluster), levels = names(CLUSTER_COL))
set.seed(1); vk <- v[sample(nrow(v)), ]                                 # random drawing order: no cluster painted over another
km <- ggplot(vk, aes(disp_1, disp_2)) + geom_point(aes(color = cl), size = .6, stroke = 0) + scale_color_manual(values = CLUSTER_COL) +
  ggrepel::geom_text_repel(data = rp_dat, aes(label = lab), color = "black", size = tx, family = FONT, fontface = "bold",
                           seed = 6, box.padding = .08, point.padding = .02, min.segment.length = .05, segment.size = .25,
                           max.overlaps = Inf, force = 2, max.iter = 20000) +
  umap_theme()
tile <- function(g) {                                                  # highest expression drawn last
  d <- data.frame(disp_1 = v$disp_1, disp_2 = v$disp_2, rel = E[v$cell_key, g] / mx[[g]]); d <- d[order(d$rel), ]
  ggplot(d, aes(disp_1, disp_2, color = rel)) + geom_point(size = .42, stroke = 0) +
    scale_color_gradient(low = "#E6E6E6", high = RAMP_HI, limits = c(0, 1)) +
    scale_x_continuous(limits = xl) + scale_y_continuous(limits = yl) + umap_theme()
}
tiles <- lapply(GENES, tile)
gp <- function(face = "plain") gpar(fontfamily = FONT, fontsize = PT, fontface = face, col = "black")
# map size and left edges; gene tiles in 8 columns between GX0 and GX1, each TILE square under a STRIP for its name (in)
MAP_S <- 2.20; MAP1_X <- .06; MAP2_X <- MAP1_X + MAP_S + .22; GX0 <- .06; GX1 <- 7.26; NC <- 8; TILE <- .87; STRIP <- .15
W_A <- 7.3; H_A <- 5.26 - .60                                          # content width and height (in)
map_vp <- function(x, top) viewport(x = unit(x, "in"), y = unit(top - .22, "in"), width = unit(MAP_S, "in"), height = unit(MAP_S, "in"),
                                    just = c("left", "top"))
# ggrepel places the cluster numbers by a force simulation run when the map is drawn. The simulation amplifies
# rounding in the last bits of the device coordinates, which depend on where the map sits on the device, so the
# numbers are placed with the panel filling a device of its own size (4380 x 2832 px), the geometry of the published
# panel; the placed map is then drawn on the page.
km_placed <- local({
  page <- dev.cur()
  png(tempfile(fileext = ".png"), width = 4380, height = 2832, units = "px", res = DPI, type = "cairo", bg = "white")
  pushViewport(viewport(x = unit(0, "in"), y = unit(0, "in"), width = unit(W_A, "in"), height = unit(H_A, "in"), just = c("left", "bottom")))
  pushViewport(map_vp(MAP1_X, convertHeight(unit(1, "npc"), "in", valueOnly = TRUE)))
  g <- forceGrob(ggplotGrob(km))
  invisible(dev.off()); dev.set(page)
  g
})
draw_a <- function() {
  top <- convertHeight(unit(1, "npc"), "in", valueOnly = TRUE)
  grid.text("13 clusters, numbered as in B", unit(MAP1_X + MAP_S / 2, "in"), unit(top - .04, "in"), just = c("centre", "top"), gp = gp("bold"))
  grid.text("Four-marker classifier calls", unit(MAP2_X + MAP_S / 2, "in"), unit(top - .04, "in"), just = c("centre", "top"), gp = gp("bold"))
  pushViewport(map_vp(MAP1_X, top)); grid.draw(km_placed); popViewport()
  pushViewport(map_vp(MAP2_X, top)); print(cm, newpage = FALSE); popViewport()
  kx <- MAP2_X + MAP_S + .25; ky <- top - .22 - MAP_S / 2 + .17          # class key, right of the class map
  for (i in seq_along(CALLS)) {
    yy <- ky - (i - 1) * .17
    grid.points(unit(kx, "in"), unit(yy, "in"), pch = 16, size = unit(5, "pt"), gp = gpar(col = CALLS[[i]]))
    grid.text(CLASS_KEY[[names(CALLS)[i]]], unit(kx + .12, "in"), unit(yy, "in"), just = "left", gp = gp())
  }
  # gene tiles: 2 rows x 8; slot 16 holds the expression key
  g_top <- top - .22 - MAP_S - .06
  cp <- (GX1 - GX0) / NC; rp <- TILE + STRIP + .04
  slot <- function(k) { r <- (k - 1) %/% NC; c <- (k - 1) %% NC; list(x = GX0 + c * cp + (cp - TILE) / 2, y = g_top - r * rp) }
  for (k in seq_along(GENES)) {
    s <- slot(k)
    grid.text(GENES[k], unit(s$x + TILE / 2, "in"), unit(s$y - .04, "in"), just = c("centre", "top"), gp = gp("italic"))
    pushViewport(viewport(x = unit(s$x, "in"), y = unit(s$y - STRIP - .04, "in"), width = unit(TILE, "in"), height = unit(TILE, "in"), just = c("left", "top")))
    print(tiles[[k]], newpage = FALSE); popViewport()
  }
  s <- slot(16); kl <- c("Expression,", "fraction of each", "gene's maximum"); bw <- .78; bx <- s$x + (TILE - bw) / 2
  y0 <- s$y - .22
  for (i in seq_along(kl)) grid.text(kl[i], unit(s$x + TILE / 2, "in"), unit(y0 - (i - 1) * .125, "in"), just = c("centre", "top"), gp = gp())
  yb <- y0 - 3 * .125 - .08
  grid.raster(matrix(grDevices::colorRampPalette(c("#E6E6E6", RAMP_HI))(200), nrow = 1), x = unit(bx, "in"), y = unit(yb, "in"),
              width = unit(bw, "in"), height = unit(.08, "in"), just = c("left", "top"), interpolate = FALSE)
  for (t in c(0, .5, 1)) grid.text(format(t), unit(bx + bw * t, "in"), unit(yb - .11, "in"), just = c("centre", "top"), gp = gp())
}

draw_panel <- function() place(draw_a, 120, 48, 4380, 2832, W_A, H_A)
