# Figure S2C: marker-free map of the fitting cells, colored by cluster (numbers as in S2D) and by four-marker class.
source(file.path(dirname(PLOT_SCRIPT), "cd_data.R"))
cen <- do.call(rbind, lapply(split(v, v$cluster), function(z) {
  m <- c(median(z$disp_1), median(z$disp_2)); i <- which.min((z$disp_1 - m[1])^2 + (z$disp_2 - m[2])^2); z[i, c("cluster", "disp_1", "disp_2")] }))
umap_theme <- theme_void() + theme(legend.position = "none", plot.margin = margin(0, 0, 0, 0))
set.seed(1); vk <- v[sample(nrow(v)), ]                    # drawing order of the cells
xr <- range(v$disp_1); yr <- range(v$disp_2); padx <- .07 * diff(xr); pady <- .07 * diff(yr)
# cluster numbers: placed in map units scaled to 0-1, each started beside its anchor cell (away from the map centre),
# then pushed apart from other numbers and off cells of other clusters, and kept inside the map; a short line joins a
# number to its anchor when the number moved away from it
place_labels <- function(cen, v, min_d = .09, cell_d = .055, iters = 1500) {
  sx <- function(x) (x - xr[1] + padx) / (diff(xr) + 2 * padx); sy <- function(y) (y - yr[1] + pady) / (diff(yr) + 2 * pady)
  ax <- sx(cen$disp_1); ay <- sy(cen$disp_2); px <- sx(v$disp_1); py <- sy(v$disp_2)
  dx <- ax - .5; dy <- ay - .5; l <- pmax(sqrt(dx^2 + dy^2), 1e-6)
  lx <- ax + .08 * dx / l; ly <- ay + .08 * dy / l
  for (it in seq_len(iters)) {
    mx <- 0 * lx; my <- 0 * ly
    for (i in seq_along(lx)) {
      for (j in seq_along(lx)) if (i != j) {                      # numbers apart (labels about 0.07 wide, 0.06 tall)
        ex <- lx[i] - lx[j]; ey <- (ly[i] - ly[j]) * 1.25; d <- sqrt(ex^2 + ey^2)
        if (d < min_d) { mx[i] <- mx[i] + (min_d - d) * ex / max(d, 1e-6) / 2; my[i] <- my[i] + (min_d - d) * ey / max(d, 1e-6) / 2 }
      }
      ex <- lx[i] - px; ey <- (ly[i] - py) * 1.25; d <- sqrt(ex^2 + ey^2); k <- d < cell_d      # off every cell
      if (any(k)) { vx <- sum(ex[k] / pmax(d[k], 1e-6)); vy <- sum(ey[k] / pmax(d[k], 1e-6)); vn <- max(sqrt(vx^2 + vy^2), 1e-6)
        mx[i] <- mx[i] + .004 * vx / vn; my[i] <- my[i] + .004 * vy / vn }
      mx[i] <- mx[i] + .003 * (ax[i] - lx[i]); my[i] <- my[i] + .003 * (ay[i] - ly[i])   # weak pull to the anchor
    }
    lx <- pmin(pmax(lx + mx, .05), .95); ly <- pmin(pmax(ly + my, .05), .95)
  }
  ux <- function(s) s * (diff(xr) + 2 * padx) + xr[1] - padx; uy <- function(s) s * (diff(yr) + 2 * pady) + yr[1] - pady
  data.frame(cluster = cen$cluster, ax = cen$disp_1, ay = cen$disp_2, lx = ux(lx), ly = uy(ly),
             far = sqrt((lx - ax)^2 + (ly - ay)^2) > .045)
}
LBL_C <- place_labels(cen, v)
seg_C <- LBL_C[LBL_C$far, ]
# the segment stops short of the number: it starts 55% of the way from the number to the anchor
seg_C$sx <- seg_C$lx + .55 * (seg_C$ax - seg_C$lx); seg_C$sy <- seg_C$ly + .55 * (seg_C$ay - seg_C$ly)
plot_C1 <- ggplot(vk, aes(disp_1, disp_2)) + geom_point(aes(colour = factor(cluster)), size = .45, stroke = 0) +
  scale_colour_manual(values = CL_COL) +
  geom_segment(data = seg_C, aes(x = sx, y = sy, xend = ax, yend = ay), inherit.aes = FALSE, linewidth = .25) +
  geom_text(data = LBL_C, aes(lx, ly, label = cluster), inherit.aes = FALSE, size = MM, family = FONT, fontface = "bold") +
  coord_cartesian(xlim = xr + c(-padx, padx), ylim = yr + c(-pady, pady), expand = FALSE) + umap_theme
plot_C2 <- ggplot(v[order(match(v$presentation_class, ORD)), ], aes(disp_1, disp_2, colour = presentation_class)) +
  geom_point(size = .45, stroke = 0) + scale_colour_manual(values = CLS_COL) +
  coord_cartesian(xlim = xr + c(-padx, padx), ylim = yr + c(-pady, pady), expand = FALSE) + umap_theme

draw_C <- function(x0, y0, s, gap = .12) {     # two square maps side by side; s = side (in)
  ym <- y0 + LINE + .03
  for (k in 1:2) {
    xk <- x0 + (k - 1) * (s + gap)
    txt(c("Clusters", "Classification")[k], xk + s / 2, y0, just = c("centre", "top"), face = "bold")
    rect_at(xk, ym, s, s, fill = NA, col = "#BDBDBD", lwd = .4)
    draw_at(list(plot_C1, plot_C2)[[k]], xk, ym, s, s)
  }
  ym + s
}

draw_panel <- function() {
  letter("C", L - .1, S2_C_TOP)
  stopifnot(draw_C(L + .02, S2_C_TOP + .02, 1.13, gap = .1) == S2_C_BOTTOM)
}
