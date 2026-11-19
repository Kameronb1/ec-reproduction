# Figure S2D: one row per cluster of S2C: its fitting groups, its four-marker classes and the reference genes marking
# it; keys for C and D below.
source(file.path(dirname(PLOT_SCRIPT), "cd_data.R"))
share <- function(col, levels) unclass(table(factor(cc$cluster, levels = summ$cluster), factor(cc[[col]], levels = levels))) / summ$n_cells
SG <- share("arm", names(GROUP_COL)); SC <- share("presentation_class", rev(ORD))
stopifnot(all(abs(rowSums(SG) - 1) < 1e-12), all(abs(rowSums(SC) - 1) < 1e-12))

PITCH_D <- .155
draw_D <- function(x0, y0, w) {
  name_w <- max(vapply(sprintf("%s (%d)", summ$ident, summ$n_cells), str_w, 0)) + .32
  gene_w <- max(vapply(summ$genes, str_w, 0, face = "italic"))
  GAP_D <- .32
  bar_w <- min(1.5, (w - name_w - gene_w - .1 - 2 * GAP_D) / 2)
  xg <- x0 + name_w; xc <- xg + bar_w + GAP_D; xt <- xc + bar_w + GAP_D
  txt("Cluster (cells)", x0, y0, face = "bold"); txt("Fitting group", xg, y0, face = "bold")
  txt("Classification", xc, y0, face = "bold"); txt("Reference genes marking the cluster", xt, y0, face = "bold")
  yt <- y0 + LINE + .05
  for (i in seq_len(NCL)) {
    yc <- yt + (i - .5) * PITCH_D
    if (i %% 2 == 1) rect_at(x0, yc - PITCH_D / 2, w, PITCH_D, fill = "#F2F2F2")
    txt(summ$cluster[i], x0 + .12, yc, just = c("right", "centre"))
    grid.rect(unit(x0 + .17, "in"), unit(H - yc, "in"), unit(.08, "in"), unit(.08, "in"), gp = gpar(fill = CL_COL[[as.character(summ$cluster[i])]], col = NA))
    txt(sprintf("%s (%d)", summ$ident[i], summ$n_cells[i]), x0 + .25, yc, just = c("left", "centre"))
    for (sp in list(list(m = SG, x = xg, cols = GROUP_COL), list(m = SC, x = xc, cols = CLS_COL[rev(ORD)]))) {
      x1 <- cumsum(sp$m[i, ]); x0s <- x1 - sp$m[i, ]
      for (j in seq_along(x1)) if (sp$m[i, j] > 0) rect_at(sp$x + bar_w * x0s[j], yc - .05, bar_w * sp$m[i, j], .1, fill = sp$cols[[j]])
    }
    grid.text(summ$genes[i], unit(xt, "in"), unit(H - yc, "in"), just = c("left", "centre"), gp = gp_t(face = "italic"))
  }
  yb <- yt + NCL * PITCH_D
  for (xx in c(xg, xc)) { seg_at(xx, yb, xx + bar_w, yb, lwd = .6)
    for (v0 in c(0, 50, 100)) { seg_at(xx + bar_w * v0 / 100, yb, xx + bar_w * v0 / 100, yb + .03, lwd = .6)
      txt(paste0(v0, "%"), xx + bar_w * v0 / 100, yb + .045, just = c("centre", "top")) } }
  yb + .045 + LINE
}

draw_panel <- function() {
  letter("D", L - .1, S2_D_TOP)
  bottom <- draw_D(L + .15, S2_D_TOP + .02, R - L - .15)
  # shared key for C and D: fitting groups (D) and the three classes (C, D)
  yk <- bottom + .16
  key_row(c(list(k_text("Fitting group", "bold")), lapply(names(GROUP_COL), function(a) k_square(GROUP_LAB[[a]], GROUP_COL[[a]]))), L, yk, width = R - L)
  key_row(c(list(k_text("Classification", "bold")), lapply(ORD[3:1], function(a) k_square(CLS_LAB[[a]], CLS_COL[[a]]))), L, yk + 1.6 * LINE, width = R - L)
}
