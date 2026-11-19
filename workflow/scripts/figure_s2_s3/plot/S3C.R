# Figure S3C: percent of SC-SICs and of other SC endocrine cells classified as EC (filled) and above the cutoff
# (outline) in each SC-islet study.
g0 <- read_table("classified_by_study.tsv")
ST_G <- c("Augsornworawat 2020", "Balboa 2022", "Du 2022", "Augsornworawat 2023 (nuclei)", "Zhu 2023", "Siehler 2024")
g0 <- g0[g0$class %in% c("SC_SIC", "OTHER_MATURE_SC_ENDOCRINE"), ]
stopifnot(setequal(g0$study, ST_G), nrow(g0) == 12L, all(abs(g0$pct - round(100 * g0$classified_ec / g0$n, 2)) < .006),
          all(abs(g0$pct_candidates - round(100 * g0$candidates / g0$n, 2)) < .006), all(g0$classified_ec <= g0$candidates))
PITCH_G <- .17
draw_G <- function(x0, y0, w) {
  lw <- max(vapply(ST_G, str_w, 0)) + .08
  gap <- .35; aw <- (w - lw - gap) / 2; XMAX <- c(SC_SIC = 100, OTHER_MATURE_SC_ENDOCRINE = 15)
  XL <- c(SC_SIC = x0 + lw, OTHER_MATURE_SC_ENDOCRINE = x0 + lw + aw + gap)
  txt("SC-SICs", XL[["SC_SIC"]], y0, face = "bold"); txt("Other SC endocrine cells", XL[["OTHER_MATURE_SC_ENDOCRINE"]], y0, face = "bold")
  yt <- y0 + LINE + .05
  for (k in names(XL)) for (v0 in seq(0, XMAX[[k]], by = if (k == "SC_SIC") 50 else 5))
    seg_at(XL[[k]] + aw * v0 / XMAX[[k]], yt, XL[[k]] + aw * v0 / XMAX[[k]], yt + 6 * PITCH_G, col = "#DADADA", lwd = .4)
  for (i in seq_along(ST_G)) {
    yc <- yt + (i - .5) * PITCH_G
    if (i %% 2 == 1) rect_at(x0, yc - PITCH_G / 2, w, PITCH_G, fill = "#F2F2F2")
    txt(ST_G[i], x0 + lw - .06, yc, just = c("right", "centre"))
    for (k in names(XL)) {
      r <- g0[g0$study == ST_G[i] & g0$class == k, ]; col <- if (k == "SC_SIC") SC_SIC else OTHER_SC
      rect_at(XL[[k]], yc - .045, aw * r$pct_candidates / XMAX[[k]], .09, fill = "white", col = col, lwd = .6)
      rect_at(XL[[k]], yc - .045, aw * r$pct / XMAX[[k]], .09, fill = col, col = col, lwd = .6)
      txt(sprintf("%.0f%%", r$pct), XL[[k]] + aw * r$pct_candidates / XMAX[[k]] + .04, yc, just = c("left", "centre"))
    }
  }
  yb <- yt + 6 * PITCH_G
  for (k in names(XL)) { seg_at(XL[[k]], yb, XL[[k]] + aw, yb, lwd = .6)
    for (v0 in seq(0, XMAX[[k]], by = if (k == "SC_SIC") 50 else 5)) {
      seg_at(XL[[k]] + aw * v0 / XMAX[[k]], yb, XL[[k]] + aw * v0 / XMAX[[k]], yb + .03, lwd = .6)
      txt(paste0(v0, "%"), XL[[k]] + aw * v0 / XMAX[[k]], yb + .045, just = c("centre", "top")) } }
  txt("Cells of each population (%)", x0 + lw + (w - lw) / 2, yb + .045 + 1.15 * LINE, just = c("centre", "top"))
  yk <- yb + .045 + 2.45 * LINE + .07
  key_row(list(k_multi("Classified as EC", c(SC_SIC, OTHER_SC))), x0, yk, width = w)
  key_row(list(list(label = "Above cutoff, not classified", w_sym = .27, draw = function(x, y) {
    grid.rect(unit(x, "in"), unit(H - y, "in"), unit(.09, "in"), unit(.09, "in"), just = c("left", "centre"), gp = gpar(fill = "white", col = SC_SIC, lwd = .6))
    grid.rect(unit(x + .11, "in"), unit(H - y, "in"), unit(.09, "in"), unit(.09, "in"), just = c("left", "centre"), gp = gpar(fill = "white", col = OTHER_SC, lwd = .6)) })),
    x0, yk + 1.5 * LINE, width = w)
  yk + 1.5 * LINE + .07
}

draw_panel <- function() {
  letter("C", L - .1, S3_C_TOP)
  invisible(draw_G(L + .15, S3_C_TOP + .02, 5.3))
}
