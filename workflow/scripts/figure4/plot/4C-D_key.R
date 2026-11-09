# Figure 4C-D key, under both panels. Line 1: C's ridges (start cells: outline; graft or continued-culture cells: filled),
# the cutoff for classifying a cell as EC, then the population colors. Line 2: the four markers, and how D draws genes
# and their median.

KEY1 <- 6.42; KEY2 <- 6.58                                     # tops of the two key lines
# a small ridge glyph (.20 x .075 in, its foot .035 in under the glyph line), then its text; returns the next x
ridge_key <- function(x, y, filled, label) {
  t <- seq(0, 1, length.out = 40); yb <- H - (y + .045 + .035)
  grid.polygon(unit(x + .2 * c(t, rev(t)), "in"), unit(yb + c(.075 * exp(-((t - .55) / .2)^2), rep(0, 40)), "in"),
               gp = if (filled) gpar(fill = "#888888", alpha = .55, col = NA) else gpar(fill = NA, col = "#333333", lwd = .22 * ggplot2::.pt))
  text_at(label, x + .26, y - .005); x + .26 + str_in(label) + .20 }

draw_panel <- function() shifted_up(ROW_UP[["cd"]], {
  x <- ridge_key(.10, KEY1, FALSE, "Start cells"); x <- ridge_key(x, KEY1, TRUE, "Graft / continued-culture cells")
  key_line(x, KEY1, "#555555", lwd = .3 * ggplot2::.pt, lty = "22", len = .16, butt = TRUE); text_at("Call cutoff", x + .21, KEY1 - .005)
  population_key(x + .21 + str_in("Call cutoff") + .30, KEY1)
  y <- KEY2 - .005; lead <- "Four markers: "; genes <- "TPH1, SLC18A1, FEV, LMX1A"   # gene names in italic
  text_at(lead, .10, y); xg <- .10 + str_in(lead)
  text_at(genes, xg, y, face = "italic"); text_at("D: dots = genes; black tick = median", xg + str_in(genes, face = "italic") + .30, y)
})
