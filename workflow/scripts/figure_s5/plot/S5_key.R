# Figure S5 key, shared by S5A and S5B, below S5B: the three point kinds, the adjusted p = 0.05 line, and bold and plain
# term names (bold: also drawn in Figure 2E).

source(file.path(dirname(PLOT_SCRIPT), "enrichment_panels.R"))
KEY_GLYPH <- "grey35"
ring_glyph <- function(label) list(label = label, w_sym = .16, draw = function(x, y) {
  grid.points(unit(x + .05, "in"), unit(H - y - .04, "in"), pch = 1, size = unit(1, "char"), gp = gpar(col = "black", fontsize = RING$under * .pt + RING$under_stroke * .stroke / 2, lwd = RING$under_stroke * .stroke / 2))
  grid.points(unit(x + .05, "in"), unit(H - y - .04, "in"), pch = 1, size = unit(1, "char"), gp = gpar(col = KEY_GLYPH, fontsize = RING$size * .pt + RING$stroke * .stroke / 2, lwd = RING$stroke * .stroke / 2)) })
dot_glyph <- function(label) list(label = label, w_sym = .16, draw = function(x, y)
  grid.points(unit(x + .05, "in"), unit(H - y - .04, "in"), pch = 21, size = unit(1, "char"), gp = gpar(col = "black", fill = KEY_GLYPH, fontsize = DOT$size * .pt + DOT$stroke * .stroke / 2, lwd = DOT$stroke * .stroke / 2)))
cross_glyph <- function(label) list(label = label, w_sym = .16, draw = function(x, y)
  grid.points(unit(x + .05, "in"), unit(H - y - .04, "in"), pch = 13, size = unit(1, "char"), gp = gpar(col = KEY_GLYPH, fontsize = RING$size * .pt + RING$stroke * .stroke / 2, lwd = RING$stroke * .stroke / 2)))
thr_glyph <- function(label) list(label = label, w_sym = .12, draw = function(x, y)
  grid.lines(unit(x + .045, "in"), unit(H - y - c(-.01, .09), "in"), gp = gpar(col = "grey45", lwd = .3 * .pt, lty = "22")))
text_glyph <- function(glyph, face, label) { w <- meas(glyph, face) + .05
  list(label = label, w_sym = w, draw = function(x, y) text_at(glyph, x, y - .005, face = face)) }
KEY1 <- list(ring_glyph("Group alone"), dot_glyph("Group with its mixed genes"), cross_glyph("No group gene in the term"))
KEY2 <- list(thr_glyph("Adjusted p = 0.05"), text_glyph("Term", "bold", "also in Fig. 2E"), text_glyph("Term", "plain", "only here"))

draw_panel <- function() { key_row(KEY1, X0, PW, yK, gap = .22); key_row(KEY2, X0, PW, yK + .17, gap = .22) }
