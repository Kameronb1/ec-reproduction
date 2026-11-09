# Figure 4E-F key, under both panels. Line 1: the arrow marks E and F share (open circle at the start, arrow to the graft
# or continued culture) and the population colors; line 2: E's dashed native intestinal EC line and its pale middle-half bars,
# one per population, stacked as in E. F's ticks get two lines of their own, starting at F's panel edge (ef_layout() in helpers.R).

EF <- ef_layout()
K_START_X <- .10 + CIRCLE_R                                    # the open circle's ink starts at .10 in
K_POP_X <- K_START_X + (.12 + str_in("Start") + .2 + .23 + str_in("Graft / continued culture")) + .30
E_BAR_NOTE <- "middle half of genes (E)"
F_KEY <- c(lead = "Intestinal studies:", fetal = "fetal", adult = "adult", outside = "not used to choose the genes")

# open circle 'Start', arrow 'Graft / continued culture'
start_graft_key <- function(x, y) {
  grid.points(x = unit(x, "in"), y = unit(H - y - .045, "in"), pch = 21, size = unit(CIRCLE_PT, "pt"), gp = gpar(col = "#555555", fill = "white", lwd = CIRCLE_LWD))
  text_at("Start", x + .12, y - .005)
  x2 <- x + .12 + str_in("Start") + .2
  grid.lines(x = unit(c(x2, x2 + .18), "in"), y = unit(rep(H - y - .045, 2), "in"), arrow = ARW, gp = gpar(col = "#555555", fill = "#555555", lwd = 1.3))
  text_at("Graft / continued culture", x2 + .23, y - .005) }
# a tick drawn as on F's axis (up from its foot), then its name; returns the next x
tick <- function(b, x, y) {
  col <- if (b == "outside") "#555555" else TICK_COL[[b]]; len <- TICK_LEN[[if (b == "outside") "outside" else "chose"]]
  grid.lines(unit(rep(x + .03, 2), "in"), unit(H - (y + .085) + c(0, len), "in"), gp = gpar(col = col, lwd = .6 * ggplot2::.pt, lineend = "butt"))
  text_at(F_KEY[[b]], x + .10, y - .005, 7); x + .10 + str_in(F_KEY[[b]], 7) + .16 }

draw_panel <- function() shifted_up(ROW_UP[["ef"]], {
  # the two lines, left-aligned with each other, centered as a block under E (page margin, where E's row labels start, to
  # its last facet's right edge); the block must end .30 in before F's tick key
  end1 <- K_POP_X + sum(vapply(POP_LAB[c("SC_SIC", "SC_BETA", "SC_ALPHA")], function(l) str_in(l) + .26 + .16, 0)) - .16
  end2 <- .10 + .26 + str_in("native intestinal EC (E)") + .20 + .26 + str_in(E_BAR_NOTE)
  eR <- max(vapply(EF$bxE, function(b) b$left_in + b$width_in, 0)); xf <- EF$bxF[[1]]$left_in
  dx <- (.10 + eR) / 2 - (.10 + max(end1, end2)) / 2
  stopifnot(max(end1, end2) + dx <= xf - .30)
  start_graft_key(K_START_X + dx, EF_KEY1); population_key(K_POP_X + dx, EF_KEY1)
  key_line(.10 + dx, EF_KEY2, "#555555", lwd = .4 * ggplot2::.pt, lty = "dashed", butt = TRUE); text_at("native intestinal EC (E)", .10 + dx + .26, EF_KEY2 - .005)
  xb <- .10 + dx + .26 + str_in("native intestinal EC (E)") + .20
  for (k in names(POP_COL)) grid.lines(unit(c(xb, xb + .2), "in"), unit(rep(H - EF_KEY2 - .045 - OFF3[[k]] * .13, 2), "in"),
                                       gp = gpar(col = POP_COL[[k]], alpha = .3, lwd = 1 * ggplot2::.pt, lineend = "butt"))
  text_at(E_BAR_NOTE, xb + .26, EF_KEY2 - .005)
  text_at(F_KEY[["lead"]], xf, EF_KEY1 - .005, 7)
  x <- xf + str_in(F_KEY[["lead"]], 7) + .10; for (b in c("fetal", "adult")) x <- tick(b, x, EF_KEY1)
  tick("outside", xf - .03 + .3 * ggplot2::.pt / 96, EF_KEY2)          # its ink starting at F's panel edge
})
