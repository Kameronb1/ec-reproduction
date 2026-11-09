# Figure 4C: EC identity of SC-SIC and SC-beta cells at the start and end of each row. Left, two plots: the per-cell
# distribution of the four-marker score and of the SC-SIC 20-gene classifier's score minus its fold's threshold,
# one ridge per population and end (start: outline; end: filled), each sample weighted equally, with the cutoff for
# classifying a cell as EC.
# Right, two plots: expression change (log2) of the four markers and of the 20 SC-SIC genes, one bar per population.
# Rows as in A, with grey bands on alternate rows and the row labels at the left.

CUT4 <- 0.4645893053562607                                     # the four-marker classifier's cut-off on its score (Figure 2A)
C_OFF <- c(SC_SIC = -.04, SC_BETA = .46); C_RIDGE_H <- .42     # ridge baselines within a row; height of the taller ridge
XL1A <- c(-.95, 1.85); XL1B <- c(-1.65, 1.4)                   # the score axes: the full range of the cells' scores

# Per-cell scores of each row's start and end cells. A cell's weight is 1 / its sample's cells / the samples at its end,
# so the weighted share of cells classified as EC is the row end's mean share over samples: the check against the
# shares in rows_called.tsv (graft rows: SoupX-corrected results, rows_called_soupx.tsv), which carries 3 decimals.
SX <- read_table("rows_called_soupx.tsv")
calls <- function(rule) {
  x <- RC[RC$rule == rule, ]; k <- match(paste(rep(ROWS$id, each = 2), c("SC_SIC", "SC_BETA")), paste(x$id, x$population))
  d <- data.frame(id = x$id[k], measure = x$population[k], from = x$from[k], to = x$to[k])
  s <- match(paste(d$id, d$measure, rule), paste(SX$id, SX$population, SX$rule)); g <- !is.na(s)
  d$from[g] <- SX$from[s[g]]; d$to[g] <- SX$to[s[g]]; d }
PCG <- calls("FOUR_MARKER_GATED"); PCL <- calls("CLASSIFIER")
CC <- read_table("C_cells.tsv.gz")
stopifnot(all(CC$four_marker_score[CC$four_marker_gated] >= CUT4))                     # every cell classified as EC clears the cutoff
chk <- do.call(rbind, lapply(split(CC, paste(CC$id, CC$population, CC$end)), function(x) {
  P <- function(d) d[[x$end[1]]][d$id == x$id[1] & d$measure == x$population[1]]
  c(100 * sum(x$weight * x$four_marker_gated) - P(PCG), 100 * sum(x$weight * x$classifier_call) - P(PCL)) }))
stopifnot(nrow(chk) == 4L * nrow(ROWS), max(abs(chk)) < 1e-2)
# Ridges: weighted Gaussian kernel densities (bandwidth bw.nrd0) over the axis range; each pair of ridges (start, end) is
# scaled to its taller one, C_RIDGE_H rows high.
ridges <- function(value, xl) do.call(rbind, lapply(split(CC, paste(CC$id, CC$population)), function(x) {
  yb <- rowy(x$id[1]) + C_OFF[[x$population[1]]]
  dd <- lapply(c(from = "from", to = "to"), function(e) { z <- x[x$end == e, ]
    density(z[[value]], weights = z$weight / sum(z$weight), bw = bw.nrd0(z[[value]]), from = xl[1], to = xl[2], n = 256) })
  dm <- max(dd$from$y, dd$to$y)
  do.call(rbind, lapply(names(dd), function(e) data.frame(population = x$population[1], end = e, g = paste(x$id[1], x$population[1], e),
    x = c(dd[[e]]$x, rev(dd[[e]]$x)), y = c(yb - C_RIDGE_H * dd[[e]]$y / dm, rep(yb, 256))))) }))
RC1A <- ridges("four_marker_score", XL1A); RC1B <- ridges("classifier_margin", XL1B)
# Bars: expression change of the four markers and the 20 SC-SIC genes per row and population; the axis covers every bar.
CB <- read_table("C_bars.tsv"); CB <- CB[CB$id %in% ROWS$id & CB$population %in% c("SC_SIC", "SC_BETA"), ]
CB$y <- rowy(CB$id) + c(SC_SIC = -.18, SC_BETA = .18)[CB$population]; CB$col <- POP_COL[CB$population]
BV <- c(CB$four_marker, CB$sc_sic_20_gene)
XL2 <- c(min(-3.9, floor(min(BV) * 2) / 2 - .2), max(1.2, ceiling(max(BV) * 2) / 2 + .2))

# The four plots without their marks set the layout: the bar plots' panels are .86 in wide, the score plots share the
# rest, and the last panel ends .25 in before D's.
f1a <- colplot(XL1A, c(-.5, 0, .5, 1, 1.5), c("-0.5", "0", "0.5", "1", "1.5"), lab_y = TRUE) + Y_PAD + theme(plot.margin = margin(2, 9, 0, 2))
f1b <- colplot(XL1B, c(-1, 0, 1), c("-1", "0", "1")) + theme(plot.margin = margin(2, 9, 0, 6))
f2 <- colplot(XL2, seq(ceiling(XL2[1] / 2) * 2, floor(XL2[2] / 2) * 2, 2), function(x) ifelse(x == 0, "0", sprintf("%+g", x)), xline = 0,
              minor = seq(ceiling(XL2[1]), floor(XL2[2]), 1)) + theme(plot.margin = margin(2, 4, 0, 4), panel.grid.minor.x = element_line(color = "#E4E4E4", linewidth = .2))
NP <- vapply(list(f1a, f1b, f2, f2), function(p) nonpanel_w(ggplotGrob(p)), 0)
C_RIGHT <- D_X + side_w(ggplotGrob(d_frame()), "left") + side_w(ggplotGrob(f2), "right") - .25
P1 <- (C_RIGHT - .06 - sum(NP) - 2 * .86) / 2
C_W <- c(NP[1] + P1, NP[2] + P1, NP[3] + .86, NP[4] + .86); C_X <- .06 + c(0, cumsum(C_W)[1:3])
C_PLACE <- row_place(align_h(list(f1a, f1b, f2, f2, d_frame()))[[1]], CD_TOP, CD_GROW)

ridge_layers <- function(R, cut) list(geom_vline(xintercept = cut, color = "#555555", linetype = "22", linewidth = .3),
  geom_polygon(data = R[R$end == "to", ], aes(x, y, group = g, fill = population), alpha = .55, color = NA),
  geom_polygon(data = R[R$end == "from", ], aes(x, y, group = g), fill = NA, color = "#333333", linewidth = .22), scale_fill_manual(values = POP_COL))
pC <- list(f1a + ridge_layers(RC1A, CUT4), f1b + ridge_layers(RC1B, 0),
           f2 + geom_segment(data = CB, aes(x = 0, xend = four_marker, y = y, yend = y, color = I(col)), linewidth = .9),
           f2 + geom_segment(data = CB, aes(x = 0, xend = sc_sic_20_gene, y = y, yend = y, color = I(col)), linewidth = .9))
HEADERS <- list(c("Four-marker", "classifier"), c("SC-SIC 20-gene", "classifier"), "Four markers", "20 SC-SIC genes")

draw_panel <- function() shifted_up(ROW_UP[["cd"]], {
  g <- align_h(lapply(pC, nobg))
  bx <- lapply(1:4, function(i) panels_in(g[[i]], C_X[i], C_PLACE[["y"]], C_W[i], C_PLACE[["h"]])[[1]])
  band_rects(bx, label_x = .06)
  for (i in 1:4) draw_at(g[[i]], C_X[i], C_PLACE[["y"]], C_W[i], C_PLACE[["h"]])
  # x-axis titles under the plots (the bars' one spans both bar plots); column headers over them
  y_ax <- C_PLACE[["y"]] + C_PLACE[["h"]] + CD_XTITLE_GAP
  cd_xtitle("Four-marker score", bx[[1]], bx[[1]], y_ax); cd_xtitle("Classifier score − cutoff", bx[[2]], bx[[2]], y_ax)
  cd_xtitle("Expression change (log2)", bx[[3]], bx[[4]], y_ax)
  for (i in 1:4) col_header(HEADERS[[i]], bx[[i]])
  section_labels(g[[1]], bx[[1]], CD_SUB_BASE)
})
