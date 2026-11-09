# Figure 2E: every EC and other terminal EEC cell (primary study labels), its four-marker score (x) against its
# 20-gene score (y), one plot per study, with both cutoffs. Corner numbers count cells called by both classifiers (top
# right), by the 20-gene classifier only (top left), by the four-marker classifier only (bottom right) and by neither
# (bottom left): EC cells in the study's color, other terminal EEC in grey (no grey count for 'neither'). A cell past
# a cutoff but not classified as EC by that classifier (it failed the classifier's neighborhood check) is drawn as an
# open circle, so the corner counts can be read against the dots. As a check, the cell table must reproduce exactly
# the per-study counts of cells classified as EC (E_concordance.tsv and primary_only_training_calls.tsv).

PLOT_Y <- 7.872; PLOT_H <- 1.64; KEY_Y <- 9.627    # top and height of the plots, top of the key

calls_E <- read_table("E_concordance.tsv")
calls_E <- calls_E[calls_E$seg != "extra", c("family", "fam", "class", "rule", "called", "n")]
# 20-gene classifier: genes from D, weights and threshold fitted on the primary study labels of the other studies
po <- read_table("primary_only_training_calls.tsv")
po <- po[po$variant == "D_GENES_PRIMARY_ONLY", ]
k <- calls_E$rule == "panel"
calls_E$called[k] <- po$called[match(paste(calls_E$family[k], calls_E$class[k]), paste(po$study, po$class))]
# McCauley 2023 HIO grafts
hio_E <- read_table("hio_E_hio_summary.tsv")
hio_E <- hio_E[hio_E$panel_size == 20, ]
calls_E <- rbind(calls_E, data.frame(family = "SF_GSE214852_HIO_GRAFT", fam = "hio",
                                     class = c("POSITIVE", "POSITIVE", "NEGATIVE", "NEGATIVE"), rule = c("marker", "panel"),
                                     called = c(hio_E$ec_rule_called, hio_E$ec_called, hio_E$other_rule_called, hio_E$other_called),
                                     n = c(hio_E$ec_n, hio_E$ec_n, hio_E$other_n, hio_E$other_n)))

cells_E <- read_table("E_cells.tsv")
cells_E$call4 <- as.logical(cells_E$call4); cells_E$call20 <- as.logical(cells_E$call20)
for (i in seq_len(nrow(calls_E))) {
  d <- cells_E[cells_E$fam == calls_E$fam[i] & cells_E$class == c(POSITIVE = "EXPLICIT_EC", NEGATIVE = "TERMINAL_NON_EC_EEC")[[calls_E$class[i]]], ]
  stopifnot(nrow(d) == calls_E$n[i], sum(if (calls_E$rule[i] == "marker") d$call4 else d$call20) == calls_E$called[i])
}
E_STUDY_COL <- c(gut_cell_atlas = "#43A5B6", hickey = "#43A5B6", egozi = "#43A5B6", gse238274 = "#CCBB44", hio = "#6943AC")
E_OTHER <- "#9A9A9A"
E_XR <- c(-1.5, 2.6); E_YR <- range(cells_E$score20)
quad <- function(d) c(both = sum(d$call4 & d$call20), only20 = sum(!d$call4 & d$call20), only4 = sum(d$call4 & !d$call20),
                      none = sum(!d$call4 & !d$call20))
scatter_E <- function(f, first) {
  d <- cells_E[cells_E$fam == f, ]; d <- d[order(d$class == "EXPLICIT_EC"), ]
  d$col <- ifelse(d$class == "EXPLICIT_EC", E_STUDY_COL[[f]], E_OTHER)
  d$open <- (d$score4 >= EC_CUTOFF & !d$call4) | (d$score20 >= d$thr20 & !d$call20)
  thr <- unique(d$thr20); stopifnot(length(thr) == 1)
  qe <- quad(d[d$class == "EXPLICIT_EC", ]); qo <- quad(d[d$class != "EXPLICIT_EC", ])
  yr <- E_YR + c(-.22, .22) * diff(E_YR); dy <- .075 * diff(yr); ins <- .03 * diff(yr)
  x <- c(E_XR[2] - .1, E_XR[1] + .1, E_XR[2] - .1, E_XR[1] + .1); h <- c(1, 0, 1, 0)
  ann <- rbind(data.frame(x = x, y = c(yr[2] - ins, yr[2] - ins, yr[1] + ins, yr[1] + ins), h = h, v = c(1, 1, 0, 0), lab = qe, col = E_STUDY_COL[[f]]),
               data.frame(x = x, y = c(yr[2] - ins - dy, yr[2] - ins - dy, yr[1] + ins + dy, yr[1] + ins + dy), h = h, v = c(1, 1, 0, 0), lab = qo, col = E_OTHER))
  pc <- function(q) { p <- 100 * q / sum(q); ifelse(q == 0, "0%", ifelse(p < 1, "<1%", ifelse(p > 99 & q < sum(q), ">99%", paste0(round(p), "%")))) }
  ann$lab <- c(pc(qe), pc(qo))
  ggplot(d, aes(score4, score20)) +
    geom_point(data = d[!d$open, ], aes(color = col), size = .35, stroke = 0, alpha = .75) +
    geom_point(data = d[d$open, ], aes(color = col), shape = 1, size = .55, stroke = .2) +
    geom_vline(xintercept = EC_CUTOFF, linetype = "dashed", linewidth = .3) +
    geom_hline(yintercept = thr, linetype = "dashed", linewidth = .3) +
    geom_text(data = ann, aes(x, y, label = lab, hjust = h, vjust = v, color = col), size = 6.5 / ggplot2::.pt,
              family = "Arial", fontface = "bold", inherit.aes = FALSE) +
    scale_color_identity() +
    scale_x_continuous(limits = E_XR, breaks = -1:2, expand = c(0, 0)) +
    scale_y_continuous(limits = yr, expand = c(0, 0)) +
    labs(x = "Four-marker score", y = if (first) "20-gene score" else NULL, title = STUDY[f], subtitle = CONTEXT[f]) +
    theme_fig2() +
    theme(panel.border = element_rect(fill = NA, color = "black", linewidth = .3), axis.line.x = element_blank(),
          axis.ticks.y = element_line(color = "black", linewidth = .3), panel.grid.major = element_blank(),
          plot.title = element_text(size = 7.5, face = "bold", margin = margin(0, 0, 0, 0)),
          plot.subtitle = element_text(size = 7, margin = margin(0, 0, 1, 0)),
          axis.title.x = element_text(margin = margin(1, 0, 0, 0))) +
    (if (first) NULL else theme(axis.text.y = element_blank()))
}
plot_E <- wrap_plots(lapply(seq_along(FAMS), function(i) scatter_E(FAMS[i], i == 1)), nrow = 1)

draw_panel <- function() {
  draw_at(plot_E, .08, PLOT_Y, 7.32, PLOT_H)
  y <- KEY_Y; left <- .10; width <- 7.30
  labs <- c("EC (primary study label)", "Other terminal EEC", "Cutoffs", "Past a cutoff, not called (neighborhood check)",
            "Numbers: percent of each label's cells called by both, by one or by neither classifier")
  ww <- vapply(labs, str_in, 0) + c(.36, .14, .30, .14, 0); gap <- .16
  x <- left + (width - sum(ww[1:4]) - 3 * gap) / 2
  for (j in 1:3) grid.points(unit(x + .045 + (j - 1) * .11, "in"), unit(H - y - .05, "in"), pch = 16,
                             gp = gpar(col = c("#43A5B6", "#CCBB44", "#6943AC")[j], fontsize = 5))
  text_at(labs[1], x + .36, y); x <- x + ww[1] + gap
  grid.points(unit(x + .045, "in"), unit(H - y - .05, "in"), pch = 16, gp = gpar(col = E_OTHER, fontsize = 5))
  text_at(labs[2], x + .14, y); x <- x + ww[2] + gap
  grid.lines(unit(c(x, x + .22), "in"), unit(rep(H - y - .05, 2), "in"), gp = gpar(lty = 2, lwd = .3 * ggplot2::.pt))
  text_at(labs[3], x + .3, y); x <- x + ww[3] + gap
  grid.points(unit(x + .045, "in"), unit(H - y - .05, "in"), pch = 1, gp = gpar(col = E_OTHER, fontsize = 5, lwd = .5))
  text_at(labs[4], x + .14, y)
  text_at(labs[5], left + (width - str_in(labs[5])) / 2, y + .16)       # second key line
}
