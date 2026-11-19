# Figure S3B: sensitivity and specificity of the four-marker classifier after thinning each test study's counts to 25%
# or 50% (three draws each; 100% = full data).
dp <- rbind(read_table("four_marker_depth.tsv"), read_table("hio_four_marker_depth.tsv"))
ST_F <- c("SF_ELMENTAITE2021_GUT_ATLAS", "SF_EGOZI2021_DEVELOPMENTAL_EEC", "SF_HICKEY2023_ADULT_EEC", "SF_GSE214852_HIO_GRAFT")
stopifnot(setequal(unique(dp$study), ST_F), all(table(dp$study) == 7L))
long_F <- do.call(rbind, lapply(c("sensitivity", "specificity"), function(m)
  data.frame(study = dp$study, fraction = dp$fraction, replicate = dp$replicate, metric = m, value = dp[[m]])))
long_F$st <- factor(long_F$study, levels = ST_F); long_F$pct <- 100 * long_F$fraction
long_F$metric <- factor(long_F$metric, levels = c("sensitivity", "specificity"),
                        labels = c("EC labels classified", "Other terminal EEC not classified"))
mean_F <- aggregate(value ~ st + metric + pct, long_F, mean)
# one plot per measure, stacked, each with its own row title; the two share the facet columns and the x axis
mk_F <- function(m, x_axis) {
  d <- long_F[long_F$metric == m, ]; dm <- mean_F[mean_F$metric == m, ]
  ggplot(d, aes(pct, value)) +
    geom_line(data = dm, colour = CLS_COL[["HIGH_CONFIDENCE_EC"]], linewidth = .5) +
    geom_point(shape = 21, fill = CLS_COL[["HIGH_CONFIDENCE_EC"]], colour = "black", stroke = EDGE, size = 1.1) +
    facet_grid(. ~ st) +
    scale_x_continuous(breaks = c(25, 50, 100), labels = pct_lab, limits = c(15, 110)) +
    scale_y_continuous(limits = c(0, 1), breaks = c(0, .5, 1), labels = c("0", "0.5", "1")) +
    labs(x = if (x_axis) "Counts kept in the scored study" else NULL, y = "Proportion") + theme_p2() +
    theme(panel.spacing.x = unit(.2, "in"), axis.line.y = element_line(colour = "black", linewidth = .3),
          axis.ticks.y = element_line(colour = "black", linewidth = .3),
          axis.text.x = if (x_axis) element_text(colour = "black", size = PT) else element_blank(),
          axis.title.y = element_text(margin = margin(r = 3)))
}
plot_F1 <- mk_F("EC labels classified", FALSE); plot_F2 <- mk_F("Other terminal EEC not classified", TRUE)
HEAD_F <- list(c("Elmentaite 2021", ""), c("Egozi 2021", ""), c("Hickey 2023 (nuclei)", ""), c("McCauley 2023", ""))
draw_F <- function(x0, y0, w, h) {               # h: height of each of the two plot rows
  facet_heads(plot_F1, x0, y0, w, h, HEAD_F, 4, .2)
  o <- panel_offsets(plot_F1, w, h)
  y1 <- y0 + LINE + .04
  txt("EC labels classified as EC (sensitivity)", x0 + o$l, y1, face = "bold")
  draw_at(plot_F1, x0, y1 + LINE + .02, w, h)
  y2 <- y1 + LINE + .02 + h + .1
  txt("Other terminal EEC not classified as EC (specificity)", x0 + o$l, y2, face = "bold")
  o2 <- panel_offsets(plot_F2, w, h); h2 <- h + ((h - (o2$b - o2$t)) - (h - (o$b - o$t)))   # equal panel heights in both rows
  draw_at(plot_F2, x0, y2 + LINE + .02, w, h2)
  o3 <- panel_offsets(plot_F2, w, h2); stopifnot(abs(o3$l - o$l) < 1e-6, abs((o3$b - o3$t) - (o$b - o$t)) < 1e-6)   # same columns, same panel height
  yk <- y2 + LINE + .02 + h2 + .13
  key_row(list(k_point("One draw", CLS_COL[["HIGH_CONFIDENCE_EC"]]), k_line("Mean of three draws (100%: full data)", CLS_COL[["HIGH_CONFIDENCE_EC"]])),
          x0, yk, width = w)
  yk + .07
}

draw_panel <- function() {
  letter("B", L - .1, S3_B_TOP)
  stopifnot(draw_F(L + .05, S3_B_TOP + .02, R - L - .05, .8) == S3_B_BOTTOM)
}
