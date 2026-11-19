# Figure S3A: percent of EC-labeled cells and of other terminal EEC classified as EC across cutoffs in each test
# study, with the five leave-one-fitting-group-out refits as diamonds and the range of their positions shaded.
CUT <- jsonlite::fromJSON(result_path("fitted_constants.json"))$ec_core_cutoff
d <- read_table("S3A_classified_across_cutoffs.tsv")
FAM_E <- c("gut_cell_atlas", "egozi", "hickey", "gse214852_hio")
SRC_E <- c(gut_cell_atlas = EC_GUT, egozi = EC_GUT, hickey = EC_GUT, gse214852_hio = EC_HIO)
HEAD_E <- list(TITLE_B$gut_cell_atlas, TITLE_B$egozi, TITLE_B$hickey, TITLE_B$hio)
s2 <- d[d$kind == "curve", c("family", "class", "cutoff", "post_rate")]
s2$st <- factor(s2$family, levels = FAM_E); s2$col <- ifelse(s2$class == "Primary-labeled EC", SRC_E[s2$family], OTHER_EEC)
s2$grp <- paste(s2$family, s2$class)
r2 <- d[d$kind == "refit", c("family", "class", "cutoff", "post_rate")]; r2$st <- factor(r2$family, levels = FAM_E)
stopifnot(nrow(r2) == 20L, all(table(r2$family) == 5L), setequal(unique(d$family), FAM_E))
BAND <- range(r2$cutoff)                                   # the five refit positions
plot_E <- ggplot(s2[order(s2$grp, s2$cutoff), ], aes(cutoff, post_rate)) +
  annotate("rect", xmin = BAND[1], xmax = BAND[2], ymin = -Inf, ymax = Inf, fill = "#E9E9E9") +
  geom_vline(xintercept = CUT, linetype = "22", linewidth = .35) +
  geom_step(aes(colour = col, group = grp), linewidth = .5, direction = "hv") +
  geom_point(data = r2, aes(cutoff, post_rate), shape = 23, fill = "white", colour = "black", size = 1.3, stroke = .35) +
  scale_colour_identity() + facet_grid(. ~ st) +
  scale_x_continuous(breaks = c(-.5, 0, .5, 1, 1.5), labels = c("-0.5", "0", "0.5", "1", "1.5")) +
  scale_y_continuous(breaks = c(0, 50, 100), labels = pct_lab) +
  coord_cartesian(xlim = c(-.9, 1.5), ylim = c(-2.5, 103), expand = FALSE) +
  labs(x = "Four-marker score cutoff", y = "Classified as EC\n(% of cells with that label)") + theme_p2() +
  theme(panel.spacing.x = unit(.2, "in"), axis.line.y = element_line(colour = "black", linewidth = .3),
        axis.ticks.y = element_line(colour = "black", linewidth = .3), axis.title.y = element_text(lineheight = .95, margin = margin(r = 3)))

draw_E <- function(x0, y0, w, h) {
  facet_heads(plot_E, x0, y0, w, h, HEAD_E, 4, .2)
  draw_at(plot_E, x0, y0 + 2 * LINE + .05, w, h)
  yk <- y0 + 2 * LINE + .05 + h + .13
  key_row(list(k_vdash(sprintf("Cutoff used, %.2f", CUT)), k_square("Range of the five refitted cutoffs", "#E9E9E9"),
               k_point("Refit leaving out one fitting group", "white", pch = 23)), x0, yk, width = w)
  key_row(list(k_line("EC, primary intestine", EC_GUT), k_line("EC, HIO graft", EC_HIO),
               k_line("Other terminal EEC", OTHER_EEC)), x0, yk + 1.5 * LINE, width = w)
  yk + 1.5 * LINE + .07
}

draw_panel <- function() {
  letter("A", L - .1, .12)
  stopifnot(draw_E(L + .05, .16, R - L - .05, 1.45) == S3_A_BOTTOM)
}
