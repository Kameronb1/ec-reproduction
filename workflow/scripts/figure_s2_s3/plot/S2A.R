# Figure S2A: four-marker score histogram of the 1,164 fitting cells, colored by markers detected, with the
# two-component mixture and the cutoff.
FIT <- jsonlite::fromJSON(result_path("fitted_constants.json"))
CUT <- FIT$ec_core_cutoff
stopifnot(FIT$n_cells == 1164L)
hist_A <- read_table("A_histogram.tsv"); hist_A$markers <- factor(hist_A$markers, levels = 4:0)
curv_A <- read_table("A_mixture_curves.tsv")
stopifnot(sum(hist_A$cells) == 1164L)
CURVE_COL <- c(Mixture = "#555555", Lower = "#447AA5", Upper = "#CC7636")
FIT_STUDIES <- c("Elmentaite 2020", "Burclaff 2022", "Wang 2020", "Yu 2021a")
ymaxA <- max(hist_A$cells, curv_A$y)
plot_A <- ggplot(hist_A, aes(x, cells, fill = markers)) +
  geom_col(width = .094) +
  geom_line(data = curv_A, aes(x, y, colour = curve, linetype = curve, group = curve), inherit.aes = FALSE, linewidth = .4) +
  geom_vline(xintercept = CUT, linetype = "22", linewidth = .4) +
  annotate("text", x = CUT + .07, y = ymaxA * 1.02, label = sprintf("Cutoff %.2f", CUT), hjust = 0, vjust = 1, size = MM, family = FONT) +
  annotate("text", x = .62, y = ymaxA * .86, label = "Fitting studies", hjust = 0, vjust = 1, size = MM, family = FONT, fontface = "bold") +
  annotate("text", x = .62, y = ymaxA * .78, label = paste(FIT_STUDIES, collapse = "\n"), hjust = 0, vjust = 1, size = MM,
           family = FONT, lineheight = .95) +
  scale_fill_manual(values = MARK_COL) + scale_colour_manual(values = CURVE_COL) +
  scale_linetype_manual(values = c(Mixture = "solid", Lower = "22", Upper = "22")) +
  scale_x_continuous(breaks = -1:2, expand = expansion(mult = c(.03, .04))) +
  scale_y_continuous(expand = expansion(mult = c(0, .04))) +
  labs(x = "Four-marker score", y = "Fitting cells (1,164)") + theme_p2() +
  theme(panel.grid.major.x = element_blank(), axis.line.y = element_line(colour = "black", linewidth = .3),
        axis.ticks.y = element_line(colour = "black", linewidth = .3))

draw_A <- function(x0, y0, w, h_plot) {
  draw_at(plot_A, x0, y0, w, h_plot)
  y <- y0 + h_plot + .16
  txt("Markers detected, of four", x0 + .05, y, just = c("left", "centre"))
  grid.text("TPH1, SLC18A1, FEV, LMX1A", x = unit(x0 + .05, "in"), y = unit(H - y - LINE, "in"), just = c("left", "centre"),
            gp = gp_t(face = "italic"))
  key_row(lapply(as.character(0:4), function(k) k_square(k, MARK_COL[[k]])), x0 + .05, y + 2.1 * LINE, gap = .12)
  yy <- y + 3.4 * LINE
  for (i in 1:3) {
    key_row(list(k_line(c("Two-component Gaussian mixture", "Lower component", "Upper component")[i], CURVE_COL[i],
                        lty = c(1, 2, 2)[i])), x0 + .05, yy + (i - 1) * 1.1 * LINE)
  }
  yy + 2.2 * LINE + .06                                     # bottom of the panel
}

draw_panel <- function() {
  letter("A", L - .1, .12)
  stopifnot(draw_A(L, .32, 2.25, 2.55) == S2_A_BOTTOM)
}
