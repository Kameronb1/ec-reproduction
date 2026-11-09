# Figure 2A: marker-score histogram of the 1,199 fitting cells, colored by how many of the four EC markers each
# cell has a read of, with the two-component Gaussian mixture, its components and the classifier's cutoff.
# Key below the plot.

PLOT_Y <- .139; PLOT_H <- 1.965; KEY_Y <- 2.174    # top and height of the plot, top of the key

hist_A <- read_table("A_histogram.tsv")
hist_A$markers <- factor(hist_A$markers, levels = 4:0)
curves_A <- read_table("A_mixture_curves.tsv")
fitting_studies <- sort(SN[c("SF_GSE185224_ADULT_REGIONAL_INTESTINE", "SF_ELMENTAITE2020_DEVGUT",
                             "SF_GSE125970_ADULT_INTESTINE", "SF_YU2021_ENDODERM:fetal"), "author_year"])
MARKER_COLS <- c("0" = "#ECE8DF", "1" = "#D8C5AA", "2" = "#CB9B72", "3" = "#B86E46", "4" = "#8E3D29")
CURVE_COLS <- c(Mixture = "#555555", Lower = "#447AA5", Upper = "#CC7636")

plot_A <- ggplot(hist_A, aes(x, cells, fill = markers)) +
  geom_col(width = .094) +
  geom_line(data = curves_A[curves_A$curve == "Mixture", ], aes(x, y, color = curve, linetype = curve),
            inherit.aes = FALSE, linewidth = .4) +
  geom_line(data = curves_A[curves_A$curve != "Mixture", ], aes(x, y, color = curve, linetype = curve, group = curve),
            inherit.aes = FALSE, linewidth = .4) +
  geom_vline(xintercept = EC_CUTOFF, linetype = "dashed", linewidth = .4, color = "black") +
  annotate("text", x = EC_CUTOFF - .07, y = 262, label = sprintf("Cutoff\n%.4f", EC_CUTOFF), hjust = 1, vjust = 1,
           size = TM, color = "black", family = "Arial") +
  annotate("text", x = .60, y = 262, label = "Fitting studies", hjust = 0, vjust = 1, size = TM, fontface = "bold",
           color = "black", family = "Arial") +
  annotate("text", x = .60, y = 240, label = paste(fitting_studies, collapse = "\n"), hjust = 0, vjust = 1,
           lineheight = .95, size = TM, color = "black", family = "Arial") +
  scale_fill_manual(values = MARKER_COLS) +
  scale_color_manual(values = CURVE_COLS) +
  scale_linetype_manual(values = c(Mixture = "solid", Lower = "dashed", Upper = "dashed")) +
  scale_x_continuous(breaks = c(-1, 0, 1, 2), expand = expansion(mult = c(.03, .04))) +
  scale_y_continuous(breaks = c(0, 100, 200), expand = expansion(mult = c(0, .07))) +
  labs(x = "Four-marker score", y = "Fitting cells") +
  theme_fig2() + theme(panel.grid.major.x = element_blank(), axis.line.y = element_line(color = "black", linewidth = .3))

draw_panel <- function() {
  draw_at(plot_A, .10, PLOT_Y, 2.20, PLOT_H)
  # key: marker-count colors, then the three curves
  y <- KEY_Y
  text_at("Markers detected (of 4)", .18, y)
  grid.text("TPH1, SLC18A1, FEV, LMX1A", x = unit(.18, "in"), y = unit(H - (y + .115), "in"), just = c("left", "top"),
            gp = gpar(fontfamily = "Arial", fontsize = PT, fontface = "italic", col = "black"))
  for (i in 0:4) key(as.character(i), MARKER_COLS[as.character(i)], .20 + i * .40, y + .26)
  curve_labels <- c("Gaussian mixture fit (2 components)", "Lower Gaussian component", "Upper Gaussian component")
  for (i in 1:3) {
    yy <- y + .42 + (i - 1) * .125
    grid.lines(x = unit(c(.20, .40), "in"), y = unit(c(H - yy - .04, H - yy - .04), "in"),
               gp = gpar(col = CURVE_COLS[i], lty = c(1, 2, 2)[i], lwd = .7))
    text_at(curve_labels[i], .46, yy)
  }
}
