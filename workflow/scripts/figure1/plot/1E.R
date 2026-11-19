# Figure 1E: four-marker score of SC-islet cells, one row per study. In each row the SC-SIC, SC-beta and other
# endocrine cells (mostly SC-alpha and -delta) each get a density scaled to its own peak; the dashed line is the cutoff
# of the four-marker EC classifier fitted on intestinal cells; at the right, the share of SC-SIC and of SC-beta cells
# above it. The scores and the cutoff come from the classifier unit (Supplementary Figures S2 and S3).
suppressPackageStartupMessages(library(ggridges))

E_bc <- read_table("sc_islet_scores_by_cell.tsv.gz", figure = "figure_s2_s3")
E_cutoff <- read_table("sc_islet_cutoff.tsv", figure = "figure_s2_s3")$cutoff
E_top <- read_table("sc_islet_above_cutoff.tsv")                                             # rows top to bottom
E_COLS <- c("SC-SIC" = "#AA3377", "SC-beta" = "#4477AA", "Other endocrine" = "#8C8C8C")
E_bc$who <- factor(ifelse(E_bc$class == "SC_SIC", "SC-SIC", ifelse(E_bc$class == "SC_BETA", "SC-beta", "Other endocrine")), rev(names(E_COLS)))
stopifnot(setequal(unique(E_bc$study), E_top$study))
E_bc$study <- factor(E_bc$study, rev(E_top$study))                                             # the first level is the bottom row
E_pct <- data.frame(study = factor(levels(E_bc$study), levels(E_bc$study)), sic = rev(E_top$sc_sic_pct), beta = rev(E_top$sc_beta_pct))
E_SCALE <- 1.2                                                                               # ridge height in rows; ridges may overlap the row above
E_TOP <- length(levels(E_bc$study)) + E_SCALE + 0.05                                         # room for the top row's curves
E_XLIM <- c(-1, ceiling(max(E_bc$score) * 20) / 20)                                          # the axis ends where the data end
E_RCOL <- 1.15                                                                              # in: room right of the plot for the two columns
E_PLOT_H <- 1.50
theme_fig1e <- function() theme_minimal(base_family = "Arial", base_size = PT) + theme(
  text = element_text(family = "Arial", colour = "black", size = PT),
  axis.text = element_text(colour = "black", size = PT), axis.title = element_text(colour = "black", size = PT),
  strip.text = element_text(colour = "black", size = PT, face = "bold", hjust = 0, margin = margin(1, 0, 2, 0)),
  panel.grid.minor = element_blank(), panel.grid.major = element_line(colour = "#E4E4E4", linewidth = 0.2),
  axis.line.x = element_line(colour = "black", linewidth = 0.3), axis.ticks.x = element_line(colour = "black", linewidth = 0.3),
  axis.ticks.length = unit(1.5, "pt"), legend.position = "none", plot.margin = margin(1, 3, 1, 1))
pE <- ggplot(E_bc, aes(score, study, fill = who, colour = who, group = interaction(study, who))) +
  geom_vline(xintercept = E_cutoff, linetype = "22", colour = "grey20", linewidth = 0.4) +
  geom_density_ridges(aes(height = after_stat(ndensity)), stat = "density", trim = TRUE, scale = E_SCALE, alpha = 0.25, linewidth = 0.3, bw = 0.06) +
  scale_fill_manual(values = E_COLS, guide = "none") + scale_colour_manual(values = E_COLS, guide = "none") +
  scale_x_continuous(limits = E_XLIM, breaks = seq(-1, 1.5, .5), labels = c("-1", "-0.5", "0", "0.5", "1", "1.5"), expand = c(0, 0)) +
  scale_y_discrete(expand = expansion(add = c(0.1, 0))) + coord_cartesian(ylim = c(1, E_TOP), clip = "off") +
  labs(x = "Four-marker score", y = NULL) +
  theme_fig1e() + theme(panel.grid.major.y = element_line(colour = "#E4E4E4", linewidth = 0.2), panel.grid.major.x = element_blank(),
                        plot.margin = margin(14, E_RCOL * 72, 1, 1))                      # top: room for the column headers
E_KEY <- list(square_item("SC-SIC", E_COLS[["SC-SIC"]]), square_item("SC-beta", E_COLS[["SC-beta"]]),
              square_item("Other endocrine cells (mostly SC-alpha, -delta)", E_COLS[["Other endocrine"]]),
              dash_item(sprintf("Cutoff from intestinal cells, %.2f", E_cutoff)))

# x0, y0 = panel top-left (in); w = width; plot_h = plot height; PH = row height. The percentage columns are drawn level
# with each study name; the key row follows below the plot, centred on the panel.
draw_E <- function(x0, y0, w, plot_h, PH) {
  pushViewport(viewport(x = unit(x0, "in"), y = unit(PH - y0, "in"), width = unit(w, "in"), height = unit(plot_h, "in"), just = c("left", "top")))
  g <- ggplotGrob(pE); grid.draw(g)
  # the panel's box inside this viewport (in), from the gtable
  pl <- g$layout[g$layout$name == "panel", ]
  wd <- vapply(seq_along(g$widths), function(i) convertWidth(g$widths[i], "in", valueOnly = TRUE), 0)
  ht <- vapply(seq_along(g$heights), function(i) convertHeight(g$heights[i], "in", valueOnly = TRUE), 0)
  p_left <- sum(wd[seq_len(pl$l - 1)]); p_w <- w - sum(wd[-pl$l])
  p_top <- sum(ht[seq_len(pl$t - 1)]); p_h <- plot_h - sum(ht[-pl$t])
  popViewport()
  yr <- ggplot_build(pE)$layout$panel_params[[1]]$y.range
  y_of <- function(v) y0 + p_top + (yr[2] - v) / diff(yr) * p_h                 # row inches from the top for data y
  x_sic <- x0 + p_left + p_w + .40; x_beta <- x0 + p_left + p_w + .85           # column centres (in)
  yy <- y_of(seq_len(nrow(E_pct)))
  for (i in seq_len(nrow(E_pct))) {
    tx(paste0(E_pct$sic[i], "%"), x_sic, yy[i], PH, PT, col = E_COLS[["SC-SIC"]], just = c("centre", "centre"))
    tx(paste0(E_pct$beta[i], "%"), x_beta, yy[i], PH, PT, col = E_COLS[["SC-beta"]], just = c("centre", "centre")) }
  y_cols <- min(yy) - .17                                                        # column headers (centre)
  tx("SC-SIC", x_sic, y_cols, PH, 7, TRUE, col = E_COLS[["SC-SIC"]], just = c("centre", "centre"))
  tx("SC-beta", x_beta, y_cols, PH, 7, TRUE, col = E_COLS[["SC-beta"]], just = c("centre", "centre"))
  y_span <- y_cols - .145                                                        # the spanning header and its rule
  tx("% above cutoff", (x_sic + x_beta) / 2, y_span, PH, 7, just = c("centre", "centre"))
  grid.segments(unit(x_sic - .2, "in"), unit(PH - (y_span + .075), "in"), unit(x_beta + .2, "in"), unit(PH - (y_span + .075), "in"),
                gp = gpar(col = "grey45", lwd = .5))
  if (y_span - .06 < .02) stop("panel E: the percentage headers run above the row top")
  if (x_beta + .25 > x0 + w) stop("panel E: the percentage columns run past the panel")
  key_row(E_KEY, x0, w, y0 + plot_h + .10, PH)
  panel_heading("E", .15, .02, PH)
}

draw_panel <- function() draw_after_72dpi(function() in_row("E", function(PH) draw_E(.15, .02, W - .30, E_PLOT_H, PH)))
