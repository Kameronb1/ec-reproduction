# Figure 2B: four-marker classifier call rate among the EC cells (primary study label) of each donor, day or HIO
# sample group, in five study plots; a dashed line marks each study's pooled rate. Key below.

PLOT_Y <- .064; PLOT_H <- 2.64; KEY_Y <- 2.824     # top and height of the plot, top of the key

units_B <- read_table("B_units.tsv")[c("family", "unit", "n", "called", "rate", "label")]
# McCauley 2023 HIO grafts: cells classified as EC among the Endocrine_1 (EC) cells of each HIO sample group
hio_B <- read_table("hio_BC_classifier_calls_by_label.tsv")
hio_B <- hio_B[hio_B$author_label == "Endocrine_1", ]
hio_B <- hio_B[order(hio_B$sample), ]
units_B <- rbind(units_B, data.frame(family = "hio", unit = hio_B$sample, n = hio_B$n, called = hio_B$rule_called,
                                     rate = 100 * hio_B$rule_called / hio_B$n,
                                     label = paste0(hio_B$rule_called, "/", hio_B$n)))
POINT_MAX <- 2.9                                    # largest point (mm), for the largest row
B_SHAPE <- c(gse238274 = 17, hio = 15)              # triangles for Singh 2024 days, squares for HIO samples, else circles
B_COL <- c(gse238274 = "#CCBB44", hio = "#6943AC")  # else teal

study_plot_B <- function(f, show_x) {
  d <- units_B[units_B$family == f, ]
  if (f == "egozi") d$unit <- gsub("_", " ", d$unit, fixed = TRUE)
  if (f == "gut_cell_atlas") d$unit <- sub(" \\([^)]*\\)$", "", d$unit)       # atlas donor ID only
  if (f == "gse238274") { d <- d[order(d$unit), ]; d$unit <- sub("^D", "Day ", d$unit) }
  d$row <- factor(d$unit, levels = rev(d$unit))
  pooled <- 100 * sum(d$called) / sum(d$n)
  ggplot(d, aes(rate, row)) +
    geom_vline(xintercept = pooled, color = "#777777", linetype = "dashed", linewidth = .35) +
    geom_point(aes(size = n), shape = if (f %in% names(B_SHAPE)) B_SHAPE[[f]] else 16,
               color = if (f %in% names(B_COL)) B_COL[[f]] else "#43A5B6", stroke = 0) +
    geom_text(aes(label = label), x = 108, hjust = 0, size = 7 / ggplot2::.pt, color = "black", family = "Arial") +
    scale_size_area(max_size = POINT_MAX, limits = c(0, max(units_B$n))) +
    scale_x_continuous(limits = c(-2, 153), breaks = c(0, 50, 100), labels = function(x) paste0(x, "%"), expand = c(0, 0)) +
    labs(x = NULL, y = NULL, title = STUDY[f], subtitle = CONTEXT[f]) +
    theme_fig2() +
    theme(panel.grid.major.y = element_blank(), axis.text.y = element_text(size = 7, color = "black"),
          plot.title = element_text(size = 8, face = "bold", margin = margin(0, 0, 1, 0)),
          plot.subtitle = element_text(size = PT, margin = margin(0, 0, 2, 0))) +
    (if (show_x) NULL else theme(axis.text.x = element_blank()))
}
pb <- setNames(lapply(FAMS, function(f) study_plot_B(f, show_x = f %in% c("hickey", "hio"))), FAMS)
plot_B <- wrap_plots(wrap_plots(pb$gut_cell_atlas, pb$hickey, ncol = 1, heights = c(10, 9)),
                     wrap_plots(pb$egozi, pb$gse238274, pb$hio, ncol = 1, heights = c(6, 4, 2)), ncol = 2)

point_key <- function(n, x, y, shape = 16, col = "black")
  grid.points(x = unit(x, "in"), y = unit(H - y, "in"), pch = shape,
              gp = gpar(col = col, fontsize = POINT_MAX * sqrt(n / max(units_B$n)) * ggplot2::.pt))

draw_panel <- function() {
  g <- patchworkGrob(plot_B)
  draw_at(g, 2.46, PLOT_Y, 4.94, PLOT_H)
  # one x-axis title centered under both columns' plots
  s <- panel_span(g, 4.94)
  draw_xtitle("Called EC (% of EC cells)", 2.46 + mean(s), PLOT_Y + PLOT_H)
  left <- 2.46; width <- 4.94; y <- KEY_Y
  # key line 1: the three sample kinds
  labels <- c("Donor or sample", "Day (one induced line)", "HIO sample group")
  ww <- vapply(labels, str_in, 0) + .16; gap <- .24
  x <- left + (width - sum(ww) - 2 * gap) / 2
  for (i in 1:3) {
    point_key(max(units_B$n), x + .055, y + .06, c(16, 17, 15)[i], c("#43A5B6", "#CCBB44", "#6943AC")[i])
    text_at(labels[i], x + .16, y)
    x <- x + ww[i] + gap
  }
  # key line 2: point sizes and the pooled-rate line
  title <- "EC cells (primary study label):"; sizes <- c(10, 50, 150); pl <- "Pooled rate for the study"
  sw <- vapply(as.character(sizes), str_in, 0) + .16
  total <- str_in(title) + .17 + sum(sw) + .17 * 3 + .13 + .16 + str_in(pl)
  x <- left + (width - total) / 2
  text_at(title, x, y + .19)
  x <- x + str_in(title) + .17
  for (i in 1:3) {
    point_key(sizes[i], x + .055, y + .24)
    text_at(as.character(sizes[i]), x + .16, y + .19)
    x <- x + sw[i] + .17
  }
  x <- x + .13
  grid.lines(x = unit(rep(x + .055, 2), "in"), y = unit(H - c(y + .16, y + .26), "in"),
             gp = gpar(col = "#777777", lty = 2, lwd = .35 * ggplot2::.pt))
  text_at(pl, x + .16, y + .19)
}
