# Figure 2D: dot plot of the 20 EC-enriched genes (rows, in panel rank order) by dataset (columns): color = mean
# expression in EC minus other EEC; size = share of the dataset's EC cells with at least one read of the gene.
# Color bar and size key below.

PLOT_Y <- 1954 / 600; PLOT_H <- 2249 / 600; KEY_Y <- 7.255   # top and height of the plot (whole pixels), top of the key

effects_D <- read_table("D_effect_heatmap.tsv")[c("family_id", "gene_symbol", "effect", "panel_rank")]
# McCauley 2023 HIO grafts: the two HIO sample groups' effects are averaged, as donors and samples are in the tissue
# studies
hio_D <- read_table("hio_D_hio_effects.tsv")
hio_D <- aggregate(effect ~ gene_symbol, hio_D[hio_D$sample != "all", ], mean)
effects_D <- rbind(effects_D, data.frame(family_id = "SF_GSE214852_HIO_GRAFT", gene_symbol = hio_D$gene_symbol,
                                         effect = hio_D$effect,
                                         panel_rank = effects_D$panel_rank[match(hio_D$gene_symbol, effects_D$gene_symbol)]))
# Columns: four datasets labeled by the four-marker classifier and three with primary study labels (the seven gene-ranking
# datasets, from six studies: Elmentaite 2020 and 2021 are one study), then Singh 2024 and McCauley 2023. x positions
# leave small gaps between the groups.
cols_D <- data.frame(
  family_id = c("SF_GSE185224_ADULT_REGIONAL_INTESTINE", "SF_ELMENTAITE2020_DEVGUT", "SF_GSE125970_ADULT_INTESTINE",
                "SF_YU2021_ENDODERM", "SF_ELMENTAITE2021_GUT_ATLAS", "SF_HICKEY2023_ADULT_EEC",
                "SF_EGOZI2021_DEVELOPMENTAL_EEC", "SF_INTESTINAL_GSE238274", "SF_GSE214852_HIO_GRAFT"),
  x = c(1, 2, 3, 4, 5.3, 6.3, 7.3, 8.6, 9.6))
cols_D$label <- SN[sub("^SF_YU2021_ENDODERM$", "SF_YU2021_ENDODERM:fetal", cols_D$family_id), "author_year"]
genes_D <- unique(effects_D$gene_symbol[order(effects_D$panel_rank)])
# detection: share of EC cells with >= 1 read (McCauley 2023: mean of the two HIO sample groups, as for the color);
# the detection table's effects must equal the drawn ones
det_D <- read_table("D_detection.tsv")
k <- match(paste(effects_D$family_id, effects_D$gene_symbol), paste(det_D$family_id, det_D$gene_symbol))
stopifnot(!anyNA(k), max(abs(effects_D$effect - det_D$effect[k])) < 1e-8)
effects_D$det_ec <- det_D$det_ec[k]
D_MAX <- 3.6                                                  # dot diameter (mm) for 100% of EC cells
effects_D$x <- cols_D$x[match(effects_D$family_id, cols_D$family_id)]
effects_D$y <- match(effects_D$gene_symbol, rev(genes_D))
D_LO <- "#2F6DA8"; D_HI <- "#8E3D29"; D_MID <- "#F2F2F2"; D_LIM <- c(-2, 6)       # values beyond 6 are drawn at 6
ny <- length(genes_D)
# column-group brackets above the dots
groups_D <- data.frame(x0 = c(.6, 4.9), x1 = c(4.4, 10.0), lab = c("Four-marker\nclassifier calls", "Primary study\nlabels"))
group2_D <- data.frame(x0 = .6, x1 = 7.7, lab = "Gene-ranking datasets")

plot_D <- ggplot(effects_D, aes(x, y, fill = effect)) +
  geom_point(aes(size = det_ec), shape = 21, color = "black", stroke = .15) +
  geom_segment(data = groups_D, aes(x = x0, xend = x1, y = ny + .85, yend = ny + .85), inherit.aes = FALSE,
               linewidth = .35, color = "black") +
  geom_text(data = groups_D, aes(x = (x0 + x1) / 2, y = ny + 1.1, label = lab), inherit.aes = FALSE, vjust = 0,
            lineheight = .9, size = 7 / ggplot2::.pt, family = "Arial", color = "black") +
  geom_segment(data = group2_D, aes(x = x0, xend = x1, y = ny + 2.95, yend = ny + 2.95), inherit.aes = FALSE,
               linewidth = .35, color = "black") +
  geom_text(data = group2_D, aes(x = (x0 + x1) / 2, y = ny + 3.2, label = lab), inherit.aes = FALSE, vjust = 0,
            size = 7 / ggplot2::.pt, family = "Arial", color = "black") +
  scale_fill_gradient2(low = D_LO, mid = D_MID, high = D_HI, midpoint = 0, limits = D_LIM, oob = scales::squish) +
  scale_size_area(max_size = D_MAX, limits = c(0, 1)) +
  scale_x_continuous(breaks = cols_D$x, labels = cols_D$label, limits = c(.5, 10.35), expand = c(0, 0)) +
  scale_y_continuous(breaks = seq_along(genes_D), labels = rev(genes_D), limits = c(.5, ny + 4.1), expand = c(0, 0)) +
  labs(x = NULL, y = NULL) +
  theme_fig2() +
  theme(panel.grid = element_blank(), panel.grid.major = element_blank(), panel.grid.minor = element_blank(),
        axis.line = element_blank(), axis.line.x = element_blank(), axis.ticks = element_blank(),
        axis.ticks.x = element_blank(),
        axis.text.x = element_text(angle = 55, hjust = 1, vjust = 1, size = PT, color = "black"),
        axis.text.y = element_text(face = "italic", size = 8, color = "black"), plot.margin = margin(1, 2, 0, 0),
        plot.background = element_blank())      # no background: the plot's top reaches up beside panel B's key

draw_panel <- function() {
  draw_at(plot_D, 5.15, PLOT_Y, 2.24, PLOT_H)
  # color bar: title, ramp (a 240-pixel raster stretched across the bar), tick labels, centered under the dots
  left <- 5.15; width <- 2.24; y <- KEY_Y
  title <- "EC minus other EEC (log expression)"
  text_at(title, left + (width - str_in(title)) / 2, y)
  bw <- 1.50; x0 <- left + (width - bw) / 2
  v <- seq(D_LIM[1], D_LIM[2], length.out = 240)
  ramp <- ifelse(v < 0, grDevices::rgb(grDevices::colorRamp(c(D_LO, D_MID))(pmin(1, pmax(0, (v - D_LIM[1]) / (0 - D_LIM[1])))) / 255),
                 grDevices::rgb(grDevices::colorRamp(c(D_MID, D_HI))(pmin(1, pmax(0, v / D_LIM[2]))) / 255))
  grid.raster(matrix(ramp, nrow = 1), x = unit(x0, "in"), y = unit(H - y - .15, "in"), width = unit(bw, "in"),
              height = unit(.09, "in"), just = c("left", "top"), interpolate = FALSE)
  for (t in c(-2, 0, 2, 4, 6))
    grid.text(if (t == 6) "6+" else as.character(t), x = unit(x0 + bw * (t - D_LIM[1]) / diff(D_LIM), "in"),
              y = unit(H - y - .27, "in"), just = c("centre", "top"), gp = gpar(fontfamily = "Arial", fontsize = PT, col = "black"))
  # size key under the color bar
  y <- y + .40; ttl <- "EC cells detecting it:"; sz <- c(.25, .5, 1); lb <- c("25%", "50%", "100%")
  tw <- str_in(ttl) + .08 + sum(vapply(lb, str_in, 0) + .22) - .07
  x <- left + (width - tw) / 2; text_at(ttl, x, y); x <- x + str_in(ttl) + .08
  for (i in 1:3) {
    grid.points(unit(x + .06, "in"), unit(H - y - .05, "in"), pch = 21, gp = gpar(col = "black", fill = "#8C8C8C", lwd = .15 * ggplot2::.pt, fontsize = D_MAX * sqrt(sz[i]) * ggplot2::.pt))
    text_at(lb[i], x + .15, y); x <- x + str_in(lb[i]) + .22
  }
}
