# Figure S12B: serotonin- minus dopamine-neuron gene score (adult brain) per sample: SC-SIC beside SC-beta and SC-alpha
# of each SC-islet sample (grey line: their range), with native intestinal EC, islet beta and brain neuron references. Brain values
# beyond the axis breaks are drawn at 0.3 scale. Open circles: single nuclei; filled: whole cells.

suppressPackageStartupMessages(library(data.table))

STUDY <- c(GSE151117 = "Augsornworawat 2020", GSE167880 = "Balboa 2022", GSE185036 = "Du 2022", GSE199636 = "Augsornworawat 2023",
           GSE202497 = "Zhu 2023", GSE226346 = "Siehler 2024")
ORDER_3A <- c("Augsornworawat 2020", "Balboa 2022", "Zhu 2023", "Siehler 2024", "Du 2022", "Augsornworawat 2023")   # Figure 3A order
SAMPLE <- c("S7.d18" = "day 18", "S7.d20" = "day 20", "S7.d25" = "day 25", s6_week2 = "S6 week 2", sc_islet_1 = "SC-islet 1", sc_islet_2 = "SC-islet 2")
SAMPLES <- list("Balboa 2022" = c("day 18", "day 20", "day 25"), "Augsornworawat 2023" = c("S6 week 2", "SC-islet 1", "SC-islet 2"))
SHOWN <- c("Augsornworawat 2023" = "Augsornworawat 2023 (nuclei)")
REF <- c("Primary gut EC", "HIO EC", "Induced EC", "Adult islet beta"); BR <- c("Brain serotonin neurons", "Brain dopamine neurons")

# one row per sample and group: SC-islet samples first (sample named where a study has several), then the references
U <- as.data.table(read_table("s9_brain/raphe_units.tsv"))
sc <- U[grepl("^SC-", grp)][, `:=`(cell = sub(" \\(nuclei study\\)$", "", grp), sample = SAMPLE[sub("^[^|]+\\|[^_]+__([^|]+)\\|.*$", "\\1", group)])]
sc[is.na(sample), sample := ""]
Bt <- rbind(sc[, .(study, sample, cell, v = ser_minus_da)], U[grp %in% REF, .(study, sample = "", cell = grp, v = ser_minus_da)],
            U[grp %in% BR, .(study, sample = "", cell = grp, v = ser_minus_da)])
Bt[, nuclei := study %in% NUC_STUDY]
Bt[, sname := STUDY[study]]
Bt[, rowkey := ifelse(is.na(sname), cell, paste(sname, sample))]
# rows, top to bottom: brain references, intestinal EC and islet references, then SC-islet samples by study (bold study name over
# plain sample names when a study has several samples)
lay <- data.table(rowkey = c(BR, REF), label = NAMES[c(BR, REF)], face = "plain")
for (s in ORDER_3A) {
  sh <- if (s %in% names(SHOWN)) SHOWN[[s]] else s
  if (is.null(SAMPLES[[s]])) lay <- rbind(lay, data.table(rowkey = paste(s, ""), label = sh, face = "bold"))
  else lay <- rbind(lay, data.table(rowkey = NA_character_, label = sh, face = "bold"),
                    data.table(rowkey = paste(s, SAMPLES[[s]]), label = SAMPLES[[s]], face = "plain"))
}
NS <- nrow(lay); lay[, y := NS - seq_len(NS) + 1]
Bt[, y := lay$y[match(rowkey, lay$rowkey)]]
Y_SEP <- c(NS - 2 + .5, NS - 6 + .5)                    # below the brain rows; below the reference rows
XB <- 1.17; fx <- function(x) ifelse(abs(x) > 2, sign(x) * (1.22 + (abs(x) - 3.3) * .3), x)   # brain values at 0.3 scale past the breaks
Bt[, x := fx(v)]
Bt[, `:=`(col = factor(cell, c("SC-SIC", "SC-beta", "SC-alpha", REF, BR)), shp = factor(ifelse(nuclei, "Single nuclei", "Whole cells"), names(SHP)))]
seg <- Bt[!is.na(sname), .(lo = min(v), hi = max(v)), by = y]
Bp <- Bt[order(nuclei)]                                 # nucleus points drawn last, on top
BCOL <- c("SC-SIC" = SIC, "SC-beta" = BETA, "SC-alpha" = ALPHA, "Primary gut EC" = "#43A5B6", "HIO EC" = "#6943AC",
          "Induced EC" = "#CCBB44", "Adult islet beta" = "#1C3A63", "Brain serotonin neurons" = "#222222", "Brain dopamine neurons" = "#888888")
XL <- c(-1.80, 1.80); YL <- c(.4, NS + .6)
# the row labels sit in the left plot margin, sized to the widest label measured on a 72 dpi device
width_72 <- function(s, bold) { pushViewport(viewport(gp = gpar(fontfamily = FONT, fontsize = 7.5, fontface = if (bold) "bold" else "plain")))
  on.exit(popViewport()); convertWidth(stringWidth(s), "in", valueOnly = TRUE) }
png(tempfile(fileext = ".png"), type = "cairo"); LW <- max(mapply(function(s, f) width_72(s, f == "bold"), lay$label, lay$face)); invisible(dev.off())
plot_B <- ggplot(Bt, aes(x, y)) +
  geom_vline(xintercept = 0, colour = "grey70", linewidth = .3) +
  geom_hline(yintercept = Y_SEP, colour = "grey80", linewidth = .3) +
  annotate("rect", xmin = c(-XB, XB) - .025, xmax = c(-XB, XB) + .025, ymin = -Inf, ymax = Inf, fill = "white") +
  annotate("segment", x = rep(c(-XB, XB), each = 2) + c(-.04, 0, -.04, 0), xend = rep(c(-XB, XB), each = 2) + c(0, .04, 0, .04),
           y = YL[1] - .35, yend = YL[1] + .25, linewidth = .35) +
  geom_segment(data = seg, aes(x = lo, xend = hi, y = y, yend = y), inherit.aes = FALSE, colour = "grey75", linewidth = .4) +
  geom_point(data = Bp, aes(colour = col, shape = shp), size = 1.5, stroke = .6, fill = "white", position = position_jitter(height = .12, width = 0, seed = 1)) +
  annotate("text", x = XL[1] - .06, y = lay$y, label = lay$label, fontface = lay$face, hjust = 1, vjust = .5, size = 7.5 / .pt, family = FONT, colour = "black") +
  scale_colour_manual(values = BCOL, breaks = c("SC-SIC", "SC-beta", "SC-alpha"), name = NULL) +
  scale_shape_manual(values = c("Whole cells" = 16, "Single nuclei" = 21), name = NULL, drop = FALSE) +   # 21 + white fill: an open circle that hides what is under it
  scale_x_continuous(breaks = fx(c(-5, -4, -1, -.5, 0, .5, 1, 4, 5)), labels = mlab(c(-5, -4, -1, -.5, 0, .5, 1, 4, 5)), expand = c(0, 0)) +
  scale_y_continuous(breaks = NULL) + coord_cartesian(xlim = XL, ylim = YL, expand = FALSE, clip = "off") +
  labs(x = bquote("Serotonin- minus dopamine-neuron genes (adult brain, log"[2] * ")"), y = NULL) + panel_theme("x", legend = "bottom") +
  guides(colour = guide_legend(order = 1, nrow = 1, override.aes = list(size = 1.8)), shape = guide_legend(order = 2, nrow = 1, override.aes = list(size = 1.8))) +
  theme(legend.box = "vertical", legend.box.just = "left", legend.justification = "left", legend.spacing.y = unit(0, "pt"), legend.margin = margin(0, 0, 0, 0),
        axis.text.y = element_blank(), axis.title.x = element_text(hjust = 1), plot.margin = margin(2, 6, 2, 72 * (LW + .10)))

draw_panel <- function() plot_panel(plot_B, 120, 1572, 3.72, 3.30)
