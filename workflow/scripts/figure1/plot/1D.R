# Figure 1D: transcripts for making and storing serotonin in the same cell. For each study and cell group (groups of
# at least 30 cells), the share of cells with TPH1, DDC and SLC18A1 detected after SoupX correction: one dot per
# study (single-nucleus studies white with a ring in the group's color), a grey line spans
# the studies. The last column stacks, per cell group, the share of cells with all three, two, one or none of the genes
# (mean of the studies' shares). Row labels give each group's number of studies and cells; key below.

KD <- serotonin_rows()
KGRP <- c("Intestinal EC", "Induced EC", "HIO graft EC", "Intestinal, other EEC", "SC-SIC", "SC-alpha", "SC-beta", "SC-delta",
          "Adult alpha", "Adult beta", "Adult delta")
KBLOCK <- setNames(c(rep("Intestine", 4), rep("SC-islet", 4), rep("Adult islet", 3)), KGRP)
GCOL <- c("Intestinal EC" = COL[["gut_ec"]], "Induced EC" = COL[["induced_ec"]], "HIO graft EC" = COL[["hio_ec"]], "Intestinal, other EEC" = "#8A9AA6",
          "SC-SIC" = COL[["sc_sic"]], "SC-beta" = COL[["beta"]], "SC-alpha" = COL[["alpha"]], "SC-delta" = COL[["delta"]],
          "Adult beta" = COL[["beta"]], "Adult alpha" = COL[["alpha"]], "Adult delta" = COL[["delta"]])
MEAS <- c(TPH1 = "TPH1\n(makes 5-HTP)", DDC = "DDC\n(makes serotonin)", SLC18A1 = "SLC18A1 / VMAT1\n(EC storage)",
          count = "How many of the three\n(mean of studies)")
# long table: one row per study, group and transcript; row label "group  (studies; cells)"
KL <- do.call(rbind, lapply(c("TPH1", "DDC", "SLC18A1"), function(k)
  data.frame(study = KD$study, grp = KD$grp, cells = KD$cells, measure = MEAS[[k]], share = KD[[k]])))
KL$measure <- factor(KL$measure, levels = MEAS)
klab <- vapply(KGRP, function(g) { d <- KD[KD$grp == g, ]; sprintf("%s  (%d; %s)", g, nrow(d), trimws(format(sum(d$cells), big.mark = ","))) }, "")
KL$row <- factor(klab[KL$grp], levels = rev(klab))
KL$block <- factor(KBLOCK[KL$grp], levels = c("Intestine", "SC-islet", "Adult islet"), labels = c("Intestine", "SC-islet", "Adult\nislet"))
KL$assay <- factor(ifelse(KL$study %in% NUC, "single-nucleus", "single-cell"), levels = c("single-cell", "single-nucleus"))
KR <- do.call(rbind, lapply(split(KL, list(KL$row, KL$measure), drop = TRUE), function(e)     # range across studies
  data.frame(row = e$row[1], measure = e$measure[1], block = e$block[1], lo = min(e$share), hi = max(e$share))))
KL$fillkey <- ifelse(KL$assay == "single-nucleus", "nucleus", KL$grp)

# last column: share of cells with 3, 2, 1, 0 of the genes, mean of the drawn studies; stacked from 0 (all three first).
# The gene-count table must cover the same cells as the co-detection table, with the same all-three share.
KC <- read_table("serotonin_gene_count.tsv")
KC$grp <- display_group(KC$grp)
KC <- KC[paste(KC$study, KC$grp) %in% paste(KD$study, KD$grp), ]
k <- match(paste(KC$study, KC$grp), paste(KD$study, KD$grp))
stopifnot(nrow(KC) == nrow(KD), abs(KC$n3 - KD$all_three[k]) <= 1e-12, KC$cells == KD$cells[k])
NCOL <- c(n3 = "#333333", n2 = "#7A7A7A", n1 = "#B9B9B9", n0 = "#EBEBEB")
KM <- do.call(rbind, lapply(split(KC, KC$grp), function(e)
  data.frame(grp = e$grp[1], n3 = mean(e$n3), n2 = mean(e$n2), n1 = mean(e$n1), n0 = mean(e$n0))))
KB <- do.call(rbind, lapply(seq_len(nrow(KM)), function(i) { v <- unlist(KM[i, names(NCOL)]); hi <- cumsum(v); lo <- hi - v
  data.frame(grp = KM$grp[i], key = names(v), mid = (lo + hi) / 2, wd = v) }))
KB$row <- factor(klab[KB$grp], levels = levels(KL$row))
KB$block <- factor(KBLOCK[KB$grp], levels = c("Intestine", "SC-islet", "Adult islet"), labels = levels(KL$block))
KB$measure <- factor(MEAS[["count"]], levels = MEAS)

plot_D <- ggplot(KL, aes(share, row)) +
  geom_tile(data = KB, aes(x = mid, y = row, width = wd, fill = key), height = .66, inherit.aes = FALSE) +
  geom_linerange(data = KR, aes(xmin = lo, xmax = hi, y = row), inherit.aes = FALSE, colour = "grey85", linewidth = .25) +
  geom_point(aes(fill = fillkey, shape = assay), colour = "black", size = 1.45, stroke = .2, show.legend = c(fill = FALSE, shape = TRUE)) +
  geom_point(data = KL[KL$assay == "single-nucleus", ], aes(colour = grp), shape = 1, size = .95, stroke = .45, show.legend = FALSE) +
  scale_colour_manual(values = GCOL, guide = "none") +
  scale_fill_manual(values = c(GCOL, "nucleus" = "white", NCOL), breaks = names(NCOL), labels = c("all 3", "2", "1", "0 of the genes"), name = NULL) +
  scale_shape_manual(values = c("single-cell" = 21, "single-nucleus" = 21), name = NULL,
                     labels = c("single-cell" = "one study (single-cell)", "single-nucleus" = "one study (single-nucleus)")) +
  guides(shape = guide_legend(order = 1, override.aes = list(fill = c("grey45", "white"), colour = "black", size = 1.6, stroke = .3)),
         fill = guide_legend(order = 2)) +
  facet_grid(block ~ measure, scales = "free_y", space = "free_y") +
  scale_x_continuous(limits = c(0, 1.0001), breaks = c(0, .5, 1), labels = c("0", "50", "100%"), expand = expansion(add = c(.03, .05))) +
  labs(x = "Cells with the transcripts detected after ambient-RNA correction (%)", y = NULL) +
  theme_classic(base_size = 7, base_family = "Arial") +
  theme(strip.background = element_blank(), strip.text.x = element_text(face = "bold", size = 7, lineheight = .9),
        strip.text.y = element_text(size = 6.5, colour = SEC), axis.text.y = element_text(size = 6.5, colour = "black"),
        axis.text.x = element_text(size = 6.5, colour = SEC),
        panel.spacing.x = unit(.10, "in"), panel.spacing.y = unit(.05, "in"), panel.grid.major.x = element_line(colour = "grey92", linewidth = .25),
        legend.position = "bottom", legend.justification = "center", legend.location = "plot", legend.box = "horizontal",
        legend.text = element_text(size = 6.5, colour = SEC), legend.key.size = unit(.10, "in"), legend.spacing.x = unit(.06, "in"),
        legend.margin = margin(0, 0, 0, 0), legend.box.margin = margin(-3, 0, 0, 0),
        plot.margin = margin(26 - 169 / 600 * 72.27, 9, 2, 10.5))    # top margin: 26 pt less 169 px (in points of 1/72.27 in)

# The plot fills row D. A 7 pt header names the row labels' counts, in the empty cell left of the first strip: its
# right edge on the labels' (tick 1.75 pt + label margin 1.4 pt left of the axis), 2 pt above the panel.
draw_D <- function() in_row(ROW_D, function(PH) {
  gt <- ggplotGrob(plot_D)
  gt <- gtable::gtable_add_grob(gt, textGrob("(studies; cells)", x = unit(1, "npc") - unit(3.15, "pt"), y = unit(2, "pt"),
                                             hjust = 1, vjust = 0, gp = gpar(fontfamily = "Arial", fontsize = 7)),
                                t = min(gt$layout$t[grepl("^strip-t", gt$layout$name)]),
                                l = min(gt$layout$l[grepl("^axis-l", gt$layout$name)]), clip = "off", name = "row-label-header")
  grid.draw(gt)
})

draw_panel <- function() draw_after_72dpi(draw_D)
