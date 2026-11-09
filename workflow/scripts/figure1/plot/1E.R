# Figure 1E: islet hormone transcripts in the cells that make the serotonin enzyme. Of the cells with TPH1 detected
# after SoupX correction, the share that also carry INS only, GCG only, SST only, two or more of them, or none. Intestinal
# EC is the mean of its three studies; induced EC and HIO graft EC are one study each; SC-SIC is one row per study (a
# white circle after a single-nucleus study's name), sorted by the share with any hormone. Row labels give the TPH1+ cells;
# key below.

EH <- read_table("hormones_in_tph1_cells.tsv")
EH$grp <- display_group(EH$grp)
HCAT <- c(INS = "INS only", GCG = "GCG only", SST = "SST only", two_plus = "Two or more", none = "No islet hormone")
HCOL <- c(INS = COL[["beta"]], GCG = COL[["alpha"]], SST = COL[["delta"]], two_plus = "#333333", none = "#EBEBEB")
# The TPH1+ cells must be panel D's (its cells times its TPH1 share), and each row's categories must sum to 1.
KD <- serotonin_rows()
k <- match(paste(EH$study, EH$grp), paste(KD$study, KD$grp))
stopifnot(!anyNA(k), EH$tph1_cells == round(KD$TPH1[k] * KD$cells[k]), abs(rowSums(EH[, names(HCAT)]) - 1) <= 1e-9)
gut <- EH[EH$grp == "Intestinal EC", ]
sic <- EH[EH$grp == "SC-SIC", ]; sic <- sic[order(-(1 - sic$none)), ]
ER <- rbind(
  data.frame(lab = sprintf("Intestinal EC, %d studies  (%s)", nrow(gut), format(sum(gut$tph1_cells), big.mark = ",")), block = "Intestine",
             t(colMeans(gut[, names(HCAT)]))),
  do.call(rbind, lapply(c("Induced EC", "HIO graft EC"), function(g) { e <- EH[EH$grp == g, ]
    data.frame(lab = sprintf("%s  (%s)", g, format(e$tph1_cells, big.mark = ",")), block = "Intestine", e[, names(HCAT)]) })),
  data.frame(lab = sprintf("%s%s  (%s)", sic$study, ifelse(sic$study %in% NUC, " \u25cb", ""), format(sic$tph1_cells, big.mark = ",", trim = TRUE)),
             block = "SC-SIC", sic[, names(HCAT)]))
# stacked bars: each category's tile spans its cumulative share
EB <- do.call(rbind, lapply(seq_len(nrow(ER)), function(i) { v <- unlist(ER[i, names(HCAT)]); hi <- cumsum(v); lo <- hi - v
  data.frame(lab = ER$lab[i], block = ER$block[i], key = names(HCAT), mid = (lo + hi) / 2, wd = v) }))
EB$lab <- factor(EB$lab, levels = rev(ER$lab)); EB$block <- factor(EB$block, levels = c("Intestine", "SC-SIC"))
EB$key <- factor(EB$key, levels = names(HCAT))

plot_E <- ggplot(EB, aes(x = mid, y = lab)) + geom_tile(aes(width = wd, fill = key), height = .66) +
  facet_grid(block ~ ., scales = "free_y", space = "free_y") +
  scale_fill_manual(values = HCOL, labels = HCAT, name = NULL) +
  scale_x_continuous(breaks = c(0, .25, .5, .75, 1), labels = c("0", "25", "50", "75", "100%")) +
  coord_cartesian(xlim = c(0, 1), expand = FALSE, clip = "off") +
  labs(x = "TPH1+ cells with islet hormone transcripts detected after ambient-RNA correction (%)", y = NULL) +
  theme_classic(base_size = 7, base_family = "Arial") +
  theme(strip.background = element_blank(), strip.text.y = element_text(size = 6.5, colour = SEC, angle = 0, hjust = 0, margin = margin(0, 0, 0, 4)),
        axis.text.y = element_text(size = 6.5, colour = "black"), axis.text.x = element_text(size = 6.5, colour = SEC),
        panel.spacing.y = unit(.05, "in"), panel.grid.major.x = element_line(colour = "grey92", linewidth = .25),
        legend.position = "bottom", legend.justification = "center", legend.location = "plot",
        legend.text = element_text(size = 6.5, colour = SEC), legend.key.size = unit(.10, "in"), legend.spacing.x = unit(.06, "in"),
        legend.margin = margin(0, 0, 0, 0), legend.box.margin = margin(-3, 0, 0, 0), plot.margin = margin(14, 9, 2, 10.5))

# The plot fills row E. A 9 pt row above the panels holds a 7 pt header over the row labels, its right edge on the
# labels' (tick 1.75 pt + label margin 1.4 pt left of the axis), 2 pt above the panel.
draw_E <- function() in_row(ROW_E, function(PH) {
  gt <- ggplotGrob(plot_E)
  top <- min(gt$layout$t[grepl("^panel", gt$layout$name)])
  gt <- gtable::gtable_add_rows(gt, unit(9, "pt"), top - 1)
  gt <- gtable::gtable_add_grob(gt, textGrob("(TPH1+ cells)", x = unit(1, "npc") - unit(3.15, "pt"), y = unit(2, "pt"), hjust = 1, vjust = 0,
                                             gp = gpar(fontfamily = "Arial", fontsize = 7)),
                                t = top, l = min(gt$layout$l[grepl("^axis-l", gt$layout$name)]), clip = "off", name = "row-label-header")
  grid.draw(gt)
})

draw_panel <- function() draw_after_72dpi(draw_E)
