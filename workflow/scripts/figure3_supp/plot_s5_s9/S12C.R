# Figure S12C: each of the 81 serotonin-neuron genes (higher in adult brain serotonin than dopamine neurons), SC-SIC minus
# SC-beta (x) against SC-SIC minus intestinal EC (y); shaded: 2-fold or more above both. Those genes are labeled, and TPH2 and
# SLC6A4 in grey, at fixed label positions (s9c_label_positions.tsv), with a leader where a label sits away from its point.

suppressPackageStartupMessages(library(data.table))

G <- as.data.table(read_table("s9_brain/ser_genes_in_sc_sic.tsv"))
G[, both := sic_minus_beta >= 1 & sic_minus_gut >= 1]
KL <- c(sprintf("2-fold or more above both (%d)", sum(G$both)), sprintf("Other serotonin-neuron genes (%d)", sum(!G$both)))
G[, k := factor(ifelse(both, KL[1], KL[2]), KL)]
NAMED <- c("TPH2", "SLC6A4")                            # defining serotonin-neuron genes the legend cites, labeled in grey
POS <- read.delim(file.path(dirname(PLOT_SCRIPT), "s9c_label_positions.tsv"))   # label centers, in data units
G[POS, on = "gene", `:=`(label_x = i.label_x, label_y = i.label_y)]
G[, labeled := both | gene %in% NAMED]
stopifnot(setequal(POS$gene, G$gene[G$labeled]))
PW_C <- 2.86; PH_C <- 1.98                              # plot size (in), so data units map to inches exactly
XR <- range(G$sic_minus_beta) + c(-1, 1) * .05 * diff(range(G$sic_minus_beta))
YR <- range(G$sic_minus_gut) + c(-1, 1) * .05 * diff(range(G$sic_minus_gut))
toX <- function(x) (x - XR[1]) / diff(XR) * PW_C; toY <- function(y) (y - YR[1]) / diff(YR) * PH_C
frX <- function(u) XR[1] + u / PW_C * diff(XR); frY <- function(u) YR[1] + u / PH_C * diff(YR)
LB <- G[labeled == TRUE, .(gene, col = ifelse(gene %in% NAMED, "grey30", "black"), px = toX(sic_minus_beta), py = toY(sic_minus_gut),
                            lx = toX(label_x), ly = toY(label_y))]
# label half-widths and half-height (in), measured on a 72 dpi device
png(tempfile(fileext = ".png"), type = "cairo"); pushViewport(viewport(gp = gpar(fontfamily = FONT, fontsize = 7, fontface = "italic")))
LB[, `:=`(hw = vapply(gene, function(s) convertWidth(stringWidth(s), "in", valueOnly = TRUE), 0) / 2 + .015,
          hh = convertHeight(stringHeight("ABC0"), "in", valueOnly = TRUE) / 2 + .015)]
popViewport(); invisible(dev.off())
# leader from the point's edge to the nearest point of its label box, drawn when the label is more than 0.05 in away
LB[, `:=`(ex = pmin(pmax(px, lx - hw), lx + hw), ey = pmin(pmax(py, ly - hh), ly + hh))]
LB[, len := sqrt((ex - px)^2 + (ey - py)^2)]
LB[, `:=`(sx = px + (ex - px) / len * .03, sy = py + (ey - py) / len * .03)]
LBl <- LB[len > .05]
gc_ <- ggplot(G, aes(sic_minus_beta, sic_minus_gut)) +
  geom_hline(yintercept = 0, colour = "grey70", linewidth = .3) + geom_vline(xintercept = 0, colour = "grey70", linewidth = .3) +
  annotate("rect", xmin = 1, xmax = Inf, ymin = 1, ymax = Inf, fill = SIC, alpha = .07) +
  geom_point(aes(colour = k), size = 1.2) +
  annotate("segment", x = frX(LBl$sx), y = frY(LBl$sy), xend = frX(LBl$ex), yend = frY(LBl$ey), colour = "grey40", linewidth = .25) +
  annotate("text", x = frX(LB$lx), y = frY(LB$ly), label = LB$gene, colour = LB$col, size = 7 / .pt, family = FONT, fontface = "italic", hjust = .5, vjust = .5) +
  scale_colour_manual(values = setNames(c(SIC, "#B5B5B5"), KL), name = NULL) +
  scale_x_continuous(labels = mlab, limits = XR, expand = c(0, 0)) + scale_y_continuous(labels = mlab, limits = YR, expand = c(0, 0)) +
  labs(x = bquote("SC-SIC minus SC-beta " * "(log"[2] * ")"), y = bquote("SC-SIC minus intestinal EC " * "(log"[2] * ")")) +
  panel_theme("both", legend = "bottom") + guides(colour = guide_legend(ncol = 1, override.aes = list(size = 1.8))) +
  theme(axis.line.y = element_line(colour = "black", linewidth = .35), axis.ticks.y = element_line(colour = "black", linewidth = .3),
        legend.justification = "left", legend.margin = margin(0, 0, 0, 0))
plot_C <- wrap_plots(gc_) + plot_layout(widths = unit(PW_C, "in"), heights = unit(PH_C, "in"))

draw_panel <- function() plot_panel(plot_C, 2472, 1572, 3.38, 3.30)
