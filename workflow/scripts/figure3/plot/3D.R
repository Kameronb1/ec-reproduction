# Figure 3D: SC-SIC against intestinal EC for every tested gene (donor-level limma-voom fit, study as a block): difference
# in log2 CPM against -log10 p. Genes with FDR >= 0.05 grey; open rings mark the four classifier genes. Seventeen
# genes are labeled, each with a thin black leader to its mark: eleven placed by ggrepel (seeded, with empty-label
# keep-out points), six at fixed positions. Two key rows below, centered under the plot panel.

dv <- read_table("panel_d_donor_fit.tsv.gz")
dv$difference <- dv$logFC; dv$nlp <- -log10(dv$p)
fdr_p <- max(dv$p[dv$FDR < 0.05])                                    # p at the FDR 0.05 line
# label placement, made at the drawn size: ggrepel's rows (the eleven gene labels at their marks, and 235 empty-label
# keep-out points over the FDR 0.05 label, the rings, the ring labels and two further areas, which steer them), the six
# fixed labels and their leaders
d_labels <- read_placement("positions_d_labels.tsv")
lab <- d_labels[d_labels$gene != "", ]; m <- match(lab$gene, dv$gene)                 # the labeled genes sit at their marks
stopifnot(!anyNA(m), max(abs(lab$difference - dv$difference[m]), abs(lab$nlp - dv$nlp[m])) < 1e-9)
d_fixed <- read_placement("positions_d_fixed_labels.tsv")
d_leaders <- read_placement("positions_d_fixed_leaders.tsv")
stopifnot(setequal(d_fixed$gene, d_leaders$gene), !anyDuplicated(d_leaders$gene))   # one leader per fixed label
x_rng <- range(dv$difference); y_rng <- range(dv$nlp)
D_XLIM <- x_rng + c(-1, 1) * (diff(x_rng) * 0.05)                                     # data range plus 5% each side,
D_YLIM <- c(y_rng[1] - 0.05 * diff(y_rng) - 0.75, y_rng[2] + 0.05 * diff(y_rng))   # and 0.75 more below for the ring labels

plot_D <- ggplot(dv, aes(difference, nlp)) +
  geom_hline(yintercept = -log10(fdr_p), linetype = "22", colour = "grey55", linewidth = 0.3) +
  geom_vline(xintercept = 0, colour = "grey75", linewidth = 0.3) +
  geom_point(data = dv[dv$call == "NO_CALL", ], colour = "grey78", size = 0.45, alpha = 0.6) +
  geom_point(data = dv[dv$call != "NO_CALL", ], aes(colour = call), size = 1.1, alpha = 0.85) +
  geom_point(data = dv[dv$label_gene, ], shape = 1, size = 2, colour = "grey10") +
  ggrepel::geom_text_repel(data = d_labels, aes(label = gene, fontface = face), size = MM, colour = "black", family = FONT,
                           segment.size = 0.2, segment.colour = "black", bg.colour = "white", bg.r = 0.12, seed = 20260925,
                           box.padding = unit(0.35, "lines"), min.segment.length = 0, force = 3, max.time = 5, max.iter = 1e5, max.overlaps = Inf) +
  annotate("text", x = max(dv$difference), y = -log10(fdr_p), label = "FDR 0.05", hjust = 1, vjust = -0.5, size = MM, colour = "black", family = FONT) +
  geom_segment(data = d_leaders, aes(x = x, y = y, xend = xend, yend = yend), inherit.aes = FALSE, colour = "black", linewidth = 0.2) +
  ggrepel::geom_text_repel(data = d_fixed, aes(x = x, y = y, label = gene, fontface = face), inherit.aes = FALSE, size = MM, family = FONT,
                           colour = "black", bg.colour = "white", bg.r = 0.12, force = 0, force_pull = 0, max.iter = 0,
                           min.segment.length = Inf, box.padding = 0, point.padding = 0, max.overlaps = Inf, seed = 1) +   # drawn where placed
  scale_colour_manual(values = c(SC_SIC_HIGHER = "#AA3377", NATIVE_EC_HIGHER = "#43A5B6"), guide = "none") +
  coord_cartesian(xlim = D_XLIM, ylim = D_YLIM, expand = FALSE) +
  labs(x = "SC-SIC minus intestinal EC (log2 CPM)", y = expression("−log10" ~ italic(p))) +
  theme_fig3()

PLOT_X <- 2.43; PLOT_W <- 2.29; PLOT_H <- 2.153
KEY_LEFT <- local({ pb <- panel_box(plot_D, PLOT_W, PLOT_H); PLOT_X + (pb$l + pb$r) / 2 - 1.15 })   # keys centered on the panel
# key item: an open ring drawn as ggplot draws shape 1 with the rings' stroke
ring_stroke <- unique(layer_data(plot_D, 5)$stroke)
ring_item <- function(label) list(label = label, w_sym = .14, draw = function(x, y)
  grid.points(x = unit(x + .04, "in"), y = unit(H - y - .04, "in"), pch = 1, size = unit(1, "char"),
              gp = gpar(col = "grey10", fontsize = 2 * .pt + ring_stroke * .stroke / 2, lwd = ring_stroke * .stroke / 2)))

draw_panel <- function() {
  draw_at(plot_D, PLOT_X, 1.727, PLOT_W, PLOT_H)
  key_row(list(dot_item("Higher in intestinal EC", "#43A5B6"), dot_item("Higher in SC-SIC", "#AA3377")), KEY_LEFT, 2.30, 4.005)
  key_row(list(dot_item("FDR ≥ 0.05", "#C7C7C7"), ring_item("Four-marker classifier gene")), KEY_LEFT, 2.30, 4.165)
}
