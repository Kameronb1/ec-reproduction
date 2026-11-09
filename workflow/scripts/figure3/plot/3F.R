# Figure 3F: the three gene groups on Figure 3D's x axis (SC-SIC minus intestinal EC). Rows = in how many of the three intestinal
# studies the gene is 2-fold up in EC: G1 and G2 genes sit in the top row, G3 genes in the bottom row and the mixed genes
# in between, each spread at random within its row (seeded). Ten genes are labeled. The y-axis title is set above the
# row labels; each label's first line (the number of studies) is bold.

GKEY <- c(G1_EC_GENES_SC_SIC_LACK = "G1", G2_EC_GENES_SC_SIC_HAVE = "G2", G3_SC_SIC_OWN_PROGRAM = "G3")
F_ROWS <- c("0 or 1 of 3", "2 of 3 (mixed)", "All 3 intestinal studies")                     # bottom to top
gr <- read_table("panel_f_gene_groups.tsv.gz"); mx <- read_table("panel_f_mixed_genes.tsv")
dv <- read_table("panel_d_donor_fit.tsv.gz")
fd <- rbind(data.frame(gene = gr$gene, g = unname(GKEY[gr$group]), row = ifelse(gr$gut_up_of_3 == 3, F_ROWS[3], F_ROWS[1])),
            data.frame(gene = mx$gene, g = "Mixed", row = F_ROWS[2]))
fd$x <- dv$logFC[match(fd$gene, dv$gene)]
fd$row <- factor(fd$row, F_ROWS); fd$g <- factor(fd$g, names(GROUP_COLS))
set.seed(20261004); fd$yj <- as.integer(fd$row) + runif(nrow(fd), -0.24, 0.24)          # spread within the row
n <- table(fd$g)
F_YLAB <- c(sprintf("0 or 1 of 3\nG3: %d", n[["G3"]]), sprintf("2 of 3\n(mixed)\n%d, no group", n[["Mixed"]]),
            sprintf("All 3\nG1: %d\nG2: %d", n[["G1"]], n[["G2"]]))
# The row labels are drawn twice over the axis text, which stays (transparent) so the layout does not move: the first
# line in bold, the lines under it plain, each with the other lines blank. The axis text ends 3.36 pt left of the panel
# (measured: drawn plain, the overlay equals the axis text to the pixel); the y range is 0.55-3.45.
F_HEAD_FACE <- "bold"
f_row_label <- function(i, first) {
  ln <- strsplit(F_YLAB[i], "\n", fixed = TRUE)[[1]]; ln[if (first) -1 else 1] <- " "
  annotation_custom(textGrob(paste(ln, collapse = "\n"), x = unit(0, "npc") - unit(3.36, "pt"), y = unit((i - .55) / 2.9, "npc"),
                             hjust = 1, gp = gpar(fontfamily = FONT, fontsize = PT, col = "black", lineheight = .85,
                                                  fontface = if (first) F_HEAD_FACE else "plain")))
}
# label centers and leaders (data units), placed at the drawn size
f_labels <- read_placement("positions_f_labels.tsv"); f_leaders <- read_placement("positions_f_leaders.tsv")

plot_F <- ggplot(fd, aes(x, yj, colour = g)) +
  geom_vline(xintercept = 0, colour = "grey55", linewidth = 0.25) +
  geom_point(size = 0.95, alpha = 0.9, stroke = 0) +
  geom_segment(data = f_leaders, aes(x = x, y = y, xend = xend, yend = yend), inherit.aes = FALSE, colour = "grey40", linewidth = 0.2) +
  geom_text(data = f_labels, aes(x = x, y = y, label = gene), inherit.aes = FALSE, size = MM, family = FONT, fontface = "italic", colour = "black") +
  scale_colour_manual(values = GROUP_COLS, guide = "none") +
  scale_y_continuous(breaks = 1:3, labels = F_YLAB, limits = c(0.55, 3.45), expand = c(0, 0)) +
  scale_x_continuous(limits = c(-15.2, 9.6), breaks = c(-10, -5, 0, 5)) +
  lapply(1:3, f_row_label, first = TRUE) + lapply(1:3, f_row_label, first = FALSE) + coord_cartesian(clip = "off") +
  labs(x = "SC-SIC minus intestinal EC (log2 CPM)", y = NULL) +
  theme_fig3() + theme(panel.grid.major.y = element_blank(), axis.text.y = element_text(lineheight = 0.85, colour = "transparent"))

draw_panel <- function() {
  text_at("Intestinal studies with the gene 2-fold up in EC", .12, 4.525)
  draw_at(plot_F, .10, 4.685, 2.74, 1.86)
}
