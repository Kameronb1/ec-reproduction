# Figure 3C: for each published gene list of serotonin-producing SC-islet cells, the share of each gene group's genes in
# the top of the list (filled), with the best other program or cluster of the same study as a tick.
T9 <- as.data.table(read_table("preprint_gene_lists.tsv"))
T9[, g := gene_set]
SRC <- c("Liu 2025, SC-SIC program", "Fiancette 2026, H1 line", "Fiancette 2026, H3 line")
T9[, src := factor(SRC[match(paste(preprint, cell_line), c("Liu 2025 (preprint) H1 and HUES8 knockout village", "Fiancette 2026 (preprint) H1",
                                                           "Fiancette 2026 (preprint) H3"))], SRC)]
stopifnot(!anyNA(T9$src), nrow(T9) == 9L)
GSIZE <- c(G1 = 28L, G2 = 34L, G3 = 84L)
if (!all(T9$genes_in_set == GSIZE[T9$g])) stop("group sizes are not 28, 34, 84")
T9[, `:=`(n = genes_in_set, top = genes_meeting_rule, other = genes_meeting_rule_there)]
T9[, `:=`(pct = 100 * top / n, opct = 100 * other / n, lab = sprintf("%d of %d", top, n), gy = factor(g, c("G3", "G2", "G1")))]
T9[, yn := as.numeric(gy)]
pc <- ggplot(T9) +
  geom_rect(aes(xmin = 0, xmax = 100, ymin = yn - .34, ymax = yn + .34, colour = I(GROUP_COL[g])), fill = "white", linewidth = .35) +
  geom_rect(aes(xmin = 0, xmax = pct, ymin = yn - .34, ymax = yn + .34, fill = I(GROUP_COL[g])), colour = NA) +
  geom_segment(aes(x = opct, xend = opct, y = yn - .44, yend = yn + .44), colour = "grey15", linewidth = .6, lineend = "butt") +
  geom_text(aes(x = 103, y = yn, label = lab), hjust = 0, size = PT / .pt, family = FONT) +
  facet_wrap(~src, nrow = 1) +
  scale_y_continuous(breaks = 1:3, labels = c("G3", "G2", "G1"), expand = expansion(add = .1)) +
  scale_x_continuous(breaks = seq(0, 100, 25), expand = expansion(add = c(1, 1))) +
  coord_cartesian(xlim = c(0, 100), clip = "off") +                # the count labels sit right of the 100% end
  labs(x = "Group genes in the top of the published list (%)", y = NULL) + panel_theme() +
  theme(panel.grid.major.x = element_line(colour = GRID_COL, linewidth = .2), panel.spacing.x = unit(.50, "in"),
        strip.text = element_text(margin = margin(2, 0, 3, 0)), axis.text.y = element_text(colour = "black", size = PT, margin = margin(0, 2, 0, 0)),
        plot.margin = margin(2, 36, 2, 2))
C_W <- 7.3; C_PLOT_H <- 1.30; C_KEY_H <- .42; C_H <- C_PLOT_H + C_KEY_H
draw_C <- function() {
  .VH <<- C_H
  pushViewport(viewport(x = 0, y = unit(C_H - C_PLOT_H, "in"), width = unit(C_W, "in"), height = unit(C_PLOT_H, "in"), just = c("left", "bottom")))
  print(pc, newpage = FALSE); popViewport()
  ky <- C_PLOT_H + .10
  key_row(list(sq_item("All genes of the group", "white", "grey45", .9), tick_item("Best other program or cluster of the same study")), 0, C_W, ky)
  key_row(list(sq_item("In the top of the list (Liu 2025: z score > 2; Fiancette 2026: first 500 markers)", "grey45")), 0, C_W, ky + .16)
}
draw_panel <- function() panel_canvas(120, 3969, C_W, C_H, draw_C)
