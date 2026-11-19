# Figure 3A: per sample, the median over each gene group's genes of log2(1 + CPM), for every population; the bar is the
# median of the samples. Populations used to define the groups are shaded.
M <- as.data.table(read_table("group_medians.tsv"))
ORD <- names(POP_COL)                                      # rows, top to bottom
stopifnot(setequal(unique(M$pop), ORD))
N_INT <- sum(ORD %in% c("Primary intestinal EC", "Other primary EEC", "HIO graft EC"))
NS <- M[g == "G1", .N, by = pop]; ROWLAB <- setNames(sprintf("%s (%d)", ORD, NS$N[match(ORD, NS$pop)]), ORD)
NR <- length(ORD)
M[, `:=`(y = NR + 1 - match(pop, ORD), facet = factor(GLAB[g], GLAB))]
set.seed(1); M[, yj := y + runif(.N, -.2, .2)]
setorder(M, nuclei)                                        # nuclei (open) drawn on top
MED <- M[, .(v = median(median_level)), by = .(y, facet)]
DEFINING <- unique(M[used_to_define == TRUE]$pop)
BAND <- CJ(y = NR + 1 - match(DEFINING, ORD), facet = factor(GLAB, GLAB))
CHIP <- data.table(facet = factor(GLAB, GLAB), fill = GROUP_COL)
SEP <- NR + 1 - c(N_INT, N_INT + 4) - .5                   # under the intestinal rows and under the SC-islet rows
XMAX <- 8.6
pa <- ggplot(M, aes(median_level, yj)) +
  geom_rect(data = BAND, aes(xmin = -Inf, xmax = Inf, ymin = y - .5, ymax = y + .5), inherit.aes = FALSE, fill = "#EDEDED") +
  geom_rect(data = CHIP, aes(xmin = -Inf, xmax = Inf, ymin = NR + .62, ymax = NR + .86, fill = I(fill)), inherit.aes = FALSE) +
  geom_vline(xintercept = seq(0, 8, 2), colour = GRID_COL, linewidth = .2) +
  geom_hline(yintercept = SEP, colour = "grey70", linewidth = .3) +
  geom_errorbar(data = MED, aes(xmin = v, xmax = v, y = y), width = .78, linewidth = .6, colour = "grey15", orientation = "y", inherit.aes = FALSE) +
  geom_point(aes(colour = pop, shape = nuclei), size = 1.3, stroke = .55, fill = "white") +
  facet_wrap(~facet, nrow = 1) +
  scale_colour_manual(values = POP_COL) + scale_shape_manual(values = c(`FALSE` = 16, `TRUE` = 21)) +
  scale_y_continuous(breaks = seq_len(NR), labels = rev(ROWLAB[ORD]), limits = c(.5, NR + .86), expand = c(0, 0)) +
  scale_x_continuous(limits = c(0, XMAX), breaks = seq(0, 8, 2), expand = expansion(mult = c(.02, .01))) +
  labs(x = expression("Median over the genes of the group, log"[2] * "(1 + CPM); one point per sample"), y = NULL) +
  panel_theme() + theme(panel.grid.major.x = element_blank(), panel.spacing.x = unit(.16, "in"), strip.text = element_text(margin = margin(1, 0, 3, 0)))
A_W <- 7.3; A_PLOT_H <- 3.0; A_H <- A_PLOT_H + .34
draw_A <- function() {
  .VH <<- A_H
  pushViewport(viewport(x = 0, y = unit(A_H - A_PLOT_H, "in"), width = unit(A_W, "in"), height = unit(A_PLOT_H, "in"), just = c("left", "bottom")))
  gt <- ggplotGrob(pa)
  s <- gt$layout[gt$layout$name == "strip-t-1-1", ]; a <- gt$layout[gt$layout$name == "axis-l-1-1", ]
  hdr <- textGrob("Population (samples)", x = unit(1, "npc") - unit(2.8, "pt"), y = unit(3, "pt"), just = c("right", "bottom"), gp = gp(PT, "bold"))
  grid.draw(gtable::gtable_add_grob(gt, hdr, t = s$t, l = a$l, clip = "off", name = "row-header"))
  popViewport()
  key_row(list(sq_item("Used to define the groups", "#EDEDED", "grey60", .4), sq_item("Not used to define the groups", "white", "grey60", .4),
               dot_item("Whole cells"), dot_item("Single nuclei", open = TRUE), tick_item("Median of samples")), 0, A_W, A_PLOT_H + .14)
}
draw_panel <- function() panel_canvas(120, 48, A_W, A_H, draw_A)
