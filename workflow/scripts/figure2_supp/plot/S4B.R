# Figure S4B: genes in all seven leave-one-dataset-out top-20 lists, redone with EC and other labels swapped: within
# whole gene-ranking datasets (all 2^7 = 128 combinations, the primary study labels among them) and within donors (200 random
# draws). Histograms of the label sets, with the count for the primary study labels (downward arrow) and its p-value.

nul <- read_table("recurrence_null.tsv")                       # size-20 lists of every completed replicate
sm <- read_table("recurrence_summary.tsv")
obs <- read_table("observed_recurrence.tsv")
o20 <- obs$genes_in_every_fold[obs$panel_size == 20]
lv <- c(FAMILY_BLOCK = sprintf("Swapped within whole datasets\n(all %d combinations,
primary study labels included)", sum(nul$swap_level == "FAMILY_BLOCK")),
        UNIT = sprintf("Swapped within donors\n(%d random draws)", sum(nul$swap_level == "UNIT")))
nul$lev <- factor(lv[nul$swap_level], levels = lv)
s20 <- sm[sm$panel_size == 20, ]; s20$lev <- factor(lv[s20$swap_level], levels = lv)
s20$lab <- sprintf("observed %d; p = %s", s20$observed, formatC(s20$empirical_p, format = "g", digits = 2))
# the observed count as a downward arrow (as in Figure S15C), stopping at a third of each facet's tallest bar
key_down <- function(data, params, size)                       # key glyph: the drawn downward arrow
  grid::segmentsGrob(.5, .95, .5, .08, arrow = params$arrow,
                     gp = grid::gpar(col = data$colour, fill = data$colour, lwd = data$linewidth * .pt, lineend = "butt", linejoin = "mitre"))
arr <- data.frame(lev = factor(lv, levels = lv), k = "Primary study labels",
                  yend = .32 * vapply(lv, function(l) max(table(nul$genes_in_every_fold[nul$lev == l])), 0))
g13 <- ggplot(nul, aes(genes_in_every_fold)) +
  geom_bar(aes(fill = "Label sets"), width = .8) +
  geom_segment(data = arr, aes(x = o20, xend = o20, y = Inf, yend = yend, color = k), inherit.aes = FALSE, linewidth = .7,
               arrow = arrow(length = unit(.055, "in"), angle = 28, type = "closed"), key_glyph = key_down) +
  scale_fill_manual(values = c("Label sets" = "#AAAAAA"), name = NULL) + scale_color_manual(values = c("Primary study labels" = PANEL), name = NULL) +
  guides(color = guide_legend(order = 1), fill = guide_legend(order = 2)) +
  geom_text(data = s20, aes(x = o20 - .3, y = Inf, label = lab), inherit.aes = FALSE, hjust = 1, vjust = 1.4, size = tx, family = FONT) +
  facet_wrap(~lev, ncol = 1, scales = "free_y") +
  scale_x_continuous(breaks = seq(0, 16, 4), limits = c(-.6, 16)) +
  labs(x = "Genes in all seven top-20 lists", y = "Label sets") + theme_kit("y", legend = "bottom") +
  theme(strip.text = element_text(lineheight = .95))

draw_panel <- function() place(g13, 2742, 48, 1758, 1884, 2.93, 3.55 - (.30 + .17))
