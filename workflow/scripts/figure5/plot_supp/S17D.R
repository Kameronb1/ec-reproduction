# Figure S17D: Figure 5D's graft comparison on promoters instead of distal elements. In each graft, within SC-SIC, SC-beta
# and SC-alpha, the change of each gene group's promoter accessibility from week 2 in vitro beyond the typical gene's
# (points, as %), against the middle 95% of 10,000 level-matched random gene sets (grey bars). Group colors in the S17A-C
# key; the bar's key below. It shares its row with S17E.

pct <- function(d) 100 * (2^d - 1)
gr <- read_table("c5_graft_relative_by_population.tsv"); gr <- gr[gr$window == "promoter", ]
gr$y <- pct(gr$group_rel); gr$lo <- pct(gr$null_q025_rel); gr$hi <- pct(gr$null_q975_rel)
gr$u <- factor(c(txp_1 = "1", txp_2 = "2", txp_3 = "3")[gr$graft], c("1", "2", "3")); gr$win <- factor(gr$population, c("SC-SIC", "SC-beta", "SC-alpha"))
pD <- ggplot(gr, aes(u, y, group = g)) + geom_hline(yintercept = 0, colour = REF_COL, linewidth = .3) +
  geom_linerange(aes(ymin = lo, ymax = hi), position = position_dodge(.8), linewidth = .85, colour = BANDC) +
  geom_point(aes(fill = g), shape = 21, colour = "black", stroke = STROKE_THIN, size = 1.2, position = position_dodge(.8)) + facet_wrap(~win, nrow = 1) +
  scale_fill_manual(values = GCOL) + labs(x = "Graft (week 2 → graft)", y = "% more opening than the
typical gene's promoter") +
  theme_panel("y") + theme(strip.text = element_text(face = "bold", hjust = .5, size = PT), plot.margin = margin(0, 2, 0, 2), panel.spacing = unit(5, "pt"))

draw_panel <- function() in_row("DE", 1, function(w, vh) {
  in_vp(0, .2, w, vh - .2, function() print(pD, newpage = FALSE))
  ky2(list(list(type = "bar", label = "Middle 95% of 10,000 random gene sets")), w / 2 + .06, .08)
})
