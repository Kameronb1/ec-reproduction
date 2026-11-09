# Figure S4C: ROC AUC (EC versus other terminal EEC) of the 20-gene classifier in each scored study against 200
# random 20-gene sets matched for expression, detection and variance, each weighted and scored the same way (random
# sets jittered vertically within their row). Rows and dotted line as in S4D. Dots have thin black outlines;
# outlined shapes (21, 23) draw smaller than solid ones at the same size value, so their sizes are set larger.

mr <- both("matched_random.tsv", "hio_matched_random.tsv")
mr$y <- unname(YR[mr$study])
real <- mr[mr$panel == "REAL", ]; rnd <- mr[mr$panel == "MATCHED_RANDOM", ]
XL19 <- c(.2, 1.2); XV19 <- 1.065                                     # x range; x of the printed AUCs
K19 <- c("20-gene classifier", sprintf("Random 20-gene sets matched for expression,\ndetection and variance (%d per study)", max(rnd$replicate)))
plot_c <- function(lo) ggplot(rnd, aes(roc_auc, y)) +                  # lo: row labels' offset from the row center
  geom_hline(yintercept = YDOT, linetype = "dotted", color = "#777777", linewidth = .3) +
  geom_point(aes(fill = K19[2], shape = K19[2], size = K19[2]), position = position_jitter(height = .18, width = 0, seed = 19),
             colour = "black", stroke = .1) +
  geom_point(data = real, aes(fill = K19[1], shape = K19[1], size = K19[1]), colour = "black", stroke = .1) +
  scale_fill_manual(values = setNames(c(PANEL, "#AAAAAA"), K19), breaks = K19, name = NULL) +
  scale_shape_manual(values = setNames(c(23, 21), K19), breaks = K19, name = NULL) +
  scale_size_manual(values = setNames(c(2.3, .6), K19), breaks = K19, name = NULL) +
  guides(fill = guide_legend(override.aes = list(size = c(2.3, 1.7)))) +
  geom_text(data = real, aes(x = XV19, label = sprintf("%.2f", roc_auc)), hjust = 0, size = tx, family = FONT) +
  geom_text(data = rowlab, aes(x = -Inf, y = y + lo, label = b), inherit.aes = FALSE, hjust = 1, fontface = "bold", size = tx, family = FONT) +
  geom_text(data = rowlab, aes(x = -Inf, y = y - lo, label = p), inherit.aes = FALSE, hjust = 1, size = tx, family = FONT) +
  scale_x_continuous(limits = XL19, breaks = seq(.2, 1, .2), expand = c(0, 0)) +
  scale_y_continuous(limits = c(.45, 5.55), breaks = NULL, expand = c(0, 0)) + coord_cartesian(clip = "off") +
  labs(x = "ROC AUC, EC versus\nother terminal EEC", y = NULL) + theme_kit("x", legend = "bottom") +
  theme(legend.direction = "vertical", legend.key.spacing.y = unit(1, "pt"), legend.justification = "right",
        plot.margin = margin(2, 6, 2, 72 * LW + 2))
# The panel reaches down C_TALL pixels further, to the foot of S4A; the rows get the height and the two-line row
# labels keep their spacing in inches: 0.2 rows at the plotting height of the panel without the extra pixels.
C_TALL <- 246
ph0 <- panel_h(plot_c(.2), 3.05 - .47)
g19 <- plot_c(.2 * ph0 / (ph0 + C_TALL / 600))

draw_panel <- function() place(g19, 2742, 1992, 1758, 1830, 2.93, 3.05 - (.30 + .17) + C_TALL / 600)
