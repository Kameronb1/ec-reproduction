# Figure S15D: ROC AUC of the 20 SC-SIC genes (SC-SIC against other SC endocrine cells, gene weights fitted on the other
# five studies) in each left-out SC-islet study (rows), beside 200 random sets of 20 genes matched for expression
# (grey dots, jittered vertically); the 20 genes' AUC is also written at the right. Dashed line: chance. Key below.

auc <- read_table("s1_marker_stability/c_matched_random.tsv.gz")
auc$st <- factor(STUDY[auc$study], levels = rev(unname(STUDY)))
real <- auc[auc$panel == "REAL", ]; rnd <- auc[auc$panel == "MATCHED_RANDOM", ]
K <- c("The 20 SC-SIC genes", sprintf("Random sets of 20 expression-matched genes (%d per study)", unique(as.vector(table(rnd$study)))))

plot_D <- ggplot(rnd, aes(roc_auc, st)) +
  geom_vline(data = data.frame(k = "Chance (0.5)"), aes(xintercept = .5, linetype = k), colour = REF, linewidth = .4) +
  geom_point(aes(color = K[2], shape = K[2], size = K[2]), position = position_jitter(height = .18, width = 0, seed = 19)) +
  geom_point(data = real, aes(color = K[1], shape = K[1], size = K[1])) +
  scale_color_manual(values = setNames(c(SCC, "#AAAAAA"), K), breaks = K, name = NULL) +
  scale_shape_manual(values = setNames(c(18, 16), K), breaks = K, name = NULL) +
  scale_size_manual(values = setNames(c(2.6, .5), K), breaks = K, name = NULL) +
  scale_linetype_manual(values = c("Chance (0.5)" = "dashed"), name = NULL) +
  guides(color = guide_legend(order = 1, override.aes = list(size = c(2.6, 1.4))), shape = guide_legend(order = 1), size = guide_legend(order = 1),
         linetype = guide_legend(order = 2, override.aes = list(linewidth = .4))) +
  geom_text(data = real, aes(x = 1.02, label = sprintf("%.2f", roc_auc)), hjust = 0, size = PT / .pt, family = FONT, colour = "black") +
  scale_x_continuous(limits = c(min(.3, floor(10 * min(auc$roc_auc, na.rm = TRUE)) / 10), 1.10), breaks = seq(0, 1, .2), labels = minus, expand = c(0, 0)) +
  labs(x = "ROC AUC in the left-out study", y = NULL) + theme_rk("x", legend = "bottom") +
  theme(legend.direction = "vertical", legend.box = "vertical", legend.box.just = "left", legend.location = "plot", legend.justification = "left",
        legend.spacing.y = unit(0, "in"), legend.margin = margin(0, 0, 0, 0), legend.key.width = unit(.16, "in"), legend.key.height = unit(.15, "in"))

draw_panel <- function() in_place(2010, 3502, 4.10, 3.35 - .47, 1728, function() print(plot_D, newpage = FALSE))
