# Figure S4D: ROC AUC, sensitivity and specificity of the 20-gene classifier and of a generic neuroendocrine gene
# set fitted the same way, in each scored study. Above the dotted line the studies whose primary study labels rank
# Fig. 2D's genes; the classifier of each study is trained on the labels of the other studies above the line.

gn <- both("generic_neuroendocrine.tsv", "hio_generic_neuroendocrine.tsv")
PK <- c("20_GENE_PANEL", "GENERIC_NEUROENDOCRINE")
NONE <- gn$panel == PK[2] & gn$ec_called == 0 & gn$other_called == 0                 # Singh 2024 and McCauley 2023
MET <- c("ROC AUC", "Sensitivity", "Specificity"); fm <- function(m) factor(m, levels = MET)
gl <- do.call(rbind, lapply(c("roc_auc", "sensitivity", "specificity"), function(k)
  data.frame(study = gn$study, panel = gn$panel, metric = c(roc_auc = "ROC AUC", sensitivity = "Sensitivity", specificity = "Specificity")[[k]], value = gn[[k]])))
gl$metric <- fm(gl$metric); gl$y <- unname(YR[gl$study]); gl$panel <- factor(gl$panel, levels = PK)
KL20 <- c(expression(paste("20-gene classifier", phantom(")"))),   # phantom: same box height as the second key, one baseline
          expression(paste("Generic neuroendocrine genes (", italic("CHGA, CHGB, SCG2, SCG3, PCSK1, PCSK2"), ")")))
none20 <- data.frame(metric = fm("Sensitivity"), y = unname(YR[gn$study[NONE]]) + .12)
dot20 <- data.frame(metric = fm("Specificity"), y = YDOT + c(.10, -.06), vj = c(0, 1),   # baseline above the line, cap top below
                    lab = c("Used to rank Fig. 2D's genes", "Not used to rank Fig. 2D's genes or to train"))
g20 <- ggplot(gl, aes(value, y)) +
  geom_hline(yintercept = YDOT, linetype = "dotted", color = "#777777", linewidth = .3) +
  geom_line(aes(group = y), color = "#BBBBBB", linewidth = .4) +
  geom_point(data = gl[gl$panel == PK[2], ], aes(color = panel, size = panel)) +    # generic set first, larger
  geom_point(data = gl[gl$panel == PK[1], ], aes(color = panel, size = panel)) +    # 20-gene classifier on top, smaller
  geom_text(data = data.frame(metric = fm("ROC AUC"), rowlab), aes(x = -Inf, y = y + .2, label = b), inherit.aes = FALSE,
            hjust = 1, fontface = "bold", size = tx, family = FONT) +
  geom_text(data = data.frame(metric = fm("ROC AUC"), rowlab), aes(x = -Inf, y = y - .2, label = p), inherit.aes = FALSE,
            hjust = 1, size = tx, family = FONT) +
  geom_text(data = none20, aes(x = .05, y = y, label = "no cells called"), inherit.aes = FALSE, hjust = 0, vjust = 0, size = tx, family = FONT) +
  geom_text(data = dot20, aes(x = Inf, y = y, label = lab, vjust = vj), inherit.aes = FALSE, hjust = 1, size = 7 / .pt, family = FONT) +
  facet_wrap(~metric, nrow = 1) +
  scale_color_manual(values = setNames(c(PANEL, "#AAAAAA"), PK), breaks = PK, labels = KL20, name = NULL) +
  scale_size_manual(values = setNames(c(1.8, 2.6), PK), breaks = PK, labels = KL20, name = NULL) +
  scale_x_continuous(limits = c(0, 1.02), breaks = c(0, .5, 1)) +
  scale_y_continuous(limits = c(.5, 5.5), breaks = NULL, expand = c(0, 0)) + coord_cartesian(clip = "off") +
  labs(x = NULL, y = NULL) + theme_kit("x", legend = "bottom") +
  theme(plot.margin = margin(2, 6, 2, 72 * LW + 2), legend.box.spacing = unit(5, "pt"), legend.margin = margin(0, 0, 0, 0))

draw_panel <- function() place(g20, 120, 3882, 4380, 1464, 7.3, 2.85 - (.30 + .17))
