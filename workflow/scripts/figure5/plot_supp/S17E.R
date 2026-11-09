# Figure S17E: CDX motif sites and CDX2 RNA in each study's EC population, both as EC minus the study's other endocrine
# cells. One row per study (bold study over plain sample type): CDX1/CDX2/HOXA10 motif-site accessibility in EC chromatin
# (chromVAR z; each study's z comes from its own peaks, so only signs compare), and CDX2 RNA (log2(1 + CPM)). McCauley 2023
# has no chromatin data. It shares its row with S17D, to its right.

TYPE_COL <- c("SC-SIC" = SC_SIC, "induced EC" = "#CCBB44", "native intestinal EC" = "#43A5B6", "HIO graft EC" = "#6943AC")
R <- data.frame(cohort = c("gse199636", "zhu", "singh", "hickey", "mccauley"),
                name = c("Augsornworawat 2023 (nuclei)", "Zhu 2023", "Singh 2024", "Hickey 2023 (nuclei)", "McCauley 2023"),
                sample_type = c("SC-SIC", "SC-SIC", "induced EC", "native intestinal EC", "HIO graft EC"))
cs <- read_table("cdx_sites/f5_s4_cdx_sites_vs_rna.tsv")
R$z <- cs$z[match(R$cohort, cs$cohort)]
RD <- read_table("cdx_sites/f5_s4_cdx2_rna_ec_minus_other.tsv"); stopifnot(setequal(RD$cohort, R$cohort))
R$rna_diff <- RD$ec_minus_other[match(R$cohort, RD$cohort)]
R$sample_type <- factor(R$sample_type, names(TYPE_COL)); R$y <- rev(seq_len(nrow(R)))
PITCH <- .255; YL <- c(.5, nrow(R) + .5); DY <- .11 / 2 / PITCH   # row pitch; study name and type .11 in apart
ZL <- c(-.15, 1.4); RL <- c(-.7, 2.8); MCY <- R$y[R$cohort == "mccauley"]
LAB_W <- str_in(R$name, PT, TRUE)                   # the label column fits the longest bold study name
lab <- ggplot() +
  annotate("text", x = 1, y = R$y + DY, label = R$name, hjust = 1, vjust = .5, fontface = "bold", size = TXT, family = FONT, colour = "black") +
  annotate("text", x = 1, y = R$y - DY, label = as.character(R$sample_type), hjust = 1, vjust = .5, size = TXT, family = FONT, colour = "black") +
  scale_x_continuous(limits = c(0, 1), expand = expansion(0)) + scale_y_continuous(limits = YL, expand = expansion(0)) +
  coord_cartesian(clip = "off") + theme_void() + theme(plot.margin = margin(2, 0, 2, 2.9), plot.background = element_blank())
th <- theme_panel("both") +
  theme(axis.text.y = element_blank(), axis.ticks.y = element_blank(), axis.title.y = element_blank(),
        axis.title.x = element_text(colour = "black", size = PT, lineheight = 1, hjust = .5, margin = margin(3, 0, 0, 0)),
        plot.title = element_text(colour = "black", size = PT, face = "bold", hjust = 0, lineheight = 1, margin = margin(0, 0, 3, 0)),
        plot.title.position = "panel")
pz <- ggplot(R[!is.na(R$z), ], aes(z, y, colour = sample_type)) +
  annotate("rect", xmin = -Inf, xmax = Inf, ymin = MCY - .5, ymax = MCY + .5, fill = "white", colour = NA) +
  annotate("segment", x = 0, xend = 0, y = MCY + .5, yend = YL[2], colour = REF_COL, linewidth = .25) +
  annotate("text", x = mean(ZL), y = MCY, label = "no chromatin data", size = TXT, family = FONT, colour = "black") +
  geom_point(aes(fill = sample_type), size = 2.2, shape = 21, colour = "black", stroke = STROKE_THIN) +
  scale_fill_manual(values = TYPE_COL, guide = "none") +
  scale_x_continuous(limits = ZL, breaks = c(0, .5, 1), labels = c("0", "0.5", "1"), expand = expansion(0)) +
  scale_y_continuous(limits = YL, breaks = R$y, expand = expansion(0)) +
  labs(x = "Accessibility of CDX motif sites,\nEC minus other endocrine cells\n(chromVAR z; compare signs only)") + th + theme(plot.margin = margin(2, 16, 2, 8))
pr <- ggplot(R, aes(rna_diff, y)) +
  geom_vline(xintercept = 0, colour = REF_COL, linewidth = .25) +
  geom_point(aes(fill = sample_type), size = 2.2, shape = 21, colour = "black", stroke = STROKE_THIN) +
  scale_fill_manual(values = TYPE_COL, guide = "none") +
  scale_x_continuous(limits = RL, breaks = c(0, 1, 2), expand = expansion(0)) +
  scale_y_continuous(limits = YL, breaks = R$y, expand = expansion(0)) +
  labs(x = "CDX2 RNA,\nEC minus other endocrine cells\n(log2(1 + CPM))") + th + theme(plot.margin = margin(2, 9, 2, 4))
pE <- lab + pz + pr + plot_layout(widths = unit(c(LAB_W, 1.25, 1), c("in", "null", "null")), heights = unit(nrow(R) * PITCH, "in")) +
  plot_annotation(theme = theme(plot.margin = margin(0, 0, 0, 0), plot.background = element_blank()))

draw_panel <- function() in_row("DE", 2, function(w, vh) in_vp(0, 0, w - .05, vh, function() print(pE, newpage = FALSE)))
