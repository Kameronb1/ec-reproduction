# Figure S6A: share of SC-SIC (left) and of other SC endocrine cells (right) classified as EC in each test SC-islet study
# (rows) by the four-marker classifier, the TPH1-only classifier and the label-trained 20-gene classifier of Figure 2E
# with its gene list (Figure 2D's intestinal genes or the 20 SC-SIC genes) and its training set (intestinal or SC-islet
# studies) swapped. Each bar is labeled with its percent and, in brackets, the further percent of cells above the cutoff
# that the neighborhood check leaves unclassified (open bar). Key below: the four-marker and TPH1-only classifiers, then the 20-gene
# classifier as gene list x training.

calls <- read_table("s2_classifier_swap/calls_by_study.tsv")
RULES <- c("MARKER_RULE_GATED", "TPH1_GATED", "GUTG_GUTT", "GUTG_SCT", "SCG_GUTT", "SCG_SCT")
calls <- calls[calls$combination %in% RULES, ]                # the six drawn classifiers (four-marker: with its neighborhood check)
# hue = gene list (teal: Figure 2D intestinal genes; wine: the 20 SC-SIC genes); shade = training set (intestinal or SC-islet
# studies); the four-marker classifier rust, the TPH1-only classifier grey
COLS <- c(MARKER_RULE_GATED = "#B86E46", TPH1_GATED = "#ABABAB", GUTG_GUTT = "#2A7F8E", GUTG_SCT = "#0D3238", SCG_GUTT = "#E3AECB", SCG_SCT = "#882255")   # rust; grey; intestinal-gene classifiers teal (= Fig. S2); SC-SIC-gene classifiers wine (= Fig. S15)
nn <- aggregate(n ~ held_out_study + class, calls[calls$combination == "MARKER_RULE_GATED", ], max)
nlab <- function(f, cl) format(nn$n[nn$held_out_study == f & nn$class == cl], big.mark = ",")
labels_y <- vapply(names(STUDY), function(f) paste0(sub(" (nuclei)", "\n(nuclei)", STUDY[[f]], fixed = TRUE), "\n",
  nlab(f, "SC_SIC"), " | ", nlab(f, "OTHER_MATURE_SC_ENDOCRINE"), " cells"), "")
calls$row <- NS + 1 - match(calls$held_out_study, names(STUDY))
off <- setNames(seq(.33, -.33, length.out = length(RULES)), RULES)   # six bars per study row
calls$y <- calls$row + off[calls$combination]
FAC <- c(SC_SIC = "SC-SIC", OTHER = "Other SC endocrine cells")
calls$facet <- factor(ifelse(calls$class == "SC_SIC", FAC[["SC_SIC"]], FAC[["OTHER"]]), FAC)
calls$cand <- 100 * calls$stage_one_called / calls$n
low <- round(calls$cand - 100 * calls$called / calls$n)
calls$lab <- paste0(round(100 * calls$called / calls$n), "%", ifelse(low >= 1, paste0("  (+", low, ")"), ""))
calls$xmax <- ifelse(calls$class == "SC_SIC", 100, 40)          # each facet's own scale end
calls$pct <- 100 * calls$called / calls$n
calls$combination <- factor(calls$combination, levels = RULES)
XMAX <- 130; BH <- .118; COUNT_PT <- 7.2
plot_A <- ggplot(calls) +
  geom_rect(aes(xmin = pct, xmax = cand, ymin = y - BH / 2 + .008, ymax = y + BH / 2 - .008, colour = combination), fill = "white", linewidth = .3) +
  geom_rect(aes(xmin = 0, xmax = pct, ymin = y - BH / 2, ymax = y + BH / 2, fill = combination)) +
  geom_blank(aes(x = xmax * 1.02)) +
  geom_hline(yintercept = seq(1.5, NS - .5, 1), colour = "#BBBBBB", linewidth = .25, linetype = "dotted") +
  facet_wrap(~facet, ncol = 2, scales = "free_x") +
  scale_fill_manual(values = COLS, guide = "none") + scale_colour_manual(values = COLS, guide = "none") +
  scale_y_continuous(breaks = NS:1, labels = unname(labels_y), limits = c(.5, NS + .5), expand = c(0, 0)) +
  scale_x_continuous(breaks = function(l) if (l[2] > 60) c(0, 25, 50, 75, 100) else c(0, 10, 20, 30), labels = function(v) paste0(v, "%"), expand = c(0, 0)) +
  coord_cartesian(clip = "off") +
  labs(x = "Cells called (% of cells with that primary study label). Note the different scales", y = NULL) +
  theme_rk("x") +
  theme(axis.text.y = element_text(lineheight = .9), strip.text = element_text(size = 8, face = "bold", hjust = 0),
        panel.spacing.x = unit(.30, "in"), plot.margin = margin(2, 14, 2, 2))

# Key band under the plot (inches): three lines. 1: four-marker, TPH1-only, outline; 2-3: the 20-gene classifier's
# gene lists, each with its intestinal- and SC-trained shades.
KEY_H <- .64; SW <- c(.17, .11)                               # key band height; swatch width, height
key_draw <- function(x0) {
  tp <- gpar(fontfamily = FONT, fontsize = PT, col = "black")
  sw <- function(x, y, col, open = FALSE) grid.rect(unit(x, "in"), unit(y, "in"), unit(SW[1], "in"), unit(SW[2], "in"), just = c("left", "centre"),
    gp = if (open) gpar(fill = "white", col = "#4D4D4D", lwd = .8) else gpar(fill = col, col = NA))
  tx <- function(s, x, y) { grid.text(s, unit(x, "in"), unit(y, "in"), just = c("left", "centre"), gp = tp); x + str_w(s) }
  item <- function(x, y, col, s, open = FALSE) { sw(x, y, col, open); tx(s, x + SW[1] + .07, y) + .22 }
  y <- KEY_H - .10 - (0:2) * .19
  x <- item(x0, y[1], COLS[["MARKER_RULE_GATED"]], "Four-marker classifier (Fig. 2A)")
  x <- item(x, y[1], COLS[["TPH1_GATED"]], "TPH1-only classifier")
  item(x, y[1], NA, "Above the cutoff, not called (neighborhood check)", open = TRUE)
  lab <- c("20-gene classifier, Fig. 2D intestinal genes:", "20-gene classifier, SC-SIC genes:"); cx <- x0 + max(str_w(lab)) + .12
  for (r in 1:2) { tx(lab[r], x0, y[r + 1])
    x <- item(cx, y[r + 1], COLS[[c("GUTG_GUTT", "SCG_GUTT")[r]]], if (r == 1) "intestinal training (= Fig. 2E)" else "intestinal training")
    item(x, y[r + 1], COLS[[c("GUTG_SCT", "SCG_SCT")[r]]], "SC-islet training") }
}
draw_A <- function() {
  gt <- ggplotGrob(plot_A)                                        # built on the open device, which measures the text
  vp_h <- convertHeight(unit(1, "npc"), "in", valueOnly = TRUE)
  pushViewport(viewport(y = unit(KEY_H, "in"), height = unit(vp_h - KEY_H, "in"), just = "bottom"))
  grid.draw(gt)
  pl <- gt$layout[grepl("^panel", gt$layout$name), ]; l1 <- min(pl$l)   # the key starts at the first facet's left edge
  x0 <- sum(vapply(seq_len(l1 - 1), function(i) convertWidth(gt$widths[i], "in", valueOnly = TRUE), 0))
  popViewport()
  pushViewport(viewport(y = unit(0, "in"), height = unit(KEY_H, "in"), just = "bottom")); key_draw(x0); popViewport()
}

draw_panel <- function() in_place(120, 48, 7.3, 7.20 + .17 - .47, 4140, draw_A)
