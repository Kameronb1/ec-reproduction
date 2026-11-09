# Figure S18B: outside checks of the Figure 5E factors in induced EC (Singh 2024) and HIO graft EC (McCauley 2023, the
# primary study label Endocrine_1 against Endocrine_2): EC minus other enteroendocrine regulon target score, one dot per
# regulon, colored as in Figure 5E; beside each study, how many of a factor's reproducible regulons keep the
# intestinal direction. Neither study is used for the Figure 5E colors.

FA <- factor_rows()
BANDS <- row_bands(FA)
singh <- read_table("identity/f5_s_gh_singh_values_ALL_ELIGIBLE.tsv"); singh <- singh[singh$part == "singh_effect", ]
mcc <- read_table("identity/gse214852_values_primary_labels.tsv")
P <- rbind(data.frame(study = "McCauley 2023 HIO graft EC", factor = mcc$factor, program = mcc$program, effect = mcc$gse214852_effect,
                      col = ifelse(mcc$program_pass, ifelse(mcc$direction == "SC-SIC higher", "pass_ec", "pass_beta"), "fail")),
           data.frame(study = "Singh 2024 induced EC", factor = singh$factor, program = singh$program, effect = singh$value, col = singh$note))
P <- P[order(P$study, P$factor, P$program, method = "radix"), ]   # regulon order sets each regulon's place inside its row
STUDIES <- c("Singh 2024 induced EC", "McCauley 2023 HIO graft EC")
HEAD <- list(c("Singh 2024", "induced EC"), c("McCauley 2023", "HIO graft EC")); names(HEAD) <- STUDIES
P$d <- FA$d[match(P$factor, FA$factor)]; P$half <- FA$half[match(P$factor, FA$factor)]
P$y <- FA$y[match(P$factor, FA$factor)] + ave(seq_len(nrow(P)), P$study, P$factor, FUN = spread)
# per study and factor: reproducible regulons whose sign matches the factor's direction, of all reproducible regulons
pp <- P[P$col != "fail", ]
CN <- aggregate(cbind(n = 1, k = as.integer(sign(pp$effect) == pp$d)) ~ study + factor, pp, sum)
CN$lab <- sprintf("%d of %d", CN$k, CN$n); CN$y <- FA$y[match(CN$factor, FA$factor)]; CN$half <- FA$half[match(CN$factor, FA$factor)]
SLIM <- tapply(abs(P$effect), P$study, max) * 1.05
HB2 <- do.call(rbind, lapply(HALF, function(h) { yy <- FA$y[FA$half == h]; data.frame(half = factor(h, HALF), y = c(min(yy) - .5, max(yy) + .5)) }))
CHEAD <- c("Colored regulons with", "the intestinal direction
in this study's EC")
XLAB <- "EC minus other enteroendocrine (regulon target score); one dot per regulon"
CNT_W <- 1.08; GAP_IN <- .22                          # count column width; gap before the second study (in)
thB <- theme_panel("x", legend = "bottom") +
  theme(strip.placement = "outside", strip.text.y.left = element_text(angle = 0, hjust = 1, face = "bold", size = PT, lineheight = .95),
        strip.text.y.right = element_blank(), panel.spacing.y = unit(.16, "in"),
        plot.title = element_text(size = 8, face = "bold", hjust = 0, margin = margin(0, 0, 1.5, 0)),
        plot.subtitle = element_text(size = PT, hjust = 0, vjust = 1, lineheight = .9, margin = margin(0, 0, 4, 0)),
        axis.title = element_blank(), plot.margin = margin(0, 0, 0, 0))
dots <- function(st, first) {
  L <- SLIM[[st]]
  g <- ggplot() +
    geom_rect(data = BANDS, aes(xmin = -Inf, xmax = Inf, ymin = ymin, ymax = ymax), fill = "#D9D9D9", colour = NA) +
    geom_vline(xintercept = 0, colour = REF_COL, linewidth = .3) + geom_blank(data = HB2, aes(x = 0, y = y)) +
    geom_point(data = P[P$study == st, ], aes(effect, y, colour = col), size = 1) +
    facet_grid(half ~ ., scales = "free_y", space = "free_y", switch = "y") +
    scale_y_continuous(breaks = FA$y, labels = parse(text = FA$expr), expand = expansion(add = c(.12, 0))) +
    scale_x_continuous(limits = c(-L, L), breaks = scales::breaks_pretty(3), expand = c(0, 0)) + coord_cartesian(clip = "off") +
    scale_colour_manual(values = REG_COL, labels = REG_LAB, breaks = names(REG_LAB), limits = names(REG_LAB), drop = FALSE, name = REG_KEY_TITLE) +
    guides(colour = guide_legend(override.aes = list(size = 1.6))) +
    labs(title = HEAD[[st]][1], subtitle = HEAD[[st]][2]) + thB      # column header: bold study over plain sample type
  if (!first) g <- g + theme(axis.text.y = element_blank(), strip.text.y.left = element_blank(), plot.margin = margin(0, 0, 0, GAP_IN * 72))
  g
}
counts <- function(st) {
  ggplot() +
    geom_rect(data = BANDS, aes(xmin = -Inf, xmax = Inf, ymin = ymin, ymax = ymax), fill = "#D9D9D9", colour = NA) +
    geom_blank(data = HB2, aes(x = .5, y = y)) +
    geom_text(data = CN[CN$study == st, ], aes(x = .5, y = y, label = lab), size = PT / .pt, family = FONT) +
    facet_grid(half ~ ., scales = "free_y", space = "free_y", switch = "y") +
    scale_y_continuous(breaks = FA$y, expand = expansion(add = c(.12, 0))) + scale_x_continuous(limits = c(0, 1), expand = c(0, 0)) +
    labs(title = CHEAD[1], subtitle = CHEAD[2]) + thB +              # column header, two plain lines
    theme(axis.text = element_blank(), axis.line.x = element_blank(), axis.ticks.x = element_blank(), panel.grid.major.x = element_blank(),
          strip.text.y.left = element_blank(), plot.title = element_text(size = PT, face = "plain", hjust = .5, margin = margin(0, 0, 1.5, 0)),
          plot.subtitle = element_text(size = PT, hjust = .5, vjust = 1, margin = margin(0, 0, 4, 0)))
}
xt <- wrap_elements(panel = textGrob(XLAB, x = .5, y = 1, just = c("centre", "top"), gp = gpar(fontfamily = FONT, fontsize = PT, col = "black"))) +
  theme(plot.margin = margin(3, 0, 0, 0), plot.background = element_blank())
pB <- dots(STUDIES[1], TRUE) + counts(STUDIES[1]) + dots(STUDIES[2], FALSE) + counts(STUDIES[2]) + xt +
  plot_layout(design = "ABCD\nEEEE", guides = "collect",
              widths = unit(c(1, CNT_W, 1, CNT_W), c("null", "in", "null", "in")),
              heights = unit(c(1, .16), c("null", "in"))) &
  theme(legend.position = "bottom", legend.title.position = "left", legend.margin = margin(0, 0, 0, 0), legend.box.spacing = unit(0, "pt"),
        legend.justification = "left",
        legend.title = element_text(size = PT, colour = "black", vjust = .5), legend.text = element_text(size = PT, colour = "black"))
pB <- pB + plot_annotation(theme = theme(plot.margin = margin(2, 8, 6, 2), legend.position = "bottom", legend.justification = "left",
                                         legend.box.spacing = unit(4, "pt"), plot.background = element_blank()))

draw_panel <- function() place(pB, 2873, 2478)
