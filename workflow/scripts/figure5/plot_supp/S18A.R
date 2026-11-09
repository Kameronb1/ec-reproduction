# Figure S18A: the intestinal evidence behind the Figure 5E identity factors, study by study. One row per factor (two
# halves: targets higher or lower in SC-SIC), one dot per regulon and per intestinal study its network was not trained
# on: EC minus other
# enteroendocrine regulon target score, colored as in Figure 5E. Elmentaite 2020 + 2021 is the mean of its three parts,
# which are named under the key and not drawn (as in Figure 5E). Keys below the plot.

FA <- factor_rows()
BANDS <- row_bands(FA)
gh <- read_table("identity/f5_s_gh_values_ALL_ELIGIBLE.tsv")
# study means, in study, factor and regulon order (which sets each regulon's place inside its factor row)
GM <- gh[gh$part == "gut_study_mean", ]; GM <- GM[order(GM$item, GM$factor, GM$program, method = "radix"), ]
STUDY <- c(ELMENTAITE = "Elmentaite 2020 + 2021", HICKEY_2023 = "Hickey 2023", EGOZI_2021 = "Egozi 2021")
SUBLINE <- c(ELMENTAITE = "mean of 3 parts", HICKEY_2023 = "", EGOZI_2021 = "")
GM$study <- factor(STUDY[GM$item], STUDY); GM$sub <- factor(SUBLINE[GM$item], unique(SUBLINE)); GM$half <- FA$half[match(GM$factor, FA$factor)]
GM$y <- FA$y[match(GM$factor, FA$factor)] + ave(seq_len(nrow(GM)), GM$item, GM$factor, FUN = spread)
LIM <- max(abs(GM$value)) * 1.05
HB <- do.call(rbind, lapply(HALF, function(h) { yy <- FA$y[FA$half == h]; data.frame(half = factor(h, HALF), x = c(-LIM, LIM), y = c(min(yy) - .5, max(yy) + .5)) }))
pA <- ggplot() +
  geom_rect(data = BANDS, aes(xmin = -Inf, xmax = Inf, ymin = ymin, ymax = ymax), fill = "#D9D9D9", colour = NA) +
  geom_vline(xintercept = 0, colour = REF_COL, linewidth = .3) + geom_blank(data = HB, aes(x, y)) +
  geom_point(data = GM, aes(value, y, colour = note), size = 1) +
  facet_grid(rows = vars(half), cols = vars(study, sub), scales = "free_y", space = "free_y", switch = "y") +
  scale_y_continuous(breaks = FA$y, labels = parse(text = FA$expr), expand = expansion(add = c(.12, 0))) +
  scale_x_continuous(limits = c(-LIM, LIM), breaks = scales::breaks_pretty(4), expand = c(0, 0)) + coord_cartesian(clip = "off") +
  scale_colour_manual(values = REG_COL, labels = REG_LAB, breaks = names(REG_LAB), name = REG_KEY_TITLE) +
  guides(colour = guide_legend(override.aes = list(size = 1.6))) +
  labs(x = "EC minus other enteroendocrine (regulon target score); one dot per regulon", y = NULL,
       caption = "Elmentaite 2020 + 2021: mean of its 2020 fetal, 2021 fetal and 2021 adult parts, counted as one study
(+) activated targets, (−) repressed targets") +
  theme_panel("x", legend = "bottom") +
  theme(strip.placement = "outside", strip.text.y.left = element_text(angle = 0, hjust = 1, face = "bold", size = PT, lineheight = .95),
        strip.text.x = element_text(hjust = 0, vjust = 1, margin = margin(0, 0, 1.5, 0)), panel.spacing.x = unit(.2, "in"),
        legend.box = "vertical", legend.box.just = "left", legend.justification = "left",
        legend.spacing.y = unit(0, "in"), legend.margin = margin(0, 0, 0, 0), legend.title = element_text(size = PT, lineheight = .9),
        plot.caption = element_text(hjust = 0, size = PT, colour = "black", lineheight = 1.1, margin = margin(2, 0, 0, 0)),
        plot.margin = margin(2, 8, 5, 2))
# "mean of 3 parts" is a plain 7.5 pt second strip line under the bold study name, with 2.5 pt more air beneath it;
# ggplot2 styles every strip layer alike, so that one label is restyled in the built grob.
restyle <- function(x) {
  if (inherits(x, "text") && identical(as.character(x$label), SUBLINE[["ELMENTAITE"]])) {
    g0 <- unclass(x$gp); g0[c("font", "fontface", "fontsize")] <- NULL
    x$gp <- do.call(gpar, c(g0, list(fontface = "plain", fontsize = PT))); return(x)
  }
  if (inherits(x, "gtable")) { x$grobs <- lapply(x$grobs, restyle); return(x) }
  if (inherits(x, "gTree")) { x$children <- do.call(gList, lapply(x$children, restyle)); return(x) }
  x
}
draw_A <- function() {
  G <- ggplotGrob(pA); st <- grep("^strip-t", G$layout$name)
  G$grobs[st] <- lapply(G$grobs[st], function(s) { s <- restyle(s); s$heights[2] <- s$heights[2] + unit(2.5, "pt"); s })
  tr <- unique(G$layout$t[st]); G$heights[tr] <- G$heights[tr] + unit(2.5, "pt")
  grid.draw(G)
}

draw_panel <- function() place(draw_A, 185, 2478)
