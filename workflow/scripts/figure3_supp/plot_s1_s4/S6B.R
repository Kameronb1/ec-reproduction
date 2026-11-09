# Figure S6B: Figure 2E's scatter for the six SC-islet studies. Every SC endocrine cell, its four-marker score (x)
# against its score from Figure 2E's 20-gene classifier (intestinal genes, intestinal training; y), one plot per study,
# both cutoffs dashed; SC-SIC drawn over the other cells. Corner numbers, as in Figure 2E: percent of each label's cells
# classified as EC by both classifiers (top right), by the 20-gene classifier only (top left), by the four-marker
# classifier only (bottom right) and by neither (bottom left); SC-SIC in magenta, other cells in grey; whole percent,
# "<1%" and ">99%" unless exactly none or all. Each classifier also applies its neighborhood check, so a dot past a cutoff
# can be unclassified; these are the same classifications as A's four-marker bar and its 20-gene bar with intestinal
# genes and intestinal training. Augsornworawat 2023 (nuclei), with many more cells, has smaller, fainter dots. Key below.

cells <- read_table("s2_classifier_swap/sc_cells.tsv.gz")
cells$call4 <- as.logical(cells$call4); cells$call20 <- as.logical(cells$call20)
# check: the cell table gives back A's counts for the two classifiers drawn here
calls <- read_table("s2_classifier_swap/calls_by_study.tsv")
for (s0 in names(STUDY)) for (cl in c("SC_SIC", "OTHER")) {
  x <- cells[cells$study == s0 & cells$class == cl, ]
  a <- calls[calls$held_out_study == s0 & calls$class == if (cl == "SC_SIC") "SC_SIC" else "OTHER_MATURE_SC_ENDOCRINE", ]
  stopifnot(nrow(x) == a$n[a$combination == "GUTG_GUTT"], sum(x$call4) == a$called[a$combination == "MARKER_RULE_GATED"],
            sum(x$call20) == a$called[a$combination == "GUTG_GUTT"])
}
CUT4 <- unique(cells$cut4); THR20 <- unique(cells$thr20); stopifnot(length(CUT4) == 1, length(THR20) == 1)
B_COL <- c(SC_SIC = SIC, OTHER = "#9A9A9A")
B_XR <- c(-1.5, 2.6); B_YR <- range(cells$score20); B_YR <- B_YR + c(-.22, .22) * diff(B_YR)   # as Figure 2E
scatter_B <- function(s0, first) {
  x <- cells[cells$study == s0, ]; x <- x[order(x$class == "SC_SIC"), ]; big <- nrow(x) > 10000
  q <- function(k) { y <- x[x$class == k, ]; a4 <- y$call4; a20 <- y$call20
    c(sum(a4 & a20), sum(!a4 & a20), sum(a4 & !a20), sum(!a4 & !a20)) }      # EC by both, 20-gene only, four-marker only, neither
  qe <- q("SC_SIC"); qo <- q("OTHER"); dy <- .075 * diff(B_YR); ins <- .03 * diff(B_YR)
  xx <- c(B_XR[2] - .1, B_XR[1] + .1, B_XR[2] - .1, B_XR[1] + .1); h <- c(1, 0, 1, 0)
  ann <- rbind(data.frame(x = xx, y = c(B_YR[2] - ins, B_YR[2] - ins, B_YR[1] + ins, B_YR[1] + ins), h = h, v = c(1, 1, 0, 0), lab = qe, col = B_COL[["SC_SIC"]]),
               data.frame(x = xx, y = c(B_YR[2] - ins - dy, B_YR[2] - ins - dy, B_YR[1] + ins + dy, B_YR[1] + ins + dy), h = h, v = c(1, 1, 0, 0),
                          lab = qo, col = B_COL[["OTHER"]]))
  pc <- function(q) { p <- 100 * q / sum(q); ifelse(q == 0, "0%", ifelse(p < 1, "<1%", ifelse(p > 99 & q < sum(q), ">99%", paste0(round(p), "%")))) }
  ann$lab <- c(pc(qe), pc(qo))
  n <- table(x$class)
  ggplot(x, aes(score4, score20)) +
    geom_point(aes(colour = B_COL[class]), size = if (big) .22 else .3, stroke = 0, alpha = if (big) .3 else .6) +
    geom_vline(xintercept = CUT4, linetype = "dashed", linewidth = .3) + geom_hline(yintercept = THR20, linetype = "dashed", linewidth = .3) +
    geom_text(data = ann, aes(x, y, label = lab, hjust = h, vjust = v, colour = col), size = 6.5 / .pt, family = FONT, fontface = "bold", inherit.aes = FALSE) +
    scale_colour_identity() + scale_x_continuous(limits = B_XR, breaks = -1:2, expand = c(0, 0)) + scale_y_continuous(limits = B_YR, expand = c(0, 0)) +
    labs(x = "Four-marker score", y = if (first) "20-gene score" else NULL, title = sub(" (nuclei)", "", STUDY[[s0]], fixed = TRUE),
         subtitle = sprintf("%s%s | %s", if (grepl("nuclei", STUDY[[s0]])) "(nuclei) " else "", format(n[["SC_SIC"]], big.mark = ","),
                            format(n[["OTHER"]], big.mark = ","))) +
    theme_rk("x") +
    theme(panel.border = element_rect(fill = NA, colour = "black", linewidth = .3), axis.line = element_blank(), panel.grid.major = element_blank(),
          panel.grid.major.x = element_blank(), axis.ticks = element_line(colour = "black", linewidth = .3),
          plot.title = element_text(size = 7, face = "bold", hjust = 0, margin = margin(0, 0, 0, 0)),
          plot.subtitle = element_text(size = 6.5, hjust = 0, margin = margin(0, 0, 1, 0)),
          plot.margin = margin(2, 3, 2, 2)) + (if (first) NULL else theme(axis.text.y = element_blank()))
}
plot_B <- wrap_plots(lapply(seq_along(STUDY), function(i) scatter_B(names(STUDY)[i], i == 1)), nrow = 1)

B_KEY_H <- .36                                                    # key band under the plots (inches)
draw_B <- function() {
  vp_h <- convertHeight(unit(1, "npc"), "in", valueOnly = TRUE)
  pushViewport(viewport(y = unit(B_KEY_H, "in"), height = unit(vp_h - B_KEY_H, "in"), just = "bottom")); print(plot_B, newpage = FALSE); popViewport()
  tp <- gpar(fontfamily = FONT, fontsize = PT, col = "black"); y1 <- B_KEY_H - .09; y2 <- y1 - .17
  dot <- function(x, col) grid.points(unit(x, "in"), unit(y1, "in"), pch = 16, size = unit(.075, "in"), gp = gpar(col = col))
  tx <- function(s, x, y) grid.text(s, unit(x, "in"), unit(y, "in"), just = c("left", "centre"), gp = tp)
  dot(.45, B_COL[["SC_SIC"]]); tx("SC-SIC", .55, y1); dot(1.15, B_COL[["OTHER"]]); tx("Other SC endocrine cells", 1.25, y1)
  grid.lines(unit(c(2.75, 2.95), "in"), unit(y1, "in"), gp = gpar(lty = "dashed", lwd = .3 * .pt))
  tx("Cutoffs: four-marker (Fig. 2A); 20-gene, intestinal genes and training (= Fig. 2E)", 3.02, y1)
  tx("Under each study: SC-SIC | other cells. Numbers: percent of each label's cells called by both, by one or by neither classifier", .45, y2)
}

draw_panel <- function() in_place(120, 4248, 7.3, 2.05 + B_KEY_H - .30, 1265, draw_B)
