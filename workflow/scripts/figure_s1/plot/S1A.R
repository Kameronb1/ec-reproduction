# Supplementary Figure S1A: islet hormone transcripts. Rows: intestinal EC, HIO graft EC and SC-SIC (cells with TPH1
# detected); SC-beta, -alpha, -delta and adult alpha, beta, delta (all cells). Columns: INS, GCG, SST, and a hormone of
# another islet cell type. One dot per study (single-nucleus studies open), a line spans the studies; row labels give
# the number of studies and cells. The table is written by the Figure 1 unit (serotonin_gene_count.R).

HG <- read_table("hormones_by_group.tsv", figure = "figure1")
ROWS_A <- data.frame(grp = c("Gut EC", "HIO graft EC", "SC-SIC", "SC-beta", "SC-alpha", "SC-delta", "Adult alpha", "Adult beta", "Adult delta"),
                     base = c("tph1", "tph1", "tph1", rep("all", 6)),
                     name = c("Intestinal EC, TPH1+", "HIO graft EC, TPH1+", "SC-SIC, TPH1+", "SC-beta", "SC-alpha", "SC-delta", "Adult alpha", "Adult beta", "Adult delta"),
                     col = COL[c("gut_ec", "hio_ec", "sc_sic", "beta", "alpha", "delta", "alpha", "beta", "delta")],
                     section = c("Intestine", "Intestine", rep("SC-islet", 4), rep("Adult islet", 3)), stringsAsFactors = FALSE)
SEC_A <- c("Intestine", "SC-islet", "Adult islet")
# y positions (row units, top = largest): a section header line, then its rows
yy <- numeric(0); ys <- numeric(0); y <- 0
for (s in SEC_A) { ys[s] <- y; y <- y - 1; for (k in which(ROWS_A$section == s)) { yy[k] <- y; y <- y - 1 }; y <- y - .25 }
ROWS_A$y <- yy
DA <- do.call(rbind, lapply(seq_len(nrow(ROWS_A)), function(i) { d <- HG[HG$grp == ROWS_A$grp[i] & HG$base == ROWS_A$base[i], ]
  if (!nrow(d)) stop("panel A: no rows for ", ROWS_A$grp[i])
  do.call(rbind, lapply(c("INS", "GCG", "SST", "other"), function(m) data.frame(row = ROWS_A$name[i], grp = d$grp, base = d$base, study = d$study,
    measure = m, pct = 100 * d[[m]], cells = d$cells, y = ROWS_A$y[i], col = ROWS_A$col[i], nuclei = d$study %in% NUCLEI))) }))
MEAS_A <- c(INS = "INS", GCG = "GCG", SST = "SST", other = "Hormone of\nanother\ncell type")   # column headers, drawn in grid (genes italic)
DA$facet <- factor(MEAS_A[DA$measure], MEAS_A)
RA <- do.call(rbind, lapply(split(DA, list(DA$row, DA$facet), drop = TRUE), function(d) data.frame(facet = d$facet[1], y = d$y[1], lo = min(d$pct), hi = max(d$pct))))
ROWS_A$n_studies <- vapply(ROWS_A$name, function(r) length(unique(DA$study[DA$row == r])), 0L)
ROWS_A$n_cells <- vapply(ROWS_A$name, function(r) { d <- DA[DA$row == r & DA$measure == "INS", ]; sum(d$cells) }, 0)
ROWS_A$label <- sprintf("%s (%d; %s)", ROWS_A$name, ROWS_A$n_studies, fmt(ROWS_A$n_cells))
YL_A <- c(min(ROWS_A$y) - .6, max(ys) + .6)
SEP_A <- ys[-1] + .6                                                      # faint rules above the 2nd and 3rd sections
pA_dots <- ggplot(DA) +
  geom_hline(yintercept = SEP_A, colour = "grey85", linewidth = .25) +
  geom_segment(data = RA, aes(x = lo, xend = hi, y = y, yend = y), colour = "grey75", linewidth = .35) +
  geom_point(data = DA[!DA$nuclei, ], aes(pct, y, fill = I(col)), shape = 21, colour = "black", size = 1.55, stroke = EDGE) +
  geom_point(data = DA[DA$nuclei, ], aes(pct, y), shape = 1, colour = "black", size = 1.9, stroke = EDGE) +
  geom_point(data = DA[DA$nuclei, ], aes(pct, y, colour = I(col)), shape = 1, size = 1.2, stroke = .5) +
  facet_wrap(~facet, nrow = 1) +
  scale_x_continuous(limits = c(-4, 104), breaks = c(0, 50, 100), labels = c("0", "50", "100"), expand = c(0, 0)) +
  scale_y_continuous(limits = YL_A, expand = c(0, 0)) +
  labs(x = "Cells with the transcript detected after ambient-RNA correction (%)", y = NULL) +
  theme_pg("x") + theme(axis.text.y = element_blank(), strip.text = element_blank(), panel.spacing.x = unit(.15, "in"), plot.margin = margin(.42 * 72, 3, 1, 1))
A_LAB_W <- 1.55
# the panel boxes (in, from the plot's top-left) of a ggplot drawn w x h on the open device
panels_in <- function(g, w, h) {
  pushViewport(viewport(width = unit(w, "in"), height = unit(h, "in")))
  wd <- vapply(seq_along(g$widths), function(i) convertWidth(g$widths[i], "in", valueOnly = TRUE), 0)
  ht <- vapply(seq_along(g$heights), function(i) convertHeight(g$heights[i], "in", valueOnly = TRUE), 0); popViewport()
  nw <- unitType(g$widths) == "null"; nh <- unitType(g$heights) == "null"
  wd[nw] <- (w - sum(wd[!nw])) * as.numeric(g$widths[nw]) / sum(as.numeric(g$widths[nw]))
  ht[nh] <- (h - sum(ht[!nh])) * as.numeric(g$heights[nh]) / sum(as.numeric(g$heights[nh]))
  cw <- c(0, cumsum(wd)); ch <- c(0, cumsum(ht)); pl <- g$layout[grepl("^panel", g$layout$name), ]; pl <- pl[order(pl$l), ]
  data.frame(l = cw[pl$l], r = cw[pl$r + 1], t = ch[pl$t], b = ch[pl$b + 1]) }
# the dot plot right of a label column; row labels, section headers and column headers drawn in grid
draw_A <- function(x, y, w, h) {
  xd <- x + A_LAB_W; wdots <- w - A_LAB_W; g <- ggplotGrob(pA_dots)
  pushViewport(viewport(unit(xd, "in"), unit(H - y, "in"), unit(wdots, "in"), unit(h, "in"), just = c("left", "top"))); grid.draw(g); popViewport()
  P <- panels_in(g, wdots, h); yof <- function(v) y + P$t[1] + (YL_A[2] - v) / diff(YL_A) * (P$b[1] - P$t[1])
  for (i in seq_len(nrow(ROWS_A))) {                                       # "TPH1" in the row labels is italic
    lab <- ROWS_A$label[i]
    e <- if (grepl("TPH1", lab, fixed = TRUE)) { p <- strsplit(lab, "TPH1", fixed = TRUE)[[1]]
      parse(text = sprintf("paste(%s, italic(\"TPH1\"), %s)", deparse(p[1]), deparse(p[2])))[[1]] } else lab
    grid.text(e, unit(xd - .05, "in"), unit(H - yof(ROWS_A$y[i]), "in"), just = c("right", "centre"), gp = gpar(fontfamily = FONT, fontsize = PT, col = "black")) }
  for (s in SEC_A) tx(s, x, yof(ys[[s]]), face = "bold", just = c("left", "centre"))
  for (k in seq_along(MEAS_A)) tx(MEAS_A[[k]], xd + (P$l[k] + P$r[k]) / 2, y + P$t[1] - .04, face = if (k <= 3) "bold.italic" else "bold", just = c("centre", "bottom"))
}

draw_panel <- function() draw_after_72dpi(function() {
  letter("A", LX, ROW1); draw_A(X0, ROW1, A_W, A_H)
  key_row(list(item("One study, whole cells", function(x, y) sym_dot(x, y, "grey45"), .14),
               item("One study, single nuclei", function(x, y) sym_dot(x, y, "grey45", open = TRUE), .14),
               item("Range of studies", function(x, y) sym_line(x, y, "grey75", lwd = .35), .26)), X0, A_W, KEY_A_Y)
})
