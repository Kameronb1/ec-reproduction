# Panels D and E of Figure 4 and their key (sourced by page.R).
#   D  promoter accessibility of each gene group against the typical gene promoter of the same cell type: Zhu 2023
#      (days 21, 32, 39) and Augsornworawat 2023 (weeks 2-4, months 6 and 12, graft median) on one y axis, the two
#      studies apart by a gap and named under the axis. Both tables hold the same quantity (chromatin_promoters_relative.R
#      and chromatin_promoters_culture_grafts.R): per gene, promoter (2 kb of a TSS) log2(1 + CPM) after thinning to the
#      cell type's shallowest point, minus the median of every gene with a promoter peak, then group median and quartiles.
#   E  distal accessibility change after grafting, against 10,000 random gene sets
GCOL5 <- c(G1 = "#43A5B6", G2 = "#A89A2E", G3 = "#B4553F")          # gene groups G1-G3
GCOL5_LIGHT <- vapply(GCOL5, function(c) grDevices::rgb(t(.4 * grDevices::col2rgb(c) + .6 * 255), maxColorValue = 255), "")
GL5 <- c("G1: EC genes low in SC-SICs", "G2: shared EC genes", "G3: separate SC-SIC program")
BAND5 <- "#D9D9D9"
pct5 <- function(d) 100 * (2^d - 1)                  # log2 difference -> % more
th5 <- function(grid_x = TRUE, grid_y = FALSE) theme_minimal(base_family = "Arial", base_size = 7.5) + theme(
  text = element_text(color = "black", size = 7.5), axis.text = element_text(color = "black", size = 7),
  axis.title = element_text(color = "black", size = 7), plot.title = element_text(color = "black", size = 7.5, face = "bold", margin = margin(0, 0, 2, 0)),
  strip.text = element_text(color = "black", size = 7, face = "bold", margin = margin(1, 0, 2, 0)),
  panel.grid.minor = element_blank(),
  panel.grid.major.x = if (grid_x) element_line(color = "#E4E4E4", linewidth = .2) else element_blank(),
  panel.grid.major.y = if (grid_y) element_line(color = "#E4E4E4", linewidth = .2) else element_blank(),
  axis.line.x = element_line(color = "black", linewidth = .35), axis.ticks.x = element_line(color = "black", linewidth = .3),
  axis.ticks.length = unit(1.3, "pt"), plot.margin = margin(2, 4, 2, 2), legend.position = "none", plot.title.position = "panel")
NOXL5 <- theme(axis.line.x = element_blank(), axis.ticks.x = element_blank())
zero_h5 <- geom_hline(yintercept = 0, color = "#777777", linewidth = .3)
# y titles
D_YLAB <- "Openness relative\nto the typical gene\npromoter of the same\ncell type (fold)"
E_YLAB <- "% more opening\nthan distal elements\nof the typical gene"
E_XLAB <- "Graft (week 2 → graft)"
DE_KEY <- c(median = "Median gene of the group", bars = "Middle half of the group genes", graft = "Graft median",
            random = "Middle 95% of 10,000 random gene sets")

# ---- D, Augsornworawat 2023: culture weeks 2-4, months 6 and 12, and the graft median, each population against its typical gene
POPL <- c("SC-SIC", "SC-beta", "SC-alpha"); TXP <- c("txp_1", "txp_2", "txp_3"); FLOOR <- -3.2
ce <- read_table("c5_culture_extended.tsv"); cw <- read_table("c5_culture_within_pops_wk2to4.tsv")
chk_w <- merge(ce[ce$sample %in% cw$sample, c("population", "sample", "g", "rel")], cw[, c("population", "sample", "g", "rel")], by = c("population", "sample", "g"))
if (nrow(chk_w) != nrow(cw) || max(abs(chk_w$rel.x - chk_w$rel.y)) > 0) stop("D: weeks 2-4 in the extended table differ from the weeks 2-4 table")
cgx <- read_table("c5_culture_extended_per_gene.tsv")
cs <- ce[grepl("^s6_", ce$sample), ]
chk_c <- merge(aggregate(rel ~ population + sample + g, cgx[grepl("^s6_", cgx$sample), ], median), cs[, c("population", "sample", "g", "rel")], by = c("population", "sample", "g"))
if (nrow(chk_c) != nrow(cs) || max(abs(chk_c$rel.x - chk_c$rel.y)) > 1e-9) stop("D: per-gene medians do not give back the drawn culture values")
chk_g <- merge(aggregate(rel ~ population + g, ce[ce$sample %in% TXP, ], median), ce[ce$sample == "graft_median", c("population", "g", "rel")], by = c("population", "g"))
if (nrow(chk_g) != 9L || max(abs(chk_g$rel.x - chk_g$rel.y)) > 1e-12) stop("D: graft median rows are not the median of the three grafts")
chk_gg <- merge(aggregate(rel ~ population + g + gene, cgx[cgx$sample %in% TXP, ], median), cgx[cgx$sample == "graft_median", c("population", "g", "gene", "rel")], by = c("population", "g", "gene"))
if (nrow(chk_gg) != sum(cgx$sample == "graft_median") || max(abs(chk_gg$rel.x - chk_gg$rel.y)) > 1e-12) stop("D: per-gene graft medians wrong")
if (!all(ce$drawn == !(ce$sample %in% TXP))) stop("D: drawn flag")
CX <- c(s6_week2 = 1, s6_week3 = 2.15, s6_week4 = 3.3, s6_month6 = 4.6, s6_month12 = 5.85, graft_median = 6.95); GDX <- c(G1 = -.2, G2 = 0, G3 = .2)
cu <- ce[ce$sample %in% names(CX), ]; cu$P <- factor(cu$population, POPL); cu$x <- unname(CX[cu$sample] + GDX[cu$g])
cg <- cgx[cgx$sample %in% names(CX), ]; cg$P <- factor(cg$population, POPL); cg$x <- unname(CX[cg$sample] + GDX[cg$g])
cq <- do.call(rbind, lapply(split(cg, interaction(cg$g, cg$sample, cg$P, drop = TRUE)), function(d)
  data.frame(g = d$g[1], P = d$P[1], x = d$x[1], lo = unname(quantile(d$rel, .25)), hi = unname(quantile(d$rel, .75)))))
# ---- D, Zhu 2023, on the same y axis -----------------------------------------------------------------------------------------
# Same quantity in both studies (producers chromatin_promoters_relative.R and chromatin_promoters_culture_grafts.R): per
# gene, promoter accessibility (peaks within 2 kb of a TSS, log2(1 + CPM)) after thinning each cell type to its shallowest
# time point, minus the median of every gene with a promoter peak in the same sample; per group the median gene and the
# 25th-75th percentiles. They differ in annotation (Zhu: UCSC hg19 RefSeq TSS; Augsornworawat: GENCODE v44 hg38 TSS) and
# in what is thinned (Zhu: total accessible events; Augsornworawat: peak counts).
zb <- read_table("c5_relative_zhu_hickey.tsv"); zb <- zb[zb$study == "Zhu 2023", ]
ZX <- c(`21` = 1, `32` = 2.15, `39` = 3.3)                   # day positions
if (nrow(zb) != 27L || min(zb$q25) <= FLOOR) stop("D: Zhu 2023 rows missing, or a bar reaches the axis floor")
zb$P <- factor(zb$population, POPL); zb$x <- unname(ZX[as.character(zb$day)] + GDX[zb$g])
SHIFT <- 4.0                                                 # Augsornworawat 2023 starts 1.7 units after Zhu's last day (a gap)
cu$x <- cu$x + SHIFT; cq$x <- cq$x + SHIFT; CXS <- CX + SHIFT
D_XLIM <- c(.62, 7.3 + SHIFT)
# marks under the axis, at fixed inch distances below the panel's lower edge (the "weeks" / "months" offsets of a
# reference plot of height F5_PANEL_H; the study names one line lower)
D_YLIM <- c(FLOOR, 1.75); D_PAD <- .05 * diff(D_YLIM); D_YRANGE <- diff(D_YLIM) + 2 * D_PAD   # coord_cartesian's 5% expansion
D_MARK_F5 <- c(rule = FLOOR - .92 * 4.95 / 4.45, text = FLOOR - 4.95 / 4.45)                  # their y values on the reference plot
d_mark_y <- function(h, f5_h) { inch <- ((D_YLIM[1] - D_PAD) - D_MARK_F5) * f5_h / D_YRANGE   # inches below the lower edge
  inch <- c(inch, study = inch[["text"]] + .125); (D_YLIM[1] - D_PAD) - inch * D_YRANGE / h }
D_SPANS <- data.frame(x = c(ZX[[1]], CXS[[1]], CXS[[4]]) - .25, xend = c(ZX[[3]], CXS[[3]], CXS[[5]]) + .25)
D_SPAN_LAB <- data.frame(x = c(mean(ZX), mean(CXS[1:3]) - .3, mean(CXS[4:5]) + .3), lab = c("days", "weeks", "months"))
D_STUDY <- data.frame(x = c(mean(ZX), mean(CXS[c(1, 6)])), lab = c("Zhu\n2023", "Augsornworawat\n2023"))   # two lines: on one line the two names do not fit a facet
d_plot <- function(h, f5_h) { my <- d_mark_y(h, f5_h)
  ggplot(cu, aes(x, rel, color = g)) + zero_h5 +
  geom_vline(xintercept = c(3.95, 6.43) + SHIFT, color = "#BBBBBB", linewidth = .2, linetype = "22") +
  geom_rect(data = cq, aes(xmin = x - .09, xmax = x + .09, ymin = lo, ymax = hi), fill = unname(GCOL5_LIGHT[cq$g]), colour = NA, inherit.aes = FALSE) +
  geom_rect(data = zb, aes(xmin = x - .09, xmax = x + .09, ymin = q25, ymax = q75), fill = unname(GCOL5_LIGHT[zb$g]), colour = NA, inherit.aes = FALSE) +
  geom_line(data = cu[cu$sample != "graft_median", ], aes(group = g), linewidth = .45) +
  geom_line(data = zb, aes(group = g), linewidth = .45) +
  geom_point(data = cu[cu$sample != "graft_median", ], aes(fill = g), shape = 21, colour = "black", stroke = .15, size = 1.3) +
  geom_point(data = zb, aes(fill = g), shape = 21, colour = "black", stroke = .15, size = 1.3) +
  geom_point(data = cu[cu$sample == "graft_median", ], shape = 23, fill = "white", stroke = .6, size = 1.6) +
  geom_segment(data = D_SPANS, aes(x = x, xend = xend, y = my[["rule"]], yend = my[["rule"]]), colour = "black", linewidth = .2, inherit.aes = FALSE) +
  geom_text(data = D_SPAN_LAB, aes(x = x, y = my[["text"]], label = lab), vjust = 1, size = 7 / ggplot2::.pt, family = "Arial", colour = "black", inherit.aes = FALSE) +
  geom_text(data = D_STUDY, aes(x = x, y = my[["study"]], label = lab), vjust = 1, lineheight = .9, size = 7 / ggplot2::.pt, family = "Arial", fontface = "bold", colour = "black", inherit.aes = FALSE) +
  facet_wrap(~P, nrow = 1) + scale_color_manual(values = GCOL5) + scale_fill_manual(values = GCOL5) +
  scale_x_continuous(breaks = c(unname(ZX), unname(CXS[1:5])), labels = c("21", "32", "39", "2", "3", "4", "6", "12"), limits = D_XLIM) +
  scale_y_continuous(breaks = log2(c(.125, .25, .5, 1, 2)), labels = c("×0.125", "×0.25", "×0.5", "×1", "×2")) + coord_cartesian(ylim = D_YLIM, clip = "off") +
  labs(x = NULL, y = D_YLAB) + th5(FALSE, TRUE) +
  theme(panel.spacing = unit(8, "pt"), axis.text.x = element_text(size = 7), plot.margin = margin(2, 4, 31, 2)) }

# ---- E: grafts, change beyond the typical gene's, against 10,000 random gene sets -------------------------
WIN <- if (all(as.character(read_table("c5_graft_window_agreement.tsv")$same) %in% c("TRUE", "True", "true"))) "promoter" else "distal"
gr <- read_table("c5_graft_relative_by_population.tsv"); gr <- gr[gr$window == WIN, ]; gr$y <- pct5(gr$group_rel); gr$lo <- pct5(gr$null_q025_rel); gr$hi <- pct5(gr$null_q975_rel)
gr$u <- factor(c(txp_1 = "1", txp_2 = "2", txp_3 = "3")[gr$graft], c("1", "2", "3")); gr$win <- factor(gr$population, c("SC-SIC", "SC-beta", "SC-alpha"))
pE <- ggplot(gr, aes(u, y, group = g)) + zero_h5 +
  geom_linerange(aes(ymin = lo, ymax = hi), position = position_dodge(.8), linewidth = .85, color = BAND5) +
  geom_point(aes(fill = g), shape = 21, colour = "black", stroke = .15, size = 1.2, position = position_dodge(.8)) + facet_wrap(~win, nrow = 1) +
  scale_fill_manual(values = GCOL5) + th5(FALSE, TRUE) + NOXL5 +
  labs(x = E_XLAB, y = E_YLAB) +
  theme(panel.spacing = unit(5, "pt"), axis.text.x = element_text(size = 7))

# ---- the D-E key -------------------------------------------------------------------------------------------------------------
PT_LWD5 <- .15 * ggplot2::.stroke / 2
`%||%` <- function(a, b) if (is.null(a)) b else a
row_keys5 <- function(items, left, width, y, gap = .16, centre = TRUE) {   # key items of types pt, vbar, bars3, pb
  ww <- vapply(items, function(it) str_in(it$label) + if (it$type == "bars3") .26 else .14, 0); total <- sum(ww) + gap * (length(ww) - 1)
  x <- if (centre) left + (width - total) / 2 else left
  for (i in seq_along(items)) { it <- items[[i]]; yy <- H - y - .045
    if (it$type == "pt") grid.points(unit(x + .05, "in"), unit(yy, "in"), pch = it$pch %||% 16, size = unit(it$size %||% 1.5, "mm"), gp = gpar(col = it$col, fill = it$fill %||% it$col, lwd = it$plwd %||% .8))
    if (it$type == "vbar") grid.rect(unit(x + .05, "in"), unit(yy, "in"), width = unit(it$w %||% .065, "in"), height = unit(it$h %||% .12, "in"), gp = gpar(fill = it$col, col = NA))
    if (it$type == "bars3") for (k in 1:3)
      grid.rect(unit(x + .025 + (k - 1) * .07, "in"), unit(yy, "in"), width = unit(.05, "in"), height = unit(it$h %||% .12, "in"), gp = gpar(fill = it$col[k], col = NA))
    if (it$type == "pb") { grid.rect(unit(x + .05, "in"), unit(yy, "in"), width = unit(.05, "in"), height = unit(.13, "in"), gp = gpar(fill = it$bar, col = NA))
      grid.points(unit(x + .05, "in"), unit(yy, "in"), pch = 21, size = unit(1.5, "mm"), gp = gpar(col = "black", fill = it$fill, lwd = it$plwd)) }
    text_at(it$label, x + if (it$type == "bars3") .26 else .14, y - .005)
    x <- x + ww[i] + gap }
  invisible(total) }
de_key <- function(left, width, y1, y2) {
  row_keys5(lapply(1:3, function(k) list(type = "pb", label = GL5[k], fill = GCOL5[[k]], bar = GCOL5_LIGHT[[k]], plwd = PT_LWD5)), left, width, y1, gap = .2)
  row_keys5(list(list(type = "pt", pch = 21, label = DE_KEY[["median"]], col = "black", fill = "grey75", plwd = PT_LWD5),
                 list(type = "bars3", label = DE_KEY[["bars"]], col = unname(GCOL5_LIGHT), h = .09),
                 list(type = "pt", pch = 23, label = DE_KEY[["graft"]], col = "black", fill = "white", plwd = 1.2),
                 list(type = "vbar", label = DE_KEY[["random"]], col = BAND5, w = .85 * ggplot2::.pt / 96, h = .09)), left, width, y2, gap = .2) }
