# Figure 4B: each graft's share of mature endocrine cells that are SC-SIC, SC-beta and SC-alpha (primary study labels),
# against the in vitro samples its paper compares it with: a grey bar spans the in vitro samples' range, and an arrow
# runs from their mean (open circle) to the graft (open arrowhead: a graft with fewer than 50 mature endocrine cells).
# One facet row per graft source, labeled at the left (the three Balboa 2022 3-month grafts pooled into one row);
# one column per cell type. Key below.

B_X <- 2.50; B_KEY_TOP <- 3.42; B_KEY_PITCH <- .135          # left edge of the source labels; top of the key; its line pitch
B_Y <- .10 + .40 - .04 - 119 / 600; B_H <- 2.62 + .30 + 119 / 600   # top and height of the plot
GCOL <- c(`SC-SIC` = "#AA3377", `SC-alpha` = "#117733", `SC-beta` = "#4477AA")

ab <- read_table("f4b_composition.tsv")
# The graft sources (facet rows), as label lines of bold (b) and plain (p) runs; each graft and its in vitro samples.
B_STRIP <- list(list(c(b = "Augsornworawat 2020"), c(p = "HUES8, from stage 6")), list(c(b = "Augsornworawat 2020"), c(p = "WS4corr, from stage 6")),
                list(c(b = "Augsornworawat 2023"), c(p = "HUES8 (nuclei)"), c(p = "from stage 6 wk 2")),
                list(c(b = "Balboa 2022"), c(p = "H1, from stage 7 wk 0")),
                list(c(b = "Balboa 2022"), c(p = "H1, from stage 7 wk 3")),
                list(c(b = "Du 2022", p = ", stage 6 (vs day 2)")))
names(B_STRIP) <- vapply(B_STRIP, function(L) paste(vapply(L, paste, "", collapse = ""), collapse = "\n"), "")
SRC <- data.frame(study = c("Augsornworawat2020", "Augsornworawat2020", "Augsornworawat2023", "Balboa2022", "Balboa2022", "Balboa2022", "Balboa2022", "Du2022"),
                  graft = c("M6_graft_HUES8", "M6_graft_WS4corr", "M6_graft", "M1_graft", "M3_graft_pooled", "M6_graft_CTRL", "M6_graft_STZ", "10_wpt"),
                  source = c("S6_in_vitro", "WS4_baseline", "sc_islet_1|sc_islet_2|S6_week2", "S7w0", "S7w0", "S7w3", "S7w3", "S6D2"),
                  strip = names(B_STRIP)[c(1, 2, 3, 4, 4, 5, 5, 6)])
GC <- c("SC-SIC", "SC-alpha", "SC-beta")
# one row per graft sample and cell type: its share, and the in vitro samples' mean, minimum and maximum
PBs <- do.call(rbind, lapply(seq_len(nrow(SRC)), function(i) {
  s <- ab[ab$study == SRC$study[i] & ab$context %in% strsplit(SRC$source[i], "|", fixed = TRUE)[[1]] & ab$category %in% GC, ]
  gr <- ab[ab$study == SRC$study[i] & ab$context == SRC$graft[i] & ab$category %in% GC, ]
  do.call(rbind, lapply(GC, function(cg) { sv <- s$pct[s$category == cg]; gv <- gr[gr$category == cg, ]
    data.frame(strip = SRC$strip[i], context = SRC$graft[i], row_lab = gv$row_lab, mature = gv$mature, category = cg, pct = gv$pct,
               src_mean = mean(sv), src_lo = min(sv), src_hi = max(sv)) })) }))
# graft labels: shortened sample label and its mature endocrine cells
short_lab <- function(x) { x <- sub(" \\(from S7 wk [03]\\)", "", x); x <- sub(" \\((HUES8|WS4corr)\\)", "", x)
  x <- sub("6-mo graft, ", "6-mo, ", x); x <- sub("10-week graft", "10-wk graft", x); x <- sub(" grafts, pooled", ", pooled", x); x <- sub(" graft", "", x); gsub("  #", " #", x) }
PBs$glab <- sprintf("%s (%s)", short_lab(PBs$row_lab), trimws(format(PBs$mature, big.mark = ",")))
PBs$strip <- factor(PBs$strip, unique(SRC$strip)); PBs$category <- factor(PBs$category, c("SC-SIC", "SC-beta", "SC-alpha"))
PBs$glab <- factor(PBs$glab, rev(unique(PBs[order(PBs$strip, PBs$context, PBs$row_lab, method = "radix"), "glab"])))
PBs$low <- PBs$mature < 50
B_TEXT <- c(x = "Share of mature endocrine cells (%, primary study labels)", header = "Graft (mature endocrine cells)",
            range = "in vitro samples: range, mean", graft = "graft", low = "graft with < 50 mature endocrine cells")

# marks = FALSE: the same plot without its marks, to measure the layout
b_plot <- function(marks) {
  p <- ggplot(PBs, aes(y = glab)) + geom_linerange(aes(xmin = src_lo, xmax = src_hi), linewidth = 2.6, color = "grey75") +
    geom_blank(aes(x = src_mean)) + geom_blank(aes(x = pct))
  if (marks) p <- p +                                         # shaft, head (open when < 50 cells), then the start circle on top
    geom_segment(data = PBs[!PBs$bare, ], aes(x = src_mean, xend = shaft_end, yend = glab, color = category), linewidth = .45, lineend = "butt") +
    geom_segment(data = PBs[!PBs$low, ], aes(x = hbase, xend = pct, yend = glab, color = category), linewidth = .45, arrow = ARW) +
    geom_segment(data = PBs[PBs$low, ], aes(x = hbase, xend = pct, yend = glab, color = category), linewidth = .30, arrow = ARW, arrow.fill = "white") +
    geom_point(data = PBs[!PBs$bare, ], aes(x = src_mean, color = category), shape = 21, fill = "white", size = PT_SIZE, stroke = PT_STROKE) +
    scale_color_manual(values = GCOL) + scale_fill_manual(values = GCOL, guide = "none")
  p + facet_grid(strip ~ category, scales = "free_y", space = "free_y", switch = "y") +
    scale_x_continuous(limits = c(0, 82), breaks = c(0, 25, 50, 75), expand = c(0, 0)) +
    labs(x = B_TEXT[["x"]], y = NULL) + th() +
    theme(strip.placement = "outside", strip.text.y.left = element_blank(), strip.text.x = element_text(size = 7.5, face = "bold"),
          panel.spacing.y = unit(2.5, "pt"), panel.spacing.x = unit(10, "pt"), panel.border = element_rect(fill = NA, color = "#BDBDBD", linewidth = .25),
          axis.text.y = element_text(size = 7)) }

# The source-label column is as wide as its widest line; the plot starts .06 in after it and ends at 7.44 in.
B_LINES <- unlist(B_STRIP, recursive = FALSE)
run_w <- function(r) sum(vapply(seq_along(r), function(i) str_in(r[[i]], 7, if (names(r)[i] == "b") "bold" else "plain"), 0))
B_ASC <- max(unlist(lapply(B_LINES, function(r) mapply(function(s, b) asc_in(s, 7, if (b) "bold" else "plain"), r, names(r) == "b"))))
B_DSC <- max(unlist(lapply(B_LINES, function(r) mapply(function(s, b) dsc_in(s, 7, if (b) "bold" else "plain"), r, names(r) == "b"))))
B_STRIP_W <- max(vapply(B_LINES, run_w, 0)); B_PLOT_X <- B_X + B_STRIP_W + .06; B_PLOT_W <- 7.44 - (B_X + B_STRIP_W + .06)
# inches per percentage point, so the arrows can be shortened by their heads in data units
gB0 <- ggplotGrob(b_plot(FALSE)); nnB <- unitType(gB0$widths) == "null"
B_IPU <- (B_PLOT_W - sum(convertWidth(gB0$widths[!nnB], "in", valueOnly = TRUE))) / 3 / diff(ggplot_build(b_plot(FALSE))$layout$panel_params[[1]]$x.range)
# The shaft stops where the head starts; the head sits on a 0.2 px stub; a change shorter than the head is a bare head.
PBs$s <- sign(PBs$pct - PBs$src_mean); PBs$hbase <- PBs$pct - PBs$s * STUB_IN / B_IPU
PBs$shaft_end <- PBs$pct - PBs$s * HEAD_AX0 / B_IPU; PBs$bare <- abs(PBs$pct - PBs$src_mean) <= HEAD_AX0 / B_IPU
pB <- b_plot(TRUE)

draw_panel <- function() shifted_up(ROW_UP[["ab"]], {
  g <- ggplotGrob(pB); draw_at(g, B_PLOT_X, B_Y, B_PLOT_W, B_H)
  # the source labels, right-aligned, each centered on its facet row (the first column's panels, top to bottom)
  P <- panels_in(g, B_PLOT_X, B_Y, B_PLOT_W, B_H); col1 <- vapply(P, `[[`, 0, "l") == min(vapply(P, `[[`, 0, "l"))
  P <- P[col1][order(vapply(P[col1], `[[`, 0, "t"))]
  x_r <- B_X + B_STRIP_W
  for (r in seq_along(B_STRIP)) { bx <- P[[r]]; L <- B_STRIP[[r]]; n <- length(L)
    top_ink <- bx$top_in + bx$height_in / 2 - (B_ASC + (n - 1) * 7.5 / 72 + B_DSC) / 2
    for (k in seq_len(n)) { runs <- L[[k]]; base <- top_ink + B_ASC + (k - 1) * 7.5 / 72; xr <- x_r
      for (j in rev(seq_along(runs))) { f <- if (names(runs)[j] == "b") "bold" else "plain"; text_bl(runs[[j]], xr, base, face = f, just = "right"); xr <- xr - str_in(runs[[j]], 7, f) } } }
  # the graft-label header, right-aligned with the graft labels, on the line of the cell-type names
  la <- unique(g$layout$l[grepl("^axis-l", g$layout$name)]); st <- unique(g$layout$t[grepl("^strip-t", g$layout$name)])
  x_right <- B_PLOT_X + sum(convertWidth(g$widths[seq_len(la - 1L)], "in", valueOnly = TRUE)) + max(vapply(levels(PBs$glab), str_in, 0, size = 7))
  base <- P[[1]]$top_in - (convertHeight(g$heights[st], "in", valueOnly = TRUE) / 2 + .5 / 72.27 - asc_in("SC-SIC", 7.5, "bold") / 2)
  text_bl(B_TEXT[["header"]], x_right, base, just = "right")
  # key, two lines centered under the plots: the in vitro range bar with its mean (open circle) and the graft arrow; then the
  # open-headed arrow. Each glyph's text .06 in after it, .20 in between items.
  PP <- panels_in(g, B_PLOT_X, B_Y, B_PLOT_W, B_H); kc <- (min(vapply(PP, `[[`, 0, "left_in")) + max(vapply(PP, function(p) p$left_in + p$width_in, 0))) / 2
  gl <- c(range = .20, graft = .18, low = .18); kline <- c(range = 1, graft = 1, low = 2); K <- NULL
  if (!any(PBs$low)) { gl <- gl[1:2]; kline <- kline[1:2] }          # the open-headed arrow is keyed only when drawn
  for (ln in unique(kline)) { it <- names(gl)[kline == ln]; x <- kc - (sum(gl[it] + .06 + vapply(B_TEXT[it], str_in, 0)) + .20 * (length(it) - 1)) / 2
    for (k in it) { tl <- x + gl[[k]] + .06; K <- rbind(K, data.frame(item = k, gl = x, gr = x + gl[[k]], tl = tl, top = B_KEY_TOP + (ln - 1) * B_KEY_PITCH)); x <- tl + str_in(B_TEXT[[k]]) + .20 } }
  yc <- K$top + .045; yy <- function(i) unit(rep(H - yc[i], 2), "in")
  grid.rect(x = unit(K$gl[1], "in"), y = unit(H - yc[1], "in"), width = unit(K$gr[1] - K$gl[1], "in"), height = unit(2.6 * ggplot2::.pt / 96, "in"), just = c("left", "center"), gp = gpar(fill = "grey75", col = NA))
  grid.points(x = unit((K$gl[1] + K$gr[1]) / 2, "in"), y = unit(H - yc[1], "in"), pch = 21, size = unit(CIRCLE_PT, "pt"), gp = gpar(col = "#555555", fill = "white", lwd = CIRCLE_LWD))
  grid.lines(x = unit(c(K$gl[2], K$gr[2]), "in"), y = yy(2), arrow = ARW, gp = gpar(col = "#555555", fill = "#555555", lwd = 1.3))
  if (nrow(K) == 3) grid.lines(x = unit(c(K$gl[3], K$gr[3] - HEAD_AX0), "in"), y = yy(3), gp = gpar(col = "#555555", lwd = 1.3, lineend = "butt"))
  if (nrow(K) == 3) grid.lines(x = unit(K$gr[3] - c(STUB_IN, 0), "in"), y = yy(3), arrow = ARW, gp = gpar(col = "#555555", fill = "white", lwd = .30 * ggplot2::.pt))
  for (i in seq_len(nrow(K))) text_at(B_TEXT[[K$item[i]]], K$tl[i], K$top[i] - .005)
})
