# Drawing helpers shared by the Figure 4 panels: text size, colors, theme, placement on the page, the eleven rows that
# panels A and C-F share, the arrow marks and key glyphs, and the measured layouts of the C-D and E-F rows. Sourced by
# workflow/scripts/common/draw_panel.R, which defines the page size W x H (inches) and DPI. Every position is in inches
# from the page's top-left corner.

PT <- 7.5                                                      # base text size (points)
POP_COL <- c(SC_SIC = "#AA3377", SC_BETA = "#4477AA", SC_ALPHA = "#117733")
POP_LAB <- c(SC_SIC = "SC-SIC", SC_BETA = "SC-beta", SC_ALPHA = "SC-alpha")
ARW_LEN <- .042; ARW_ANGLE <- 25                               # every arrowhead: .042 in sides, 25 degrees, closed
ARW <- arrow(length = unit(ARW_LEN, "in"), type = "closed", angle = ARW_ANGLE)
HEAD_AX0 <- ARW_LEN * cos(ARW_ANGLE * pi / 180)                # the head's length along its axis
HEAD_MIN_IN <- ARW_LEN + 2 / DPI                               # shorter changes are drawn as a bare arrowhead
STUB_IN <- 1 / (5 * DPI)                                       # a bare head sits on a 0.2 px stub
PT_SIZE <- .95; PT_STROKE <- .45                               # the open circle at an arrow's start
CIRCLE_PT <- PT_SIZE * ggplot2::.pt + PT_STROKE * ggplot2::.stroke / 2; CIRCLE_LWD <- PT_STROKE * ggplot2::.stroke / 2
CIRCLE_R <- .375 * CIRCLE_PT / 72 + CIRCLE_LWD / 96 / 2        # its outer radius
BAND_GRAY <- "#E8E8E8"                                         # bands on alternate rows (A, C-F)
BAND_ALPHA <- 1 - grDevices::col2rgb(BAND_GRAY)[1, 1] / 255    # black at this alpha renders BAND_GRAY on white

# The three rows of panels (A-B, C-D and their key, E-F and their key) are laid out at fixed places and drawn higher by
# these whole-pixel shifts, which close the space under each row's panel letters.
ROW_UP <- c(ab = 127, cd = 220, ef = 318) / 600
shifted_up <- function(dy, expr) {                             # draw expr dy inches higher
  pushViewport(viewport(x = 0, y = unit(dy, "in"), width = unit(W, "in"), height = unit(H, "in"), just = c("left", "bottom")))
  force(expr); popViewport() }

th <- function() theme_minimal(base_family = "Arial", base_size = PT) + theme(
  text = element_text(color = "black", size = PT), axis.text = element_text(color = "black", size = 7),
  axis.title = element_text(color = "black", size = 7), plot.title = element_text(color = "black", size = 7.5, face = "bold", margin = margin(0, 0, 2, 0)),
  strip.text = element_text(color = "black", size = 7.5, face = "bold", margin = margin(1, 0, 2, 0)),
  panel.grid.minor = element_blank(), panel.grid.major.y = element_blank(), panel.grid.major.x = element_line(color = "#E4E4E4", linewidth = .2),
  axis.line.x = element_line(color = "black", linewidth = .35), axis.ticks.x = element_line(color = "black", linewidth = .3),
  axis.ticks.length = unit(1.3, "pt"), plot.margin = margin(2, 4, 2, 2), legend.position = "none", plot.title.position = "panel")
nobg <- function(p) p + theme(plot.background = element_blank())   # C-F: the row bands show through the plots

# ---- text metrics on the open device (inches) ----
fm <- function(s, size, face, f) { pushViewport(viewport(gp = gpar(fontfamily = "Arial", fontsize = size, fontface = face))); v <- f(s); popViewport(); v }
str_in <- function(s, size = PT, face = "plain") fm(s, size, face, function(s) convertWidth(stringWidth(s), "in", valueOnly = TRUE))
asc_in <- function(s, size = PT, face = "plain") fm(s, size, face, function(s) convertHeight(unit(1, "strheight", s), "in", valueOnly = TRUE))
dsc_in <- function(s, size = PT, face = "plain") fm(s, size, face, function(s) convertHeight(unit(1, "strdescent", s), "in", valueOnly = TRUE))

# ---- drawing in page inches ----
text_at <- function(s, x, y, size = PT, face = "plain", col = "black")      # text with its top at y
  grid.text(s, x = unit(x, "in"), y = unit(H - y, "in"), just = c("left", "top"), gp = gpar(fontfamily = "Arial", fontsize = size, fontface = face, col = col))
text_bl <- function(s, x, y, size = 7, face = "plain", just = "left")      # text on the baseline y
  grid.text(s, x = unit(x, "in"), y = unit(H - y, "in"), just = c(just, "bottom"), gp = gpar(fontfamily = "Arial", fontsize = size, fontface = face, col = "black"))
draw_at <- function(g, x, y, w, h) { pushViewport(viewport(x = x / W, y = 1 - y / H, width = w / W, height = h / H, just = c("left", "top")))
  grid.draw(g); popViewport() }

# ---- ggplot layout ----
align_h <- function(ps) { gs <- lapply(ps, ggplotGrob); h <- do.call(unit.pmax, lapply(gs, function(g) g$heights)); lapply(gs, function(g) { g$heights <- h; g }) }
nonpanel_w <- function(g) { w <- g$widths; sum(convertWidth(w[unitType(w) != "null"], "in", valueOnly = TRUE)) }
nonpanel_h <- function(g) { h <- g$heights; sum(convertHeight(h[unitType(h) != "null"], "in", valueOnly = TRUE)) }
side_w <- function(g, side) { pc <- g$layout$l[g$layout$name == "panel"]                 # width left or right of a single panel
  k <- if (side == "left") seq_len(pc - 1L) else seq.int(pc + 1L, length(g$widths)); sum(convertWidth(g$widths[k], "in", valueOnly = TRUE)) }
panels_in <- function(g, x, y, w, h) {                     # every panel of gtable g drawn at (x, y, w, h), in page inches
  L <- g$layout[grepl("^panel", g$layout$name), ]; L <- L[order(L$t, L$l), ]
  nw <- unitType(g$widths) == "null"; nh <- unitType(g$heights) == "null"
  free_w <- w - sum(convertWidth(g$widths[!nw], "in", valueOnly = TRUE)); free_h <- h - sum(convertHeight(g$heights[!nh], "in", valueOnly = TRUE))
  span <- function(u, idx, isnull, free, conv) { if (!length(idx)) return(0); uu <- u[idx]; nn <- isnull[idx]
    (if (any(!nn)) sum(conv(uu[!nn], "in", valueOnly = TRUE)) else 0) + (if (any(nn)) sum(as.numeric(uu[nn])) * free / sum(as.numeric(u[isnull])) else 0) }
  lapply(seq_len(nrow(L)), function(i) list(name = L$name[i], t = L$t[i], l = L$l[i],
    left_in = x + span(g$widths, seq_len(L$l[i] - 1L), nw, free_w, convertWidth), top_in = y + span(g$heights, seq_len(L$t[i] - 1L), nh, free_h, convertHeight),
    width_in = as.numeric(g$widths[L$l[i]]) * free_w / sum(as.numeric(g$widths[nw])), height_in = as.numeric(g$heights[L$t[i]]) * free_h / sum(as.numeric(g$heights[nh])))) }

# ---- the eleven rows of A and C-F: eight grafts, then three longer in vitro cultures ----
ROWS <- data.frame(
  id = c("A-HUES8", "A-WS4", "P-TXP", "B-0-M1", "B-0-M3", "B-3-M6-CTRL", "B-3-M6-STZ", "H-10w", "P-w2-m6", "P-w2-m12", "B-w0-w6"),
  study = c("Augsornworawat 2020, HUES8", "Augsornworawat 2020, WS4corr", "Augsornworawat 2023, HUES8 (nuclei)", "Balboa 2022, H1", "Balboa 2022, H1",
            "Balboa 2022, H1", "Balboa 2022, H1", "Du 2022, hCiPSC", "Augsornworawat 2023, HUES8 (nuclei)", "Augsornworawat 2023, HUES8 (nuclei)",
            "Balboa 2022, H1"),
  change = c("stage 6 → 6-mo graft", "stage 6 → 6-mo graft", "stage 6 wk 2 → 6-mo graft", "stage 7 wk 0 → 1-mo graft", "stage 7 wk 0 → 3-mo grafts, pooled",
             "stage 7 wk 3 → 6-mo graft, control host", "stage 7 wk 3 → 6-mo graft, STZ host", "stage 6 (vs day 2) → 10-wk graft",
             "stage 6 wk 2 → mo 6", "stage 6 wk 2 → mo 12", "stage 7 wk 0 → wk 6"),
  short = c("HUES8, 6-mo graft", "WS4corr, 6-mo graft", "HUES8 (nuclei), 6-mo graft", "H1, 1-mo graft", "H1, 3-mo grafts (pooled)", "H1, 6-mo graft (control host)",
            "H1, 6-mo graft (STZ host)", "hCiPSC, 10-wk graft", "HUES8 (nuclei), mo 6", "HUES8 (nuclei), mo 12", "H1, wk 6"),
  y = c(1:8, 9.75, 10.75, 11.75))                               # row positions: a gap between grafts and longer culture
GRAFT <- ROWS$id[1:8]; CULT <- ROWS$id[9:11]
YR <- c(max(ROWS$y) + .6, .4); BAND_ROWS <- ROWS$y[seq(1, nrow(ROWS), by = 2)]
ROW_UNIT_IN <- 1.7958 / 11.95                                  # inches per row unit in C-F (plots 1.80 in tall over 11.95 row units)
rowy <- function(id) ROWS$y[match(id, ROWS$id)]
sep <- function() geom_hline(yintercept = 9.0, color = "#9A9A9A", linetype = "dotted", linewidth = .3)   # grafts | longer culture
yrows <- function(lab) scale_y_reverse(limits = YR, breaks = if (lab) ROWS$y else NULL, labels = if (lab) ROWS$short else NULL, expand = c(0, 0))
# Each row's start and end samples: those of the four-marker rows of rows_called.tsv.
RC <- read_table("rows_called.tsv")
ARM <- RC[RC$rule == "FOUR_MARKER_UNGATED" & RC$population == "SC_SIC", ]
REF <- setNames(strsplit(ARM$ref_samples, ";"), ARM$id); TEST <- setNames(strsplit(ARM$test_samples, ";"), ARM$id)

# The row-label column of C and E has a fixed width, that of the string "HUES8 *, wk 2 → mo 12 in vitro" at 7 pt; the
# labels are right-aligned in it, so C's and E's first panels start at the same x.
Y_PAD_IN <- str_in("HUES8 *, wk 2 → mo 12 in vitro", 7) - max(vapply(ROWS$short, str_in, 0, size = 7))
M0 <- calc_element("axis.text.y.left", th())$margin
Y_PAD <- theme(axis.text.y.left = element_text(margin = margin(t = as.numeric(M0[1]), r = as.numeric(M0[2]), b = as.numeric(M0[3]), l = as.numeric(M0[4]) + Y_PAD_IN * 72.27)))

# Grey bands on alternate rows behind the panels in boxes (from panels_in) and, with label_x, behind the row labels from
# label_x to .025 in before the first panel.
band_rects <- function(boxes, label_x = NULL) {
  b <- boxes[[1]]; yy <- function(v) b$top_in + (v - YR[2]) / (YR[1] - YR[2]) * b$height_in
  left <- vapply(boxes, `[[`, 0, "left_in"); right <- vapply(boxes, function(p) p$left_in + p$width_in, 0)
  if (!is.null(label_x)) { left <- c(label_x, left); right <- c(b$left_in - .025, right) }
  for (r in BAND_ROWS) for (i in seq_along(left)) grid.rect(x = unit(left[i], "in"), y = unit(H - yy(r - .5), "in"), width = unit(right[i] - left[i], "in"),
    height = unit(yy(r + .5) - yy(r - .5), "in"), just = c("left", "top"), gp = gpar(fill = "black", alpha = BAND_ALPHA, col = NA)) }
# 'Grafts (in mice)' (on the baseline grafts_base) and 'Continued culture (in vitro)' (centered on the gap between the two
# groups of rows), flush left at the page margin as in Figure 4A; b1 is gtable g's first panel.
section_labels <- function(g, b1, grafts_base) {
  xl <- .10
  yv <- function(v) b1$top_in + (v - YR[2]) / (YR[1] - YR[2]) * b1$height_in
  g_top <- yv(max(rowy(GRAFT)) + .5); g_bot <- yv(min(rowy(CULT)) - .5); cs <- "Continued culture (in vitro)"
  text_bl("Grafts (in mice)", xl, grafts_base, face = "bold")
  text_bl(cs, xl, (g_top + g_bot) / 2 + (asc_in(cs, 7, "bold") - dsc_in(cs, 7, "bold")) / 2, face = "bold") }

# ---- marks and key glyphs ----
small_head <- function(z, col, ycol) {                    # a change too short for an arrow: a bare head, tip at the end value
  if (!nrow(z)) return(NULL)
  z$htip <- z$to; z$hbase <- z$to - z$stub * sign(z$to - z$from)
  geom_segment(data = z, aes(x = hbase, xend = htip, y = .data[[ycol]], yend = .data[[ycol]]), color = col, linewidth = .45, arrow = ARW) }
# E, F: an arrow from the start value (open circle, drawn on top) to the end value: the shaft, then the head on a stub;
# thr (shortest arrow) and stub are per row, in the plot's data units.
arrows_ef <- function(d, col) {
  d$tip0 <- d$to - d$stub * sign(d$to - d$from)
  a <- d[abs(d$to - d$from) >= d$thr, ]; z <- d[abs(d$to - d$from) < d$thr, ]
  list(geom_blank(data = d, aes(x = from, y = y)), geom_blank(data = d, aes(x = to, y = y)),
       if (nrow(a)) geom_segment(data = a, aes(x = from, xend = to, y = y, yend = y), color = col, linewidth = .45),
       if (nrow(a)) geom_segment(data = a, aes(x = tip0, xend = to, y = y, yend = y), color = col, linewidth = .45, arrow = ARW),
       if (nrow(a)) geom_point(data = a, aes(x = from, y = y), shape = 21, fill = "white", color = col, size = PT_SIZE, stroke = PT_STROKE),
       small_head(z, col, "y")) }
# Key glyphs: y is the key line's top; glyphs are centered .045 in below it, text tops .005 in above it.
key_line <- function(x, y, col, lwd = 1.3, lty = 1, len = .2, butt = FALSE)
  grid.lines(unit(c(x, x + len), "in"), unit(rep(H - y - .045, 2), "in"), gp = if (butt) gpar(col = col, fill = col, lwd = lwd, lty = lty, lineend = "butt") else gpar(col = col, lwd = lwd, lty = lty))
population_key <- function(x, y) for (k in c("SC_SIC", "SC_BETA", "SC_ALPHA")) {   # colored line and name per population
  key_line(x, y, POP_COL[[k]]); text_at(POP_LAB[[k]], x + .26, y - .005); x <- x + str_in(POP_LAB[[k]]) + .26 + .16 }

# ---- C and D: one row of plots; C's four plots end .25 in (panel to panel) before D ----
D_X <- 5.80; RIGHT <- 7.44                                     # D's plot spans D_X to RIGHT
CD_TOP <- 4.29 - 85 / 600                                      # top of the plot panels ...
CD_GROW <- 85 / 600                                            # ... which reach this much above the row's laid-out top
CD_SUB_BASE <- 4.22 - 85 / 600                                 # baseline of the column headers' last line
XLD <- c(-6.4, 1.8)                                            # D's change axis (log2)
# A C-D plot without marks: the eleven rows (row labels only with lab_y), the grafts | longer-culture rule, optionally a
# grey line at xline; d_frame() is D's (its layout is the same once the marks are added).
colplot <- function(xlim, brk, lab, lab_y = FALSE, xline = NULL, minor = waiver()) ggplot() + sep() +
  (if (!is.null(xline)) geom_vline(xintercept = xline, color = "#777777", linewidth = .3)) +
  scale_x_continuous(limits = xlim, breaks = brk, labels = lab, minor_breaks = minor, expand = c(0, 0)) + yrows(lab_y) + labs(x = NULL, y = NULL, title = NULL) + th()
d_frame <- function() colplot(XLD, c(-6, -4, -2, 0), c("-6", "-4", "-2", "0"), xline = 0) + theme(plot.margin = margin(2, 4, 0, 2))
# Top and height of a row of plots (gtable g, after align_h) whose panels start at panel_top and hold the eleven rows,
# grown by grow inches.
row_place <- function(g, panel_top, grow) {
  h <- convertHeight(g$heights, "in", valueOnly = TRUE); pr <- unique(g$layout$t[grepl("^panel", g$layout$name)])
  c(y = panel_top - sum(h[seq_len(pr - 1L)]), h = nonpanel_h(g) + (YR[1] - YR[2]) * ROW_UNIT_IN + grow) }
# x-axis titles of C and D: their tops this far under the plots, as ggplot sets one
CD_XTITLE_GAP <- convertHeight(calc_element("axis.title.x.bottom", th())$margin[1], "in", valueOnly = TRUE)
# A column header of one or more lines, centered over panel box b, its last line on CD_SUB_BASE.
col_header <- function(lines, b) { n <- length(lines)
  for (j in seq_len(n)) text_bl(lines[j], b$left_in + b$width_in / 2, CD_SUB_BASE - (n - j) * (7.5 / 72), size = 7, just = "center") }
# An x-axis title centered under the panel boxes from b0 to b1, its top at y_ax.
cd_xtitle <- function(s, b0, b1, y_ax) { x1 <- b1$left_in + b1$width_in
  grid.text(s, x = unit((b0$left_in + x1) / 2, "in"), y = unit(H - y_ax, "in"), just = c("center", "top"),
            gp = gpar(fontfamily = "Arial", fontsize = 7, col = "black")) }

# ---- E and F: one row of plots; E's three gene-group facets and F share the row heights ----
E_BOX <- c(x = .06, w = 5.34); F_BOX <- c(x = 5.55, w = 1.89)  # the two plots' places
EF_TOP <- 7.43 - 119 / 600                                     # top of the facet strips ...
EF_GROW <- 119 / 600                                           # ... which reach this much above the row's laid-out top
EF_KEY1 <- 9.67; EF_KEY2 <- 9.84                               # tops of the E-F key lines
# E: the median difference from native intestinal EC of each population on the three gene groups, at the start (from) and end
# (to) of each row, and behind each arrow the 25th-75th percentile of the group's genes at the end point (each gene's
# mean over the end samples; the arrow's end is the mean of per-sample medians, so it need not sit mid-bar).
# F: the mean fetal–adult score of the row's start and end samples, and the ticks of Figure 3H's seven intestinal rows on the same genes.
# ef_layout() returns both plots with their marks and their places on the page.
GLAB <- c(G1 = "G1: EC genes low in SC-SIC", G2 = "G2: EC genes SC-SIC have", G3 = "G3: SC-SIC's own program")
GLAB2 <- c(G1 = "G1: EC genes\nlow in SC-SIC", G2 = "G2: EC genes\nSC-SIC have", G3 = "G3: SC-SIC's\nown program")
OFF3 <- c(SC_SIC = -.31, SC_BETA = 0, SC_ALPHA = .31)          # each population's offset within its row
TICK_COL <- c(fetal = "#E8601C", adult = "#43A5B6"); TICK_LEN <- c(chose = .08, outside = .045)
ef_layout <- function() {
  EN <- read_table("E_group_gaps.tsv")
  PE <- expand.grid(population = names(POP_COL), id = ROWS$id, group = names(GLAB), stringsAsFactors = FALSE)[, c("group", "id", "population")]
  k <- match(paste(PE$group, PE$id, PE$population), paste(EN$group, EN$id, EN$population)); PE$from <- EN$from[k]; PE$to <- EN$to[k]
  PE$y <- rowy(PE$id) + OFF3[PE$population]; PE$panel <- factor(GLAB[PE$group], GLAB)
  ES <- read_table("E_end_spread.tsv"); ES$y <- rowy(ES$id) + OFF3[ES$population]; ES$panel <- factor(GLAB[ES$group], GLAB)
  stopifnot(nrow(ES) == sum(!is.na(PE$to)), all(paste(ES$group, ES$id, ES$population) %in% paste(PE$group, PE$id, PE$population)[!is.na(PE$to)]))
  # F: Figure 3H's fetal–adult score
  AG <- read_table("age_score_by_sample.tsv")
  arm_mean <- function(s, p) { t <- AG[AG$population == p, ]; mean(t$age_score[match(s, t$sample)]) }
  PA <- data.frame(id = rep(ROWS$id, each = 2), population = c("SC_SIC", "SC_BETA"))
  PA$from <- mapply(function(id, p) arm_mean(REF[[id]], p), PA$id, PA$population, USE.NAMES = FALSE)
  PA$to <- mapply(function(id, p) arm_mean(TEST[[id]], p), PA$id, PA$population, USE.NAMES = FALSE)
  PA$y <- rowy(PA$id) + OFF3[PA$population]; PA$panel <- "Fetal-to-adult axis"
  # ticks: long for the four intestinal arms (EC by primary study labels), which chose the genes; short for the three intestinal test
  # sources (cells the four-marker classifier classifies as EC), not used to choose them
  TICKS <- read_table("fig3H_ticks.tsv"); TICKS$chose <- TICKS$kind == "intestinal arm"
  F_STRIP <- sprintf("%d genes:\n%d adult-high, %d fetal-high", AG$adult_genes[1] + AG$fetal_genes[1], AG$adult_genes[1], AG$fetal_genes[1])
  XA <- c(floor(min(TICKS$score, PA$from, PA$to)) - .3, ceiling(max(TICKS$score, PA$from, PA$to)) + .3)
  # the plots; marks = FALSE: the same plot without its marks, to measure the layout
  e_plot <- function(marks) ggplot() + sep() + geom_vline(xintercept = 0, color = "#555555", linetype = "dashed", linewidth = .4) +
    (if (marks) list(geom_segment(data = ES, aes(x = q25, xend = q75, y = y, yend = y, color = population), linewidth = 1.6, alpha = .3, lineend = "butt"),
                     scale_color_manual(values = POP_COL))
     else list(geom_blank(data = ES, aes(x = q25, y = y)), geom_blank(data = ES, aes(x = q75, y = y)))) +
    lapply(names(POP_COL), function(k) { d <- PE[PE$population == k, ]
      if (marks) arrows_ef(d, POP_COL[[k]]) else list(geom_blank(data = d, aes(x = from, y = y)), geom_blank(data = d, aes(x = to, y = y))) }) +
    facet_grid(. ~ panel, scales = "free_x", space = "free_x", labeller = as_labeller(setNames(GLAB2, GLAB))) + yrows(TRUE) +
    labs(x = "Median expression difference from native intestinal EC (log2)", y = NULL) + th() + theme(panel.spacing.x = unit(10, "pt")) + Y_PAD +
    theme(plot.margin = margin(2, 4, 0, 2))
  f_plot <- function(marks) ggplot() +
    unlist(lapply(names(TICK_COL), function(g) lapply(c(TRUE, FALSE), function(ch) geom_rug(data = TICKS[TICKS$group == g & TICKS$chose == ch, ],
      aes(x = score), sides = "b", color = TICK_COL[[g]], linewidth = .6, length = unit(TICK_LEN[[if (ch) "chose" else "outside"]], "in")))), recursive = FALSE) +
    sep() + lapply(c("SC_SIC", "SC_BETA"), function(k) { d <- PA[PA$population == k, ]
      if (marks) arrows_ef(d, POP_COL[[k]]) else list(geom_blank(data = d, aes(x = from, y = y)), geom_blank(data = d, aes(x = to, y = y))) }) +
    facet_wrap(~panel, labeller = as_labeller(setNames(F_STRIP, "Fetal-to-adult axis"))) + yrows(FALSE) +
    scale_x_continuous(limits = XA, breaks = seq(-2, 6, 2), expand = c(0, 0)) +
    labs(x = "Fetal–adult score (own expression;\nadult-high minus fetal-high genes, log2)", y = NULL) + th() + theme(plot.margin = margin(2, 4, 0, 2))
  # Places: the facet strips start at EF_TOP; E's facet widths follow their x ranges. An arrow shorter than its head
  # plus 2 px is drawn as a bare head (thr, stub: in each facet's data units).
  pE0 <- e_plot(FALSE); pF0 <- f_plot(FALSE)
  T3G <- align_h(list(pE0, pF0)); h3 <- convertHeight(T3G[[1]]$heights, "in", valueOnly = TRUE); L3 <- T3G[[1]]$layout
  pr3 <- unique(L3$t[grepl("^panel", L3$name)]); st3 <- unique(L3$t[grepl("^strip-t", L3$name)])
  y <- EF_TOP - sum(h3[seq_len(pr3 - 1L)])
  h <- nonpanel_h(T3G[[1]]) + (YR[1] - YR[2]) * ROW_UNIT_IN + EF_GROW
  strip_base <- y + sum(h3[seq_len(st3 - 1L)]) + h3[st3] / 2 - .5 / 72.27 + (asc_in("G1: EC genes", 7.5, "bold") + 7.5 * .9 * 1.2 / 72) / 2
  gE0 <- ggplotGrob(pE0); nnE <- unitType(gE0$widths) == "null"; nvE <- as.numeric(gE0$widths[nnE])
  bE0 <- ggplot_build(pE0); LE <- bE0$layout$layout
  xrE <- vapply(GLAB, function(lab) bE0$layout$panel_params[[LE$PANEL[LE$panel == lab]]]$x.range, c(0, 0))
  ipu_E <- (E_BOX[["w"]] - sum(convertWidth(gE0$widths[!nnE], "in", valueOnly = TRUE))) * nvE / sum(nvE) / (xrE[2, ] - xrE[1, ])
  gF0 <- ggplotGrob(pF0); nnF <- unitType(gF0$widths) == "null"; xrF <- ggplot_build(pF0)$layout$panel_params[[1]]$x.range
  ipu_F <- (F_BOX[["w"]] - sum(convertWidth(gF0$widths[!nnF], "in", valueOnly = TRUE))) / (xrF[2] - xrF[1])
  k <- match(GLAB[PE$group], GLAB); PE$thr <- (HEAD_MIN_IN / ipu_E)[k]; PE$stub <- (STUB_IN / ipu_E)[k]
  PA$thr <- HEAD_MIN_IN / ipu_F; PA$stub <- STUB_IN / ipu_F
  t3 <- align_h(lapply(list(e_plot(TRUE), f_plot(TRUE)), nobg))
  list(gE = t3[[1]], gF = t3[[2]], y = y, h = h, strip_base = strip_base,
       bxE = panels_in(t3[[1]], E_BOX[["x"]], y, E_BOX[["w"]], h), bxF = panels_in(t3[[2]], F_BOX[["x"]], y, F_BOX[["w"]], h)) }
