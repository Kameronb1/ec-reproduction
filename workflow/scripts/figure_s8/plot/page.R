# The Supplementary Figure S8 page, shared by its panel scripts: the tables, the measured layout and one drawing function
# per panel. Each panel script sources this file and draws its own part.
#   A  the 26 Reactome cholesterol-synthesis genes: change from the earlier in vitro samples in 8 grafts and 5 comparisons
#      kept in culture, per cell type (tiles; colors end at +-3 log2; a dot where the mean change is at least 0.5 log2 with
#      the same sign in every earlier and later sample pair), and each row's median change (bars and values). Balboa 2022
#      week 0 -> week 3 is in the table but not drawn.
#   B  promoter openness of G1-G3 against the typical promoter of the same cell type, Hickey 2023 primary intestinal EC
#      and other enteroendocrine cells, by region; Figure 4D's fold scale and y range
#   C  promoter change after grafting (Augsornworawat 2023, week 2 -> graft) against 10,000 random gene sets; Figure 4E's
#      percent scale (4E draws the distal rows of the same table)
# Positions are in inches from the page's top-left corner. The layout is measured on a cairo device of the page's size.

text_arrow <- function(s, x, y, just = "left") {             # text on the baseline y, each arrow drawn bold
  parts <- strsplit(s, "→", fixed = TRUE)[[1]]; if (endsWith(s, "→")) parts <- c(parts, "")
  runs <- c(rbind(parts, c(rep("→", length(parts) - 1L), "")))[seq_len(2L * length(parts) - 1L)]
  faces <- rep(c("plain", "bold"), length.out = length(runs)); w <- mapply(str_in, runs, 7, faces)
  if (just == "right") x <- x - sum(w)
  for (k in seq_along(runs)) { x0 <- x + sum(w[seq_len(k - 1L)])
    if (runs[k] == "→") grid.lines(x = unit(x0 + w[k] * c(.10, .90), "in"), y = unit(rep(H - y + .034, 2), "in"),
                                     arrow = arrow(angle = 25, length = unit(.032, "in"), type = "closed"),
                                     gp = gpar(col = "black", fill = "black", lwd = 1.1, linejoin = "mitre"))
    else if (nzchar(runs[k])) text_bl(runs[k], x0, y) } }
rect_in <- function(x0, x1, y0, y1, fill, col = NA, lwd = .5, alpha = 1)   # page inches, y down from the top
  grid.rect(unit(x0, "in"), unit(H - y0, "in"), width = unit(x1 - x0, "in"), height = unit(y1 - y0, "in"), just = c("left", "top"),
            gp = gpar(fill = fill, col = col, lwd = lwd, alpha = alpha, linejoin = "mitre"))
seg_in <- function(x0, y0, x1, y1, col = "black", lwd = .5, lty = 1)
  grid.segments(unit(x0, "in"), unit(H - y0, "in"), unit(x1, "in"), unit(H - y1, "in"), gp = gpar(col = col, lwd = lwd, lty = lty, lineend = "butt"))
minus <- function(s) gsub("-", "−", s, fixed = TRUE)

# ---- text -----------------------------------------------------------------------------------------------------------------
RIGHT <- 7.44; SEC_PT <- 7.5
SECTION <- c(graft = "Grafts (in mice)", cult = "Kept in culture (in vitro)")
TAB_HEAD <- c(study = "Study", line = "Line and end point", n1 = "Samples", n2 = "start → end")
A_HEAD <- c(genes = "26 cholesterol-synthesis genes (Reactome), smallest fall in grafts at left", med1 = "Median change", med2 = "(log2)")
A_KEY <- c(scale = "Change from the earlier in vitro samples, log2(1 + CPM)",
           dot1 = "Mean change of at least 0.5 log2, with the same sign", dot2 = "in every earlier and later sample pair")
GL <- c("G1: EC genes low in SC-SICs", "G2: shared EC genes", "G3: separate SC-SIC program")
BC_KEY <- c(median = "Median gene of the group", bars = "Middle half of the group genes", random = "Middle 95% of 10,000 random gene sets")
B_FACET <- c(EC = "Primary intestinal EC", `other enteroendocrine` = "Other enteroendocrine")
B_YLAB <- "Openness relative\nto the typical gene\npromoter of the same\ncell type (fold)"      # Figure 4D's y title
C_YLAB <- "% more opening\nthan the promoter\nof the typical gene"                          # Figure 4E's, for promoters
C_XLAB <- "Graft (week 2 → graft)"

# ---- A data ---------------------------------------------------------------------------------------------------------------
TT <- read_table("s2_cholesterol_tiles.tsv"); MED <- read_table("s2_cholesterol_medians.tsv")
LABS <- read.delim(file.path(dirname(PLOT_SCRIPT), "s8_rows.tsv"), check.names = FALSE, stringsAsFactors = FALSE, encoding = "UTF-8")
DROP <- "B-w0-w3"; stopifnot(DROP %in% TT$comparison, DROP %in% MED$comparison)
TT <- TT[TT$comparison != DROP, ]; MED <- MED[MED$comparison != DROP, ]; LABS <- LABS[LABS$id != DROP, ]
# rows: graft rows in Figure 4B's order; kept-in-culture rows in the shared study order (Balboa 2022, Augsornworawat
# 2023, Zhu 2023), then Veres 2019
ROWS <- data.frame(
  id = c("A-HUES8", "A-WS4", "B-0-M1", "B-0-M3", "B-3-M6-CTRL", "B-3-M6-STZ", "H-10w", "P-TXP",
         "B-w0-w6", "P-w2-m6", "P-w2-m12", "Z-D21-D39", "V-w0-w5"),
  block = rep(c("graft", "cult"), c(8, 5)),
  study = c("Augsornworawat 2020", "Augsornworawat 2020", rep("Balboa 2022", 4), "Du 2022", "Augsornworawat 2023",
            "Balboa 2022", "Augsornworawat 2023", "Augsornworawat 2023", "Zhu 2023", "Veres 2019"),
  short = c("HUES8, 6-mo graft", "WS4corr, 6-mo graft", "H1, 1-mo graft", "H1, 3-mo grafts (pooled)", "H1, 6-mo graft (control host)",
            "H1, 6-mo graft (STZ host)", "hCiPSC, 10-wk graft", "HUES8 (nuclei), 6-mo graft", "H1, wk 6", "HUES8 (nuclei), mo 6",
            "HUES8 (nuclei), mo 12", "H1, day 39", "HUES8, wk 5"), stringsAsFactors = FALSE)
stopifnot(setequal(ROWS$id, LABS$id), setequal(ROWS$id, TT$comparison), setequal(ROWS$id, MED$comparison))
ROWS$samples <- LABS$samples[match(ROWS$id, LABS$id)]                    # "start -> end" sample counts
POPS <- c("SC_SIC", "SC_BETA", "SC_ALPHA")
GO <- unique(TT$gene[order(TT$gene_column)]); NG <- length(GO); stopifnot(NG == 26)
PAL <- scales::pal_div_gradient("#2A7F8E", "grey97", "#B4553F", "Lab")      # +-3 log2 at the ends
fill_of <- function(v) PAL(scales::rescale_mid(pmin(pmax(v, -3), 3), to = c(0, 1), from = c(-3, 3), mid = 0))

# ---- B, C data ---------------------------------------------------------------------------------------------------------------
GCOL5 <- c(G1 = "#43A5B6", G2 = "#A89A2E", G3 = "#B4553F")                   # Figure 4D/E gene-group colours
GCOL5_LIGHT <- vapply(GCOL5, function(c) grDevices::rgb(t(.4 * grDevices::col2rgb(c) + .6 * 255), maxColorValue = 255), "")
BAND5 <- "#D9D9D9"; pct5 <- function(d) 100 * (2^d - 1)
hk <- read_table("c5_relative_zhu_hickey.tsv", "figure4"); hk <- hk[hk$study == "Hickey 2023", ]
RX <- c(`small intestine` = 1, colon = 2.15); GDX <- c(G1 = -.2, G2 = 0, G3 = .2)
hk$P <- factor(B_FACET[hk$population], B_FACET); hk$x <- unname(RX[hk$region] + GDX[hk$g])
D_YLIM <- c(-3.2, 1.75)                                                       # Figure 4D's y range
stopifnot(nrow(hk) == 12, !anyNA(hk$P), !anyNA(hk$x), min(hk$q25) > D_YLIM[1], max(hk$q75) < D_YLIM[2])
gr <- read_table("c5_graft_relative_by_population.tsv", "figure4"); gr <- gr[gr$window == "promoter", ]
gr$y <- pct5(gr$group_rel); gr$lo <- pct5(gr$null_q025_rel); gr$hi <- pct5(gr$null_q975_rel)
gr$u <- factor(c(txp_1 = "1", txp_2 = "2", txp_3 = "3")[gr$graft], c("1", "2", "3")); gr$win <- factor(gr$population, c("SC-SIC", "SC-beta", "SC-alpha"))
stopifnot(nrow(gr) == 27, !anyNA(gr$u), !anyNA(gr$win))

th5 <- function() theme_minimal(base_family = "Arial", base_size = 7.5) + theme(
  text = element_text(color = "black", size = 7.5), axis.text = element_text(color = "black", size = 7),
  axis.title = element_text(color = "black", size = 7), strip.text = element_text(color = "black", size = 7, face = "bold", margin = margin(1, 0, 2, 0)),
  panel.grid.minor = element_blank(), panel.grid.major.x = element_blank(), panel.grid.major.y = element_line(color = "#E4E4E4", linewidth = .2),
  axis.line.x = element_line(color = "black", linewidth = .35), axis.ticks.x = element_line(color = "black", linewidth = .3),
  axis.ticks.length = unit(1.3, "pt"), plot.margin = margin(2, 4, 2, 2), legend.position = "none", plot.background = element_blank())
zero_h5 <- geom_hline(yintercept = 0, color = "#777777", linewidth = .3)
pB <- ggplot(hk, aes(x, rel, color = g)) + zero_h5 +
  geom_rect(aes(xmin = x - .09, xmax = x + .09, ymin = q25, ymax = q75), fill = unname(GCOL5_LIGHT[hk$g]), colour = NA) +
  geom_line(aes(group = g), linewidth = .45) + geom_point(aes(fill = g), shape = 21, colour = "black", stroke = .15, size = 1.3) +
  facet_wrap(~P, nrow = 1) + scale_color_manual(values = GCOL5) + scale_fill_manual(values = GCOL5) +
  scale_x_continuous(breaks = unname(RX), labels = c("Small\nintestine", "Colon"), limits = c(.5, 2.65)) +
  scale_y_continuous(breaks = log2(c(.125, .25, .5, 1, 2)), labels = c("×0.125", "×0.25", "×0.5", "×1", "×2")) + coord_cartesian(ylim = D_YLIM) +
  labs(x = quote("Region (" * bold("Hickey 2023") * ")"), y = B_YLAB) + th5() +
  theme(panel.spacing = unit(8, "pt"), axis.text.x = element_text(size = 7, lineheight = .9))
pC <- ggplot(gr, aes(u, y, group = g)) + zero_h5 +
  geom_linerange(aes(ymin = lo, ymax = hi), position = position_dodge(.8), linewidth = .85, color = BAND5) +
  geom_point(aes(fill = g), shape = 21, colour = "black", stroke = .15, size = 1.2, position = position_dodge(.8)) + facet_wrap(~win, nrow = 1) +
  scale_fill_manual(values = GCOL5) + labs(x = C_XLAB, y = C_YLAB) + th5() +
  theme(panel.spacing = unit(5, "pt"), axis.line.x = element_blank(), axis.ticks.x = element_blank(), axis.text.x = element_text(size = 7))

# ==== layout (inches from the page's top-left), measured on a cairo device ====================================================
MEASURE <- tempfile(fileext = ".png"); open_png(MEASURE)
# A: label table, cell-type chips, heatmap, median bars and their values
TB_X_STUDY <- .14
STUDY_W <- max(vapply(unique(ROWS$study), str_in, 0, size = 7, face = "bold"), str_in(TAB_HEAD[["study"]], 7))
TB_X_LINE <- TB_X_STUDY + STUDY_W + .12
LINE_W <- max(vapply(ROWS$short, str_in, 0, size = 7), str_in(TAB_HEAD[["line"]], 7))
N_W <- max(str_in(TAB_HEAD[["n1"]], 7), str_in("start  → end", 7), vapply(ROWS$samples, str_in, 0, size = 7))
TB_N_RIGHT <- TB_X_LINE + LINE_W + .14 + N_W
CHIP0 <- TB_N_RIGHT + .09; CHIP1 <- CHIP0 + .075; HM0 <- CHIP1 + .035
VAL_W <- max(vapply(minus(MED$lab), str_in, 0, size = 7)); BAR1 <- RIGHT - VAL_W - .06; BAR0 <- BAR1 - .62
HM1 <- BAR0 - .14; GW <- (HM1 - HM0) / NG
BAR_LIM <- c(-2.4, .3); bx <- function(v) BAR0 + (v - BAR_LIM[1]) / diff(BAR_LIM) * (BAR1 - BAR0)
stopifnot(GW >= .095, min(MED$change) > BAR_LIM[1], max(MED$change) < BAR_LIM[2])
GENE_H <- max(vapply(GO, str_in, 0, size = 7, face = "italic"))
LT_A <- .10
HEAD_BASE <- .22                                                   # the gene-header line, beside the letter
GENE_BOT <- HEAD_BASE + .08 + GENE_H                               # rotated gene names end here
TAB_BASE <- GENE_BOT                                               # table heads on this baseline
RH <- .09; CH <- 3 * RH; CG <- .034                               # a cell-type row, a comparison, the gap between comparisons
SEC_GAP <- .20                                                     # section header line height
Y0 <- TAB_BASE + .07 + SEC_GAP                                     # first graft comparison's top
ROWS$top <- NA_real_; y <- Y0
for (i in seq_len(nrow(ROWS))) { if (i > 1 && ROWS$block[i] != ROWS$block[i - 1]) y <- y + SEC_GAP + .02; ROWS$top[i] <- y; y <- y + CH + CG }
A_BOT <- max(ROWS$top) + CH
SEC_BASE <- c(graft = Y0 - CG / 2 - .055, cult = ROWS$top[ROWS$block == "cult"][1] - CG / 2 - .055)
row_top <- function(id, pop) ROWS$top[match(id, ROWS$id)] + (match(pop, POPS) - 1) * RH
Y_AX <- A_BOT + .03
KY <- Y_AX + .25                                                    # A's key line (top)
CB0 <- 2.55; CB1 <- CB0 + 1.85; DOT_X <- CB1 + .32
stopifnot(DOT_X + .08 + max(str_in(A_KEY[["dot1"]]), str_in(A_KEY[["dot2"]])) < RIGHT)
A_END <- KY + .36
# B | C: one row, one key under both
GAP_ROWS <- .14
LT_B <- A_END + GAP_ROWS
BC_Y <- LT_B + .067; BC_H <- 1.95
B_X <- .06; B_W <- 3.05; C_X <- B_X + B_W + .14; C_W <- RIGHT - C_X
T2 <- align_h(list(pB, pC))
BC_KEY_TOP <- BC_Y + BC_H + .05 + c(0, .135)
PAGE_END <- BC_KEY_TOP[2] + .135
invisible(dev.off()); unlink(MEASURE)
stopifnot(PAGE_END < H - .1)

# ==== drawing ==================================================================================================================
BAND_IDS <- ROWS$id[seq(1, nrow(ROWS), by = 2)]                      # alternate comparisons shaded, unbroken from the table to the values
draw_A <- function() {
  text_at("A", .10, LT_A, 11, "bold")
  for (i in which(ROWS$id %in% BAND_IDS)) rect_in(.06, RIGHT, ROWS$top[i] - CG / 2, ROWS$top[i] + CH + CG / 2, "black", alpha = BAND_ALPHA)
  # headers
  text_bl(A_HEAD[["genes"]], (HM0 + HM1) / 2, HEAD_BASE, just = "centre")
  for (k in seq_along(GO))
    grid.text(GO[k], unit(HM0 + (k - .5) * GW, "in"), unit(H - GENE_BOT, "in"), rot = 90, hjust = 0, vjust = .5,
              gp = gpar(fontfamily = "Arial", fontsize = 7, fontface = "italic"))
  text_bl(A_HEAD[["med1"]], (BAR0 + RIGHT) / 2, TAB_BASE - 7.5 / 72, just = "centre"); text_bl(A_HEAD[["med2"]], (BAR0 + RIGHT) / 2, TAB_BASE, just = "centre")
  stopifnot((BAR0 + RIGHT) / 2 + str_in(A_HEAD[["med1"]], 7) / 2 < RIGHT + .03)
  text_bl(TAB_HEAD[["study"]], TB_X_STUDY, TAB_BASE); text_bl(TAB_HEAD[["line"]], TB_X_LINE, TAB_BASE)
  text_bl(TAB_HEAD[["n1"]], TB_N_RIGHT, TAB_BASE - 7.5 / 72, just = "right"); text_arrow(TAB_HEAD[["n2"]], TB_N_RIGHT, TAB_BASE, just = "right")
  for (k in names(SECTION)) text_bl(SECTION[[k]], .06, SEC_BASE[[k]], size = SEC_PT, face = "bold")
  # label table: a study is named on its first row in each section
  first <- c(TRUE, ROWS$study[-1] != ROWS$study[-nrow(ROWS)] | ROWS$block[-1] != ROWS$block[-nrow(ROWS)])
  for (i in seq_len(nrow(ROWS))) { b <- ROWS$top[i] + CH / 2 + asc_in("H", 7) / 2
    if (first[i]) text_bl(ROWS$study[i], TB_X_STUDY, b, face = "bold")
    text_bl(ROWS$short[i], TB_X_LINE, b); text_arrow(ROWS$samples[i], TB_N_RIGHT, b, just = "right")
    for (p in POPS) rect_in(CHIP0, CHIP1, row_top(ROWS$id[i], p), row_top(ROWS$id[i], p) + RH, POP_COL[[p]], col = "white", lwd = .3 * ggplot2::.pt) }
  # bar grid and zero line
  for (v in c(-2, -1)) seg_in(bx(v), Y0 - CG / 2, bx(v), A_BOT, col = "grey80", lwd = .5)
  # tiles (white hairline borders), dots
  x0 <- HM0 + (TT$gene_column - 1) * GW; top <- row_top(TT$comparison, TT$population)
  rect_in(x0, x0 + GW, top, top + RH, fill_of(TT$change), col = "white", lwd = .3 * ggplot2::.pt)
  d <- which(TT$dot); grid.points(unit(x0[d] + GW / 2, "in"), unit(H - top[d] - RH / 2, "in"), pch = 16, size = unit(.75 * 72.27 / 25.4 + .2, "pt"), gp = gpar(col = "grey15"))
  # median bars and values
  mt <- row_top(MED$comparison, MED$population)
  rect_in(pmin(bx(0), bx(MED$change)), pmax(bx(0), bx(MED$change)), mt + .018, mt + RH - .018, POP_COL[MED$population])
  for (j in seq_len(nrow(MED))) text_bl(minus(MED$lab[j]), RIGHT, mt[j] + RH / 2 + asc_in("0", 7) / 2, just = "right")
  seg_in(bx(0), Y0 - CG / 2, bx(0), Y_AX, lwd = .6); seg_in(bx(BAR_LIM[1]), Y_AX, bx(BAR_LIM[2]), Y_AX, lwd = .5)
  for (v in c(-2, -1, 0)) { seg_in(bx(v), Y_AX, bx(v), Y_AX + .025); text_at(minus(as.character(v)), bx(v) - str_in(minus(as.character(v)), 7) / 2, Y_AX + .04, 7) }
  # key: cell-type chips; colour scale; dot rule
  x <- .14; yc <- KY + .045
  for (p in POPS) { rect_in(x, x + .075, yc - .04, yc + .04, POP_COL[[p]]); text_at(POP_LAB[[p]], x + .12, KY - .005, 7); x <- x + .12 + str_in(POP_LAB[[p]], 7) + .18 }
  stopifnot(x < CB0 - .1)
  text_at(A_KEY[["scale"]], CB0, KY - .135, 7)
  e <- seq(CB0, CB1, length.out = 231); v <- -3 + (seq_len(230) - .5) / 230 * 6
  rect_in(e[-231], e[-1], yc - .045, yc + .045, fill_of(v), col = fill_of(v), lwd = .2)
  for (k in -3:3) { xk <- CB0 + (k + 3) / 6 * (CB1 - CB0); seg_in(xk, yc + .045, xk, yc + .07)
    s <- c("≤ −3", "−2", "−1", "0", "1", "2", "≥ 3")[k + 4]; text_at(s, xk - str_in(s, 7) / 2, yc + .08, 7) }
  grid.points(unit(DOT_X, "in"), unit(H - yc, "in"), pch = 16, size = unit(.75 * 72.27 / 25.4 + .2, "pt"), gp = gpar(col = "grey15"))
  text_at(A_KEY[["dot1"]], DOT_X + .08, KY - .005, 7); text_at(A_KEY[["dot2"]], DOT_X + .08, KY - .005 + .117, 7) }
row_keys <- function(items, left, width, y, gap = .2) {           # key row (pb, pt, bars3, vbar)
  ww <- vapply(items, function(it) str_in(it$label) + if (it$type == "bars3") .26 else .14, 0); total <- sum(ww) + gap * (length(ww) - 1)
  x <- left + (width - total) / 2; PLWD <- .15 * ggplot2::.stroke / 2
  for (i in seq_along(items)) { it <- items[[i]]; yy <- H - y - .045
    if (it$type == "pt") grid.points(unit(x + .05, "in"), unit(yy, "in"), pch = 21, size = unit(1.5, "mm"), gp = gpar(col = "black", fill = it$fill, lwd = PLWD))
    if (it$type == "vbar") grid.rect(unit(x + .05, "in"), unit(yy, "in"), width = unit(.85 * ggplot2::.pt / 96, "in"), height = unit(.09, "in"), gp = gpar(fill = BAND5, col = NA))
    if (it$type == "bars3") for (k in 1:3) grid.rect(unit(x + .025 + (k - 1) * .07, "in"), unit(yy, "in"), width = unit(.05, "in"), height = unit(.09, "in"), gp = gpar(fill = GCOL5_LIGHT[[k]], col = NA))
    if (it$type == "pb") { grid.rect(unit(x + .05, "in"), unit(yy, "in"), width = unit(.05, "in"), height = unit(.13, "in"), gp = gpar(fill = it$bar, col = NA))
      grid.points(unit(x + .05, "in"), unit(yy, "in"), pch = 21, size = unit(1.5, "mm"), gp = gpar(col = "black", fill = it$fill, lwd = PLWD)) }
    text_at(it$label, x + if (it$type == "bars3") .26 else .14, y - .005)
    x <- x + ww[i] + gap }
  invisible(c(left + (width - total) / 2, left + (width + total) / 2)) }
draw_B <- function() { text_at("B", .10, LT_B, 11, "bold"); draw_at(T2[[1]], B_X, BC_Y, B_W, BC_H) }
draw_C <- function() { text_at("C", C_X + .02, LT_B, 11, "bold"); draw_at(T2[[2]], C_X, BC_Y, C_W, BC_H) }
draw_BC_key <- function() {
  K1 <- row_keys(lapply(1:3, function(k) list(type = "pb", label = GL[k], fill = GCOL5[[k]], bar = GCOL5_LIGHT[[k]])), .10, RIGHT - .10, BC_KEY_TOP[1])
  K2 <- row_keys(list(list(type = "pt", label = BC_KEY[["median"]], fill = "grey75"), list(type = "bars3", label = BC_KEY[["bars"]]),
                      list(type = "vbar", label = BC_KEY[["random"]])), .10, RIGHT - .10, BC_KEY_TOP[2])
  stopifnot(K1[1] > 0, K1[2] < RIGHT, K2[1] > 0, K2[2] < RIGHT) }
