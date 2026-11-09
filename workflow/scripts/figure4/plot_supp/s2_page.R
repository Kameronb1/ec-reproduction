# Figure S16, shared by S16A-D (sourced by each after helpers.R): the tables, the page geometry and the four panels'
# drawing.
#   A  heatmap of the 26 cholesterol-synthesis genes of Figure 4D per comparison and cell type, each line's median change
#      as a bar, and the NFE2L1 regulon target score change (AUCell) as a second bar column on the same rows;
#   B  regulons most enriched for the shared graft decline, and
#   C  for cholesterol synthesis: one table on one shared log axis ("targets measured" is the n of "x of n");
#   D  transcription factors whose (+) regulons are decline hits in two or more networks, and NFE2L1: networks down,
#      factors across.
# Text is 7 pt Arial. Positions are inches on the page (x from the left, y down from the top: a text's baseline or a
# shape's top); the geometry is measured on the open 600-dpi cairo page device, as the panels are drawn. The page's
# panel letters are set separately and are not drawn here.

# ---- data ----------------------------------------------------------------------------------------------------------------------
TT <- read_table("s2_cholesterol_tiles.tsv")
MED <- read_text_cols("s2_cholesterol_medians.tsv", "lab")
CV <- read_table("s2_nfe2l1_bars.tsv")
# The drawn rows in order: block and two-line label (drawing-only, s2_rows.tsv). Balboa 2022 stage 7 wk 0 → wk 3 is in the
# tables but not drawn: it shares its three wk 0 samples, and the first half of its interval, with wk 0 → wk 6.
LAB <- read_drawing("s2_rows.tsv")
stopifnot(all(LAB$id %in% TT$comparison), all(LAB$id %in% MED$comparison), all(LAB$id %in% CV$comparison))
TT <- TT[TT$comparison %in% LAB$id, ]; MED <- MED[MED$comparison %in% LAB$id, ]; CV <- CV[CV$comparison %in% LAB$id, ]
# the samples at each end are the tiles' counts of start and end samples; a pooled end (one profile) counts the samples
# it pools (pooled_test)
NS <- unique(TT[c("comparison", "ref_n", "test_n")]); stopifnot(!anyDuplicated(NS$comparison), setequal(NS$comparison, LAB$id))
k <- match(LAB$id, NS$comparison); pooled <- suppressWarnings(as.integer(LAB$pooled_test)); stopifnot(all(NS$test_n[k][!is.na(pooled)] == 1L))
LAB$samples <- sprintf("%d → %d", NS$ref_n[k], ifelse(is.na(pooled), NS$test_n[k], pooled))
LAB$band <- ave(seq_along(LAB$id), LAB$block, FUN = seq_along) %% 2 == 1     # each block starts grey, then alternates
stopifnot(setequal(paste(CV$comparison, CV$population), as.vector(outer(LAB$id, POPS, paste))))
FA <- read_table("s3a_decline_screen.tsv"); FB <- read_table("s3b_cholesterol_screen.tsv")
DT <- read_text_cols("s3d_regulon_recurrence.tsv", c("tile_top", "tile_bottom"))
GO <- unique(TT$gene[order(TT$gene_column)]); NG <- length(GO)                 # genes, left to right
PAL <- scales::pal_div_gradient("#2A7F8E", "gray97", "#B4553F", "Lab")
fill_of <- function(v) PAL(scales::rescale_mid(pmin(pmax(v, -3), 3), to = c(0, 1), from = c(-3, 3), mid = 0))   # +-3 log2 at the ends
GRP_LEV <- c("NFE2L1", "SOX4, SOX11 (also fall in culture alone, where scored)", "other regulons")   # marks: red circle, grey diamond, black circle
GRP_COL <- setNames(c("#CC3311", "gray50", "black"), GRP_LEV); GRP_PCH <- setNames(c(16, 18, 16), GRP_LEV); GRP_SZ <- setNames(c(1.6, 2.15, 1.6), GRP_LEV)
grp_of <- function(tf) ifelse(tf %in% c("SOX4", "SOX11"), GRP_LEV[2], ifelse(tf == "NFE2L1", GRP_LEV[1], GRP_LEV[3]))
LINE <- 7 * 1.2 / 72                                                            # 7 pt line pitch (in)
PTMM <- 72.27 / 25.4; STROKE <- 96 / 25.4                                       # ggplot2's .pt and .stroke

# ---- drawing helpers: page inches, drawn into the current panel ----------------------------------------------------------------
# A panel's draw function sets OX, OY (its top-left on the page, in) and PH (its height, in), so every element keeps its page place.
OX <- 0; OY <- 0; PH <- 10
U <- function(x) unit(x - OX, "in"); V <- function(y) unit(PH - (y - OY), "in")
gp_t <- function(size = 7, face = "plain", col = "black") gpar(fontfamily = FONT, fontsize = size, fontface = face, col = col)
tw <- function(s, size = 7, face = "plain")                                     # widths on the device being drawn
  vapply(s, function(x) convertWidth(grobWidth(textGrob(x, gp = gp_t(size, face))), "in", valueOnly = TRUE), 0, USE.NAMES = FALSE)
tx <- function(s, x, base, size = 7, face = "plain", hjust = 0, rot = 0, col = "black")
  grid.text(s, U(x), V(base), hjust = hjust, vjust = 0, rot = rot, gp = gp_t(size, face, col))
pieces <- function(p, x, base, hjust = 0) {                                     # p: list of c(text, face), placed by measured widths
  w <- vapply(p, function(s) tw(s[1], 7, s[2]), 0); x <- x - hjust * sum(w)
  for (k in seq_along(p)) { tx(p[[k]][1], x, base, 7, p[[k]][2]); x <- x + w[k] } }
tx_arrow <- function(s, x, base, hjust = 0) {                                 # 7 pt text whose arrows are drawn, as in Figure 4A:
  parts <- strsplit(s, "→", fixed = TRUE)[[1]]; if (endsWith(s, "→")) parts <- c(parts, "")   # Arial has no arrow glyph (the
  runs <- c(rbind(parts, "→"))[seq_len(2L * length(parts) - 1L)]                # fallback font has no bold), so each is a 1.1 pt
  w <- tw(runs, 7, "bold"); w[runs != "→"] <- tw(runs[runs != "→"]); x <- x - hjust * sum(w)   # line with a closed head
  for (k in seq_along(runs)) { x0 <- x + sum(w[seq_len(k - 1L)])
    if (runs[k] == "→") grid.lines(U(x0 + w[k] * c(.10, .90)), V(rep(base - .034, 2)), arrow = arrow(angle = 25, length = unit(.032, "in"), type = "closed"),
                                     gp = gpar(col = "black", fill = "black", lwd = 1.1, linejoin = "mitre"))
    else if (nzchar(runs[k])) tx(runs[k], x0, base) } }
rct <- function(x0, x1, y0, y1, fill, col = NA, lwd = 0.5, lty = 1) grid.rect(U(x0), V(y1), width = unit(x1 - x0, "in"), height = unit(y1 - y0, "in"),
  just = c("left", "bottom"), gp = gpar(fill = fill, col = col, lwd = lwd, lty = lty, linejoin = "mitre"))
sgm <- function(x0, y0, x1, y1, col = "black", lwd = 0.5, lty = 1) grid.segments(U(x0), V(y0), U(x1), V(y1), gp = gpar(col = col, lwd = lwd, lty = lty, lineend = "butt"))
pt <- function(x, y, pch, size_mm, col) grid.points(U(x), V(y), pch = pch, size = unit(1, "char"), gp = gpar(col = col, fill = col, fontsize = size_mm * PTMM + 0.5 * STROKE / 2))
minus <- function(v) gsub("-", "−", format(v, drop0trailing = TRUE, trim = TRUE), fixed = TRUE)

# ==== geometry ===================================================================================================================
# A: heatmap columns, median bars and values, NFE2L1 bars
X <- list(lab = .13, band0 = .10, samp_r = 2.20, band1 = 2.23, chip0 = 2.27, chip1 = 2.35, hm0 = 2.39, gw = .116)
X$hm1 <- X$hm0 + NG * X$gw
X$bar0 <- X$hm1 + .12; X$bar1 <- X$bar0 + .55; X$val_r <- X$bar1 + .29
X$c0 <- X$val_r + .12; X$c1 <- 7.40
# Row labels and their bands as in Figure 4A: indented under the group headers, and a row whose study repeats the row above names only its
# change. Everything from the sample column to the NFE2L1 panel's left edge moves right as far as the widest label needs.
A_INDENT <- .10; A_CHANGE_INDENT <- .10   # change lines sit a further .10 in in, under their study
st <- paste0(LAB$bold, LAB$rest); LAB$repeat_study <- c(FALSE, st[-1] == st[-nrow(LAB)] & LAB$block[-1] == LAB$block[-nrow(LAB)])
need <- X$lab + A_INDENT + max(ifelse(LAB$repeat_study, 0, tw(LAB$bold, 7, "bold") + tw(LAB$rest)), A_CHANGE_INDENT + tw(LAB$change)) + .06 + max(tw(LAB$samples)) - X$samp_r
if (need > 0) for (k in c("samp_r", "band1", "chip0", "chip1", "hm0", "hm1", "bar0", "bar1", "val_r", "c0")) X[[k]] <- X[[k]] + need
BAR_LIM <- c(-2.3, 0.4); bx <- function(v) X$bar0 + (v - BAR_LIM[1]) / diff(BAR_LIM) * (X$bar1 - X$bar0)
# NFE2L1 axis: the drawn changes and 0, with 5% margins
LIMC <- range(c(0, CV$difference)); LIMC <- LIMC + c(-1, 1) * .05 * diff(LIMC); cx <- function(v) X$c0 + (v - LIMC[1]) / diff(LIMC) * (X$c1 - X$c0)
BRKC <- NULL                                                                    # its breaks: the first set whose labels fit
for (cand in list(c(-.06, -.03, 0, .03), c(-.05, -.025, 0, .025), c(-.06, 0, .03), c(-.05, 0), c(-.06, 0))) {
  cand <- cand[cand >= LIMC[1] & cand <= LIMC[2]]
  if (length(cand) > 1 && min(diff(cx(cand))) >= max(tw(minus(cand))) + .05) { BRKC <- cand; break } }
stopifnot(length(BRKC) > 1)
RH <- .088; CH <- 3 * RH; CG <- .024                                            # cell-type line, one comparison, gap between comparisons
BLK <- list(graft = list(head = "Grafts (in mice)", base = .935, top = .97), in_vitro = list(head = "Continued culture (in vitro)", base = 3.105, top = 3.14))
# the culture block's header sits .143 in below the last graft row, its rows .035 in below that
BLK$in_vitro$base <- BLK$graft$top + sum(LAB$block == "graft") * (CH + CG) - CG + .143; BLK$in_vitro$top <- BLK$in_vitro$base + .035
LAB$top <- NA_real_
for (b in names(BLK)) { i <- which(LAB$block == b); LAB$top[i] <- BLK[[b]]$top + (seq_along(i) - 1) * (CH + CG) }
LAB$center <- LAB$top + CH / 2; A_BOT <- max(LAB$top) + CH
row_top <- function(cmp, pop) LAB$top[match(cmp, LAB$id)] + (match(pop, POPS) - 1) * RH
Y_AX <- A_BOT + .02; KY <- Y_AX + .32
cb0 <- X$hm0 - .25; cb1 <- cb0 + 2.30
DOTL <- c("Mean rise or fall of at least 0.5 log2, and every earlier ×", "later sample pair changes the same way (needs ≥ 2", "pairs, so 1 → 1 rows carry no dot)")
stopifnot(cb1 + .38 + max(tw(DOTL)) <= 7.45)
A_END <- KY + .265 + .02
# B, C: one forest table
prep <- function(x) { x$rest <- sub(", [0-9,]+ targets\\)$", ")", substring(x$label, nchar(x$tf) + 1L)); x$grp <- grp_of(x$tf)
  x$nf <- x$tf == "NFE2L1"; x$kn <- sprintf("%d of %d", x$shared_count, x$measured_targets)       # "of n": the targets measured
  x$qb <- sub("\\^.*$", "", x$matched_q_drawn); x$qe <- sub("^.*\\^", "", x$matched_q_drawn); x }   # "6.8 × 10", "−8"
FA <- prep(FA); FB <- prep(FB); FF <- rbind(FA, FB)
stopifnot(all(grepl(", [0-9,]+ targets\\)$", FF$label)), all(FF$rest != substring(FF$label, nchar(FF$tf) + 1L)))
labw <- vapply(seq_len(nrow(FF)), function(i) tw(FF$tf[i], 7, if (FF$nf[i]) "bold.italic" else "italic") + tw(FF$rest[i], 7, if (FF$nf[i]) "bold" else "plain"), 0)
LR <- .13 + max(max(labw), tw("Regulon (network)", 7, "bold"))
FP0 <- LR + .07; FP1 <- FP0 + 1.00
XL <- c(min(1, min(FF$fold_lo)) * .85, max(FF$fold_hi) * 1.12); fx <- function(f) FP0 + (log10(f) - log10(XL[1])) / diff(log10(XL)) * (FP1 - FP0)
BRL <- c(1, 3, 10, 30, 100); BRL <- BRL[BRL <= XL[2]]; BRG <- c(2, 5, 20, 50); BRG <- BRG[BRG <= XL[2]]
WKN <- max(tw(c(FF$kn, "Targets", "in set"), 7, "bold")); WQB <- max(tw(FF$qb)); WQE <- max(tw(FF$qe)); WQ <- max(WQB + WQE, tw("Matched q", 7, "bold"))
WRK <- max(tw(c("Rank", "of 1,113"), 7, "bold"))
CKN <- FP1 + .09 + WKN / 2; CQ <- CKN + WKN / 2 + .09 + WQ / 2; CRK <- CQ + WQ / 2 + .09 + WRK / 2; F_END <- CRK + WRK / 2
XQ0 <- CQ - (WQB + WQE) / 2 + WQB                                               # the q column aligned on the "× 10"
RF <- .106; SEPG <- .03                                                          # row pitch; gap at the dotted line
Y0 <- A_END + .14
HB <- Y0 + .14; HROW <- HB + .17 + LINE                                         # B heading; column header (two lines, second at HROW)
rows_of <- function(x, y0) { x$yc <- y0 + (x$row - .5) * RF + (x$part == "extra") * SEPG; x }
FA <- rows_of(FA, HROW + .045); B_END <- HROW + .045 + nrow(FA) * RF + any(FA$part == "extra") * SEPG
HC <- B_END + .20; FB <- rows_of(FB, HC + .045); C_ROWS_END <- HC + .045 + nrow(FB) * RF + any(FB$part == "extra") * SEPG
KW <- .8 / 2.54; KEY_Y1 <- C_ROWS_END + .04 + CAP7 + .122 + .15
F_BOT <- KEY_Y1 + 2 * LINE + .03
# D: networks down, factors across
DX0 <- F_END + .22
NETS <- unique(DT$network); NTYPE <- DT$network_group[match(NETS, DT$network)]; CTF <- unique(DT$factor); TYPES <- c("SC-islet", "adult islet", "intestine")
NLR <- DX0 + .08 + max(tw(NETS))
TWD <- (7.40 - NLR - .06 - 3 * .03) / 4; XCD <- setNames(NLR + .06 + (seq_along(CTF) - .5) * TWD + (seq_along(CTF) - 1) * .03, CTF)
stopifnot(max(tw(c(DT$tile_top, DT$tile_bottom)[nzchar(c(DT$tile_top, DT$tile_bottom))], 7, "bold")) <= TWD - .03)
RD <- .225; GH <- .15
yy <- Y0 + .45; YC <- c(); GHY <- c()
for (t in TYPES) { GHY[t] <- yy + .11; yy <- yy + GH; for (n in NETS[NTYPE == t]) { YC[n] <- yy + RD / 2; yy <- yy + RD } }
D_ROWS_END <- yy; DKY <- D_ROWS_END + .15; D_BOT <- DKY + 3 * LINE + .03
SW <- c(.45, .30) / 2.54
stopifnot(DX0 + SW[1] + .06 + tw("no regulon for this factor in this network") <= 7.45, max(F_BOT, D_BOT) <= 9.95)

# ==== panels ======================================================================================================================
draw_A <- function() {
  hdr2 <- function(lines, x0, x1) { for (k in seq_along(lines)) tx(lines[k], (x0 + x1) / 2, .145 + (k - 1) * LINE, 7, "bold", hjust = .5); sgm(x0, .29, x1, .29) }
  hdr2(c(sprintf("%d cholesterol-synthesis genes of Fig. 4D,", NG), "by mean graft change (smallest fall at left)"), X$hm0, X$hm1)
  hdr2(c("Median change", "(log2)"), X$bar0, X$val_r)
  w1 <- tw("NFE2L1", 7, "bold.italic"); w2 <- tw(" regulon", 7, "bold"); xm <- (X$c0 + X$c1) / 2
  grid.text("NFE2L1", U(xm - (w1 + w2) / 2), V(.145), hjust = 0, vjust = 0, gp = gp_t(7, "bold.italic"))
  tx(" regulon", xm - (w1 + w2) / 2 + w1, .145, 7, "bold")
  tx("target score change", xm, .145 + LINE, 7, "bold", hjust = .5); sgm(X$c0, .29, X$c1, .29)
  tx("Samples", X$samp_r, .668, hjust = 1); tx_arrow("start → end", X$samp_r, .668 + LINE, hjust = 1)
  for (k in seq_along(GO)) grid.text(GO[k], U(X$hm0 + (k - .5) * X$gw), V(.80), rot = 90, hjust = 0, vjust = .5, gp = gp_t(7, "italic"))
  for (v in c(-2, -1)) sgm(bx(v), BLK$graft$top, bx(v), A_BOT, col = "gray85")
  for (v in BRKC[BRKC != 0]) sgm(cx(v), BLK$graft$top, cx(v), A_BOT, col = "gray85")
  for (b in names(BLK)) tx(BLK[[b]]$head, X$lab, BLK[[b]]$base, 7, "bold")
  for (i in seq_len(nrow(LAB))) { r <- LAB[i, ]
    if (r$band) rct(X$band0 + A_INDENT, X$band1, r$top, r$top + CH, "#E8E8E8")
    lx <- X$lab + A_INDENT
    if (r$repeat_study) tx_arrow(r$change, lx + A_CHANGE_INDENT, r$center + .035) else {
      pieces(list(c(r$bold, "bold"), c(r$rest, "plain")), lx, r$center - .0236); tx_arrow(r$change, lx + A_CHANGE_INDENT, r$center - .0236 + LINE) }
    tx_arrow(r$samples, X$samp_r, r$center + .035, hjust = 1)
    for (p in POPS) rct(X$chip0, X$chip1, row_top(r$id, p), row_top(r$id, p) + RH, POP_COL[[p]], col = "white", lwd = .3 * PTMM) }
  # heatmap tiles drawn SC-SIC lines first and genes alphabetically within a line (neighboring tiles' white borders
  # overlap, so the drawing order shows in the pixels); a dot where every earlier x later sample pair changes the same way
  o <- order(match(TT$population, POPS), TT$gene, method = "radix")
  x0 <- X$hm0 + (TT$gene_column - 1) * X$gw; top <- row_top(TT$comparison, TT$population)
  rct(x0[o], x0[o] + X$gw, top[o], top[o] + RH, fill_of(TT$change[o]), col = "white", lwd = .3 * PTMM)
  d <- which(TT$dot); pt(x0[d] + X$gw / 2, top[d] + RH / 2, 16, .75, "gray15")
  # median bars and values; NFE2L1 bars; the two axes
  mt <- row_top(MED$comparison, MED$population)
  rct(pmin(bx(0), bx(MED$change)), pmax(bx(0), bx(MED$change)), mt + .022, mt + RH - .022, POP_COL[MED$population])
  for (j in seq_len(nrow(MED))) tx(MED$lab[j], X$val_r, mt[j] + RH / 2 + CAP7 / 2, hjust = 1)
  ct <- row_top(CV$comparison, CV$population)
  rct(pmin(cx(0), cx(CV$difference)), pmax(cx(0), cx(CV$difference)), ct + .022, ct + RH - .022, POP_COL[CV$population])
  for (ax in list(list(f = bx, lim = BAR_LIM, brk = c(-2, -1, 0)), list(f = cx, lim = LIMC, brk = BRKC))) {
    sgm(ax$f(0), BLK$graft$top, ax$f(0), Y_AX, lwd = .6); sgm(ax$f(ax$lim[1]), Y_AX, ax$f(ax$lim[2]), Y_AX)
    for (v in ax$brk) { sgm(ax$f(v), Y_AX, ax$f(v), Y_AX + .03); tx(minus(v), ax$f(v), Y_AX + .03 + .025 + CAP7, hjust = .5) } }
  # key row: cell types, color bar (230 steps), dot rule
  tx("Cell type", X$lab, KY); x <- X$lab
  for (p in POPS) { rct(x, x + .08, KY + .04, KY + .14, POP_COL[[p]]); tx(POP_LABEL[[p]], x + .12, KY + .14 - .02); x <- x + .12 + tw(POP_LABEL[[p]]) + .16 }
  tx("Change from the earlier in vitro samples, log2(1 + CPM)", cb0, KY)
  e <- seq(cb0, cb1, length.out = 231); v <- -3 + (seq_len(230) - .5) / 230 * 6
  rct(e[-231], e[-1], KY + .04, KY + .14, fill_of(v), col = fill_of(v), lwd = .2)
  for (k in -3:3) { xk <- cb0 + (k + 3) / 6 * (cb1 - cb0); sgm(xk, KY + .14, xk, KY + .17)
    tx(c("≤ −3", "−2", "−1", "0", "1", "2", "≥ 3")[k + 4], xk, KY + .265, hjust = .5) }
  dx <- cb1 + .30; pt(dx, KY - .035, 16, .75, "gray15"); for (k in seq_along(DOTL)) tx_arrow(DOTL[k], dx + .08, KY + (k - 1) * LINE)
}
# one section of the forest table: heading, NFE2L1's pale red band, gridlines, the dotted line above the extra rows (the
# other NFE2L1, SOX4 and SOX11 hits), and per row the interval, mark, label, targets in set, matched q and rank
sect <- function(x, hb, head) {
  tx(head, .32, hb, 7, "bold")
  y0 <- min(x$yc) - RF / 2; yb <- max(x$yc) + RF / 2
  for (i in which(x$nf)) rct(.10, F_END + .04, x$yc[i] - RF / 2, x$yc[i] + RF / 2, adjustcolor("#CC3311", .12))
  for (b in c(BRL, BRG)) sgm(fx(b), y0, fx(b), yb, col = "gray90", lwd = .3 * PTMM)
  sgm(fx(1), y0, fx(1), yb, col = "gray60", lwd = .3 * PTMM); rct(FP0, FP1, y0, yb, NA, col = "gray80", lwd = .3 * PTMM)
  if (any(x$part == "extra")) { ys <- y0 + sum(x$part == "top") * RF + SEPG / 2; sgm(.10, ys, F_END + .04, ys, col = "gray55", lwd = .4 * PTMM, lty = "dotted") }
  for (i in seq_len(nrow(x))) { r <- x[i, ]; b <- r$yc + CAP7 / 2
    sgm(fx(r$fold_lo), r$yc, fx(r$fold_hi), r$yc, col = GRP_COL[[r$grp]], lwd = .55 * PTMM); pt(fx(r$matched_fold), r$yc, GRP_PCH[[r$grp]], GRP_SZ[[r$grp]], GRP_COL[[r$grp]])
    pieces(list(c(r$tf, if (r$nf) "bold.italic" else "italic"), c(r$rest, if (r$nf) "bold" else "plain")), LR, b, hjust = 1)
    tx(r$kn, CKN, b, hjust = .5); tx(r$qb, XQ0, b, hjust = 1); tx(r$qe, XQ0, b - .019); tx(as.character(r$rank), CRK, b, hjust = .5) }
}
draw_B <- function() {
  sect(FA, HB, "Shared graft decline (369 genes)")
  tx("Regulon (network)", LR, HROW, 7, "bold", hjust = 1)
  for (h in list(list(CKN, c("Targets", "in set")), list(CQ, c("", "Matched q")), list(CRK, c("Rank", "of 1,113"))))
    for (k in 1:2) if (nzchar(h[[2]][k])) tx(h[[2]][k], h[[1]], HROW - (2 - k) * LINE, 7, "bold", hjust = .5)
}
draw_C <- function() {
  sect(FB, HC, "Cholesterol synthesis (Reactome; the 26 genes of Fig. 4D)")
  yb <- C_ROWS_END
  for (b in BRL) tx(format(b, trim = TRUE), fx(b), yb + .04 + CAP7, hjust = .5)
  tx("Matched fold enrichment (log scale), 95% interval", (FP0 + FP1) / 2, yb + .04 + CAP7 + .122, hjust = .5)
  # key: the three regulon marks, then the dotted line
  ky <- KEY_Y1; kx <- .10
  for (k in GRP_LEV[c(1, 3)]) { sgm(kx, ky - CAP7 / 2, kx + KW, ky - CAP7 / 2, col = GRP_COL[[k]], lwd = .55 * PTMM); pt(kx + KW / 2, ky - CAP7 / 2, GRP_PCH[[k]], GRP_SZ[[k]], GRP_COL[[k]])
    tx(k, kx + KW + .07, ky, face = if (k == "NFE2L1") "italic" else "plain")
    kx <- kx + KW + .07 + tw(k, 7, "italic") + .2 }
  k <- GRP_LEV[2]; ky2 <- ky + LINE; sgm(.10, ky2 - CAP7 / 2, .10 + KW, ky2 - CAP7 / 2, col = GRP_COL[[k]], lwd = .55 * PTMM); pt(.10 + KW / 2, ky2 - CAP7 / 2, GRP_PCH[[k]], GRP_SZ[[k]], GRP_COL[[k]])
  pieces(list(c("SOX4", "italic"), c(", ", "plain"), c("SOX11", "italic"), c(" (also fall in culture alone, where scored)", "plain")), .10 + KW + .07, ky2)
  ky3 <- ky2 + LINE; sgm(.10, ky3 - CAP7 / 2, .10 + KW, ky3 - CAP7 / 2, col = "gray55", lwd = .4 * PTMM, lty = "dotted")
  pieces(list(c("below the dotted line: the other ", "plain"), c("SOX4", "italic"), c(" and ", "plain"), c("SOX11", "italic"), c(" hit regulons", "plain")), .10 + KW + .07, ky3)
}
# one tile per factor and network (dark: hit; light: tested, not a hit; outlined: no regulon for the factor in that
# network), its matched fold (top) and targets in the set / measured (bottom); networks grouped by tissue; key below
draw_D <- function() {
  FL <- DT[!duplicated(DT$factor), ]
  for (f in CTF) { g <- grp_of(f); w1 <- tw(f, 7, "bold.italic"); x1 <- XCD[[f]] - (w1 + .1) / 2
    pt(x1 + .03, Y0 + .30 - CAP7 / 2, GRP_PCH[[g]], GRP_SZ[[g]], GRP_COL[[g]])
    tx(f, x1 + .1, Y0 + .30, 7, "bold.italic")
    tx(sub("^hit in ", "", FL$factor_label[FL$factor == f]), XCD[[f]], Y0 + .30 + LINE, hjust = .5) }
  tx("Hit in", NLR, Y0 + .30 + LINE, hjust = 1)
  for (t in TYPES) tx(t, DX0, GHY[[t]], 7, "bold")
  for (n in NETS) tx(n, DX0 + .08, YC[[n]] + CAP7 / 2)
  for (i in seq_len(nrow(DT))) { r <- DT[i, ]; xc <- XCD[[r$factor]]; yc <- YC[[r$network]]
    if (r$cell == "no list") { rct(xc - TWD / 2, xc + TWD / 2, yc - (RD - .02) / 2, yc + (RD - .02) / 2, "white", col = "gray75", lwd = .25 * PTMM); next }
    rct(xc - TWD / 2, xc + TWD / 2, yc - (RD - .02) / 2, yc + (RD - .02) / 2, if (r$cell == "hit") "#A6A6A6" else "#EDEDED")
    f <- if (r$cell == "hit") "bold" else "plain"
    tx(r$tile_top, xc, yc - .012, 7, f, hjust = .5); tx(r$tile_bottom, xc, yc + .093, 7, f, hjust = .5) }
  tx("Tile: matched fold (top); targets in set /", DX0, DKY); tx("targets measured (bottom)", DX0, DKY + LINE)
  kx <- DX0; ky <- DKY + 2 * LINE
  for (j in 1:2) { rct(kx, kx + SW[1], ky - CAP7 / 2 - SW[2] / 2, ky - CAP7 / 2 + SW[2] / 2, c("#A6A6A6", "#EDEDED")[j]); tx(c("hit", "tested, not a hit")[j], kx + SW[1] + .06, ky)
    kx <- kx + SW[1] + .06 + tw(c("hit", "tested, not a hit")[j]) + .2 }
  ky <- ky + LINE; rct(DX0, DX0 + SW[1], ky - CAP7 / 2 - SW[2] / 2, ky - CAP7 / 2 + SW[2] / 2, "white", col = "gray75", lwd = .25 * PTMM)
  tx("no regulon for this factor in this network", DX0 + SW[1] + .06, ky)
}

# ==== page places (px at 600 dpi): A across the page; B and C share the left column beside D ====================================
px <- function(v) as.integer(round(v * DPI))
XD <- px(DX0 - .07)                                                              # D's left edge; B and C end there
yB <- px(Y0); yC <- px(B_END + .07)
S2 <- list(A = list(x = 0L, w = 4500L, y = 0L, h = px(A_END + .02), draw = draw_A),
           B = list(x = 0L, w = XD, y = yB, h = yC - yB, draw = draw_B),
           C = list(x = 0L, w = XD, y = yC, h = px(F_BOT + .02) - yC, draw = draw_C),
           D = list(x = XD, w = 4500L - XD, y = yB, h = px(D_BOT + .02) - yB, draw = draw_D))
stopifnot(S2$A$h <= yB, yC + S2$C$h <= 6000L, yB + S2$D$h <= 6000L)
# one panel on its own strip device, placed on the page at its whole-pixel place
s2_panel <- function(name) { p <- S2[[name]]
  in_strip(p$y, p$h, function() { OX <<- p$x / DPI; OY <<- p$y / DPI; PH <<- p$h / DPI; p$draw() }, x_px = p$x, w_px = p$w) }
