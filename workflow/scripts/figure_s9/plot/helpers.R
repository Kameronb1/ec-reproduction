# Drawing helpers of Supplementary Figure S9: text, colors, theme, keys, the four panels' plots and their places on the
# page. Sourced by workflow/scripts/common/draw_panel.R, which opens the page device (7.5 x 10 in, 600 dpi, cairo) and
# defines read_table() and result_path(). The tables are Figure 5's results (regulon screen and CDX2 knockout). Every
# position is in inches from the page's top-left corner. Panels lower on the page are placed under the measured height of
# the panels above them, so each plot script lays out the whole page and draws only its own panel (render_s9()).
#   A  regulon target scores, EC minus other enteroendocrine cells, in the three intestinal studies and, set apart, in HIO
#      graft EC (McCauley 2023), with the count of reproducible regulons keeping the intestinal direction; rows as Figure 5B
#   B  G2 and G3 among the pooled targets of the 19 candidate regulators with targets higher in SC-SIC
#   C  cell-type shares in the CDX2 wild-type and knockout SC-islets, with the published shares
#   D  the 78 reproducible scored regulons after CDX2 knockout relative to matched random gene sets; the two reproducible
#      CDX2(+) regulons named by network

if (.Platform$OS.type == "windows") stopifnot(nzchar(Sys.setlocale("LC_CTYPE", "English_United States.utf8")))

W <- 7.5; H <- 10                                   # page size (in)
pt <- 7                                             # page text size (pt); every string is 7 pt or larger, letters 11 pt
TXT <- 7 / ggplot2::.pt; FONT <- "Arial"; REF_COL <- "#777777"; EDGE <- .15   # hairline black edge of every point
GCOL <- c(G1 = "#43A5B6", G2 = "#A89A2E", G3 = "#B4553F")          # gene groups G1-G3
GL <- c("G1: EC genes low in SC-SICs", "G2: shared EC genes", "G3: separate SC-SIC program")
POP_COL <- c("SC-SIC" = "#AA3377", "SC-beta" = "#4477AA", "SC-alpha" = "#117733")
BAND <- "#D9D9D9"                                   # row bands
SC_SIC <- POP_COL[["SC-SIC"]]; SC_BETA <- POP_COL[["SC-beta"]]
G2COL <- GCOL[["G2"]]; G3COL <- GCOL[["G3"]]
DOTCOL <- c(pass_ec = SC_SIC, pass_beta = SC_BETA, fail = "grey45")    # Figure 5B colors
minus <- function(v) sub("^-", "−", format(v, trim = TRUE, drop0trailing = TRUE))
fig5 <- function(name) read_table(name, figure = "figure5")

th <- function(grid_x = TRUE, grid_y = FALSE) theme_minimal(base_family = "Arial", base_size = pt) + theme(
  text = element_text(color = "black", size = pt), axis.text = element_text(color = "black", size = 7),
  axis.title = element_text(color = "black", size = 7), plot.title = element_text(color = "black", size = 7.5, face = "bold", margin = margin(0, 0, 2, 0)),
  strip.text = element_text(color = "black", size = 7, face = "bold", margin = margin(1, 0, 2, 0)),
  panel.grid.minor = element_blank(),
  panel.grid.major.x = if (grid_x) element_line(color = "#E4E4E4", linewidth = .2) else element_blank(),
  panel.grid.major.y = if (grid_y) element_line(color = "#E4E4E4", linewidth = .2) else element_blank(),
  axis.line.x = element_line(color = "black", linewidth = .35), axis.ticks.x = element_line(color = "black", linewidth = .3),
  axis.ticks.length = unit(1.3, "pt"), plot.margin = margin(2, 4, 2, 2), legend.position = "none", plot.title.position = "panel")
NOXL <- theme(axis.line.x = element_blank(), axis.ticks.x = element_blank())
text_at <- function(s, x, y, size = pt, bold = FALSE, col = "black", just = c("left", "top"), vp_h = H)
  grid.text(s, x = unit(x, "in"), y = unit(vp_h - y, "in"), just = just, gp = gpar(fontfamily = "Arial", fontsize = size, fontface = if (bold) "bold" else "plain", col = col))
str_in <- function(s, size = pt) { pushViewport(viewport(gp = gpar(fontfamily = "Arial", fontsize = size))); v <- convertWidth(stringWidth(s), "in", valueOnly = TRUE); popViewport(); v }
heading <- function(letter, x, y) text_at(letter, x, y, 11, TRUE)   # panel letter, 11 pt bold, top-left at (x, y)

# ---- rows: the 26 candidate regulators of Figure 5B, in its order (direction, then symbol) ----------------------------------
FA <- fig5("regulons/factors.tsv"); FA <- FA[FA$selected %in% c("TRUE", "True", "true"), ]
HALF <- c("Targets higher\nin SC-SIC", "Targets lower\nin SC-SIC")
FA$half <- factor(ifelse(FA$direction == "SC-SIC higher", HALF[1], HALF[2]), HALF)
FA$tf <- sub("\\(.*$", "", FA$factor); FA <- FA[order(FA$half, FA$tf, FA$factor), ]
FA$sym <- ifelse(grepl("\\(\\+\\)$", FA$factor), "+", "−")
FA$expr <- sprintf('italic("%s")~"(%s)"', FA$tf, FA$sym)
FA$y <- as.numeric(factor(FA$factor, rev(FA$factor))); FA$d <- ifelse(FA$direction == "SC-SIC higher", 1, -1)
stopifnot(nrow(FA) == 26, sum(FA$d == 1) == 19)
PV <- fig5("regulons/program_values.tsv"); PV$pass <- PV$program_pass %in% c("TRUE", "True", "true")
PV$col <- ifelse(!PV$pass, "fail", ifelse(PV$direction == "SC-SIC higher", "pass_ec", "pass_beta"))
PV <- PV[order(match(PV$factor, FA$factor), PV$network), ]
PV$off <- ave(seq_len(nrow(PV)), PV$factor, FUN = function(i) if (length(i) == 1) 0 else seq(-.30, .30, length.out = length(i)))
PV$y <- FA$y[match(PV$factor, FA$factor)] + PV$off; PV$half <- FA$half[match(PV$factor, FA$factor)]
stopifnot(nrow(PV) == 122, sum(PV$pass) == 79)
BANDS <- transform(FA[seq_len(nrow(FA)) %% 2 == 0, c("factor", "half", "y")], ymin = y - .5, ymax = y + .5)
HB <- do.call(rbind, lapply(HALF, function(h) { yy <- FA$y[FA$half == h]; data.frame(half = factor(h, HALF), y = c(min(yy) - .5, max(yy) + .5)) }))

# ---- A: values --------------------------------------------------------------------------------------------------------------
SV <- fig5("regulons/supplement_values.tsv")
GC <- SV[SV$part == "gut_study_mean", c("factor", "program", "item", "value", "note")]   # intestinal study means (Elmentaite: mean of its parts)
ITEM <- c(ELMENTAITE = "Elmentaite 2020 + 2021", EGOZI_2021 = "Egozi 2021", HICKEY_2023 = "Hickey 2023")   # Figure 1B order
GC <- merge(GC, PV[, c("program", "y", "half", "col", "factor")], by = c("program", "factor"))
stopifnot(all(GC$note == GC$col))                           # the study-mean colors are Figure 5B's colors
for (k in names(ITEM)) { v <- PV[[paste0("gut_mean__", k)]][match(GC$program[GC$item == k], PV$program)]
  stopifnot(max(abs(v - GC$value[GC$item == k])) < 1e-12) } # and its values are Figure 5B's program_values
GC$col_id <- unname(ITEM[GC$item])
MC <- fig5("regulons/mccauley2023_values.tsv")             # McCauley 2023 HIO graft EC, published labels (Endocrine_1 vs Endocrine_2)
OC <- data.frame(study = "McCauley 2023 HIO graft EC", factor = MC$factor, program = MC$program, effect = MC$gse214852_effect,
                 col = PV$col[match(MC$program, PV$program)])
OST <- c("McCauley 2023 HIO graft EC" = "McCauley 2023")
OC <- merge(OC[, c("study", "factor", "program", "effect", "col")], PV[, c("program", "y", "half", "col")], by = "program", suffixes = c("", ".pv"))
stopifnot(nrow(OC) == 122, all(OC$col == OC$col.pv))
OC$col_id <- unname(OST[OC$study]); OC$value <- OC$effect
pp <- OC[OC$col != "fail", ]; pp$d <- FA$d[match(pp$factor, FA$factor)]
CN <- aggregate(cbind(n = 1, k = as.integer(sign(pp$value) == pp$d)) ~ col_id + factor, pp, sum)
CN$lab <- sprintf("%d of %d", CN$k, CN$n); CN$y <- FA$y[match(CN$factor, FA$factor)]; CN$half <- FA$half[match(CN$factor, FA$factor)]
CT <- aggregate(cbind(k, n) ~ col_id, CN, sum)
stopifnot(nrow(CT) == 1, CT$k[CT$col_id == "McCauley 2023"] == 78, CT$n[CT$col_id == "McCauley 2023"] == 79)
A_DOTS <- rbind(GC[, c("col_id", "factor", "program", "value", "col")], OC[, c("col_id", "factor", "program", "value", "col")])

# ---- A: drawing (one plot per column; label column, four dot columns, one count column; widths in inches) ----------------
PITCH_A <- .108; HGAP_A <- .10                              # row pitch; gap between the two row groups (in)
thA <- theme_void(base_family = FONT, base_size = 7) + theme(
  strip.text = element_blank(), panel.spacing.y = unit(HGAP_A, "in"), plot.margin = margin(0, 0, 0, 0),
  plot.title = element_text(size = 7, face = "bold", hjust = 0, margin = margin(0, 0, 1, 0)),
  plot.subtitle = element_text(size = 7, hjust = 0, vjust = 1, lineheight = .95, margin = margin(0, 0, 3, 0)), plot.title.position = "panel")
ysA <- function() list(facet_grid(half ~ ., scales = "free_y", space = "free_y"), geom_blank(data = HB, aes(x = 0, y = y), inherit.aes = FALSE),
                       scale_y_continuous(expand = c(0, 0)))
dotcol <- function(key) {
  d <- A_DOTS[A_DOTS$col_id == key, ]; d$y <- PV$y[match(d$program, PV$program)]; d$half <- PV$half[match(d$program, PV$program)]
  L <- max(abs(d$value)) * 1.08; br <- scales::breaks_pretty(3)(c(-L, L)); br <- br[abs(br) <= L * .93]
  ggplot() + geom_rect(data = BANDS, aes(xmin = -Inf, xmax = Inf, ymin = ymin, ymax = ymax), fill = BAND, colour = NA) +
    geom_vline(xintercept = br, colour = "#E4E4E4", linewidth = .2) +
    geom_vline(xintercept = 0, colour = REF_COL, linewidth = .3) + ysA() +
    geom_point(data = d, aes(value, y, fill = col), shape = 21, colour = "black", stroke = EDGE, size = 1.05) +
    scale_fill_manual(values = DOTCOL, guide = "none") +
    scale_x_continuous(limits = c(-L, L), breaks = br, labels = minus, expand = c(0, 0)) + coord_cartesian(clip = "off") +
    labs(x = NULL, y = NULL) + thA +
    theme(axis.text.x = element_text(size = 7, colour = "black", margin = margin(2, 0, 0, 0)), axis.line.x = element_line(colour = "black", linewidth = .35),
          axis.ticks.x = element_line(colour = "black", linewidth = .3), axis.ticks.length.x = unit(1.3, "pt"))
}
cntcol <- function(id) {
  ggplot() + geom_rect(data = BANDS, aes(xmin = -Inf, xmax = Inf, ymin = ymin, ymax = ymax), fill = BAND, colour = NA) + ysA() +
    geom_text(data = CN[CN$col_id == id, ], aes(x = .5, y = y, label = lab), size = TXT, family = FONT) +
    scale_x_continuous(limits = c(0, 1), expand = c(0, 0)) + coord_cartesian(clip = "off") +
    labs(x = NULL, y = NULL) + thA +
    theme(plot.subtitle = element_text(hjust = .5, size = 7, lineheight = .95, margin = margin(0, 0, 3, 0)),
          axis.text.x = element_text(size = 7, colour = NA, margin = margin(2, 0, 0, 0)))   # an invisible tick row, so the rows align
}
labcol <- ggplot() + ysA() +
  geom_text(data = FA, aes(x = 1, y = y, label = expr), parse = TRUE, hjust = 1, size = TXT, family = FONT) +
  geom_text(data = data.frame(half = factor(HALF, HALF), y = sapply(HALF, function(h) mean(FA$y[FA$half == h]))),
            aes(x = 0, y = y, label = half), hjust = 0, fontface = "bold", size = TXT, family = FONT, lineheight = .95) +
  scale_x_continuous(limits = c(0, 1), expand = c(0, 0)) + coord_cartesian(clip = "off") + labs(x = NULL, y = NULL) + thA +
  theme(axis.text.x = element_text(size = 7, colour = NA, margin = margin(2, 0, 0, 0)))
W_LAB <- 1.40; W_DOT <- 1.00; W_CNT <- .60; G1 <- .14; G_HELD <- .22; G_CNT <- .03
A_W <- c(W_LAB, G1, W_DOT, G1, W_DOT, G1, W_DOT, G_HELD, W_DOT, G_CNT, W_CNT)
sp <- function() plot_spacer() + theme(plot.margin = margin(0, 0, 0, 0))
A_ROWS_H <- (nrow(FA) + 0) * PITCH_A
pA <- wrap_plots(list(labcol, sp(), dotcol(ITEM[["ELMENTAITE"]]), sp(), dotcol(ITEM[["EGOZI_2021"]]), sp(), dotcol(ITEM[["HICKEY_2023"]]),
                      sp(), dotcol("McCauley 2023"), sp(), cntcol("McCauley 2023")),
                 nrow = 1, widths = unit(A_W, "in"), heights = unit(A_ROWS_H + HGAP_A, "in")) &
  theme(plot.background = element_blank())
A_X0 <- .10 + .02                                           # left edge of A's label column on the page (in)
X_HELD <- A_X0 + sum(A_W[1:8])                              # left edge of the held-out block
X_DOT <- A_X0 + cumsum(c(0, A_W))[c(3, 5, 7)]               # left edges of the intestinal dot columns
A_XT_C <- (X_DOT[1] + X_HELD + A_W[9]) / 2                  # x title centered under all dot columns

# ---- B: G2 and G3 among the targets of the 19 candidate regulators with targets higher in SC-SIC --------------------------
U <- SV[SV$part == "targets_G2_of_G2G3", ]
up <- FA[FA$half == HALF[1], ]; Bd <- U[match(up$factor, U$factor), ]
stopifnot(!anyNA(Bd$factor), nrow(Bd) == 19)
B <- data.frame(factor = Bd$factor, expr = up$expr, g2 = as.integer(Bd$count), n = as.integer(Bd$total), e2 = Bd$expected)
B$g3 <- B$n - B$g2; B$y <- rev(seq_len(nrow(B)))
R2 <- mean(B$e2 / B$n)                                      # G2 share of all G2 + G3 genes (34 of 118)
stopifnot(abs(R2 - 34 / 118) < .005)
TS <- fig5("regulons/tests.tsv"); TS <- TS[TS$level == "factor", ]
B$lean <- ifelse(B$g2 > B$e2, "G2", "G3")
B$q <- TS$q[match(paste(B$factor, paste0("lean_", B$lean)), paste(TS$unit, TS$test))]; stopifnot(!anyNA(B$q))
B$txt <- sprintf("%d, %d of %d", B$g2, B$g3, B$n); B$txt2 <- paste(B$lean, sprintf("%.2f", B$q))
PITCH_B <- .115; BARF <- .64
pB <- ggplot(B) +
  geom_rect(aes(xmin = 0, xmax = g2 / n, ymin = y - BARF / 2, ymax = y + BARF / 2), fill = G2COL) +
  geom_rect(aes(xmin = g2 / n, xmax = 1, ymin = y - BARF / 2, ymax = y + BARF / 2), fill = G3COL) +
  geom_vline(xintercept = R2, linetype = "22", linewidth = .45, colour = "black") +
  scale_y_continuous(breaks = B$y, labels = parse(text = B$expr), limits = c(.5, nrow(B) + .5), expand = c(0, 0)) +
  scale_x_continuous(limits = c(0, 1), breaks = c(0, .5, 1), labels = c("0%", "50%", "100%"), expand = c(0, 0)) +
  coord_cartesian(clip = "off") + labs(x = "Share of the G2 + G3 targets of each entry", y = NULL) + th(FALSE) +
  theme(axis.text.y = element_text(hjust = 1, size = 7, margin = margin(r = 3)), plot.margin = margin(0, 0, 0, 0),
        axis.title.x = element_text(size = 7, margin = margin(3, 0, 0, 0)))
B_PW <- 2.30; B_TW <- c(.08, .95, .06, .62)

# ---- C: knockout cell-type proportions ----------------------------------------------------------------------------------------
fa <- fig5("cdx2_knockout/figure5d_comparison.tsv"); fa <- fa[fa$arm == "all_stages", ]
LABK <- c(SC_EC = "SC-SIC", SC_beta = "SC-beta", SC_alpha = "SC-alpha", SC_delta = "SC-delta", ENP3 = "ENP3")
fa$lab <- factor(LABK[fa$label], LABK); fa$x <- ifelse(fa$genotype == "WT", 1, 2)
wt <- fa[fa$genotype == "WT", ]; ko <- fa[fa$genotype == "CDX2-KO", ]; ko <- ko[match(wt$label, ko$label), ]
ratio <- ko$ours[wt$label == "SC_EC"] / wt$ours[wt$label == "SC_EC"]   # knockout share over wild-type share
stopifnot(nrow(fa) == 10)
KCOL <- c(`SC-SIC` = SC_SIC, `SC-beta` = SC_BETA, `SC-alpha` = POP_COL[["SC-alpha"]], `SC-delta` = "#AE8C44", ENP3 = "grey60")
KL <- data.frame(lab = as.character(ko$lab), y = 100 * ko$ours); KL <- KL[order(KL$y), ]
KL$ly <- KL$y; for (i in seq_len(nrow(KL))[-1]) KL$ly[i] <- max(KL$y[i], KL$ly[i - 1] + 7.5)
KL$text <- ifelse(KL$lab == "SC-SIC", sprintf("SC-SIC ×%.2f", ratio), KL$lab); KLL <- KL[(KL$ly - KL$y) > 2.5, ]
pC <- ggplot(fa, aes(x, 100 * ours, colour = lab, group = lab)) +
  annotate("segment", x = .8, xend = 2.2, y = seq(0, 60, 20), yend = seq(0, 60, 20), colour = "#E4E4E4", linewidth = .2) +
  geom_line(aes(linewidth = lab == "SC-SIC")) +
  geom_segment(aes(x = x + ifelse(x == 1, -.12, .12) - .045, xend = x + ifelse(x == 1, -.12, .12) + .045, y = 100 * zhu_fig5d, yend = 100 * zhu_fig5d),
               linewidth = .35, lineend = "butt") +
  geom_point(aes(fill = lab), shape = 21, colour = "black", stroke = EDGE, size = 1.35) +
  geom_segment(data = KLL, aes(x = 2.20, xend = 2.29, y = y, yend = ly), inherit.aes = FALSE, linewidth = .2, colour = "black") +
  geom_text(data = KL, aes(x = 2.32, y = ly, label = text, fontface = ifelse(lab == "SC-SIC", "bold", "plain")), inherit.aes = FALSE,
            hjust = 0, size = TXT, family = FONT, colour = "black") +
  annotate("text", x = c(1, 2), y = -6, label = c("Wild", "CDX2"), fontface = c("plain", "italic"), vjust = 1, size = TXT, family = FONT) +
  annotate("text", x = c(1, 2), y = -6 - 4.6, label = c("type", "knockout"), vjust = 1, size = TXT, family = FONT) +
  scale_linewidth_manual(values = c(`TRUE` = 1.0, `FALSE` = .45), guide = "none") + scale_colour_manual(values = KCOL, guide = "none") +
  scale_fill_manual(values = KCOL, guide = "none") +
  scale_x_continuous(limits = c(.8, 3.25), expand = c(0, 0)) +
  scale_y_continuous(limits = c(0, 66), breaks = seq(0, 60, 20), oob = scales::oob_keep) + coord_cartesian(clip = "off") +
  labs(x = NULL, y = "% of nuclei\n(five types = 100%)") + th(FALSE, FALSE) + NOXL +
  theme(axis.text.x = element_blank(), plot.margin = margin(2, 4, 22, 2))

# ---- D: the reproducible scored regulons after knockout -----------------------------------------------------------------
KC <- fig5("cdx2_knockout/regulon_targets/context_regulons.tsv")
KC <- KC[KC$scored %in% c("True", TRUE) & KC$reproducible %in% c("True", TRUE), ]
stopifnot(nrow(KC) == 78)
KHALF <- c(h = "Targets higher\nin SC-SIC", l = "Targets lower\nin SC-SIC")
KC$g <- factor(ifelse(KC$factor_direction == "SC-SIC higher", "h", "l"), c("h", "l")); KC$half <- factor(KHALF[as.character(KC$g)], KHALF)
set.seed(7); KC$jy <- 1 + runif(nrow(KC), -.16, .16)       # vertical jitter of the dots within a row
KC$cdx <- grepl("__CDX2\\(", KC$regulon); stopifnot(sum(KC$cdx) == 2)
K_NM <- c(GSE167880 = "Balboa 2022", GSE185036 = "Du 2022") # the networks of the two reproducible CDX2(+) regulons
KC$lab <- ifelse(KC$cdx, K_NM[KC$network], ""); stopifnot(!anyNA(KC$lab)); KC$jy[KC$cdx] <- 1
KM <- aggregate(excess_over_random ~ g + half, KC, median)
KN <- aggregate(cbind(below = excess_over_random < 0, n = 1) ~ g + half, KC, sum); KN$txt <- sprintf("%d of %d", KN$below, KN$n)
XB <- range(KC$excess_over_random); XB <- XB + c(-.03, .03) * diff(XB)
pD <- ggplot(KC) +
  geom_vline(xintercept = 0, colour = REF_COL, linewidth = .3) +
  geom_point(data = KC[!KC$cdx, ], aes(excess_over_random, jy), shape = 21, fill = "grey72", colour = "black", stroke = EDGE, size = 1.1) +
  geom_segment(data = KM, aes(x = excess_over_random, xend = excess_over_random, y = 1 - .2, yend = 1 + .2), colour = "black", linewidth = .6) +
  geom_point(data = KC[KC$cdx, ], aes(excess_over_random, jy), shape = 21, fill = SC_SIC, colour = "black", stroke = EDGE, size = 2.1) +
  geom_text(data = KC[KC$cdx, ], aes(excess_over_random, jy + .25, label = lab), size = TXT, family = FONT, vjust = 0) +
  geom_text(data = KN, aes(x = XB[2] + .04 * diff(XB), y = 1, label = txt), hjust = 0, size = TXT, family = FONT) +
  geom_text(data = data.frame(half = factor(KHALF[1], KHALF)), aes(x = XB[2] + .04 * diff(XB), y = 1.95, label = "Below 0"), hjust = 0, vjust = 1,
            size = TXT, family = FONT, fontface = "bold") +
  facet_grid(half ~ ., switch = "y") +
  scale_y_continuous(breaks = NULL, limits = c(.68, 1.72), expand = c(0, 0), oob = scales::oob_keep) +
  scale_x_continuous(limits = XB, breaks = seq(-.2, .6, .2), labels = minus, expand = c(0, 0), oob = scales::oob_keep) +
  coord_cartesian(clip = "off") +
  labs(x = quote(italic("CDX2") ~ "knockout minus wild type, relative to matched random gene sets (pseudobulk, log2)"), y = NULL) +
  th(TRUE, FALSE) +
  theme(plot.margin = margin(14, 44, 2, 2), strip.placement = "outside",
        strip.text.y.left = element_text(angle = 0, hjust = 1, face = "bold", size = 7, lineheight = .95, margin = margin(0, 5, 0, 0)),
        panel.background = element_rect(fill = "#F3F3F3", colour = NA), panel.spacing.y = unit(.06, "in"),
        axis.title.x = element_text(size = 7, margin = margin(3, 0, 0, 0)))

# ---- keys (shapes in a row; labels may be plotmath expressions) ----------------------------------------------------------
kw <- function(it) { if (is.call(it$label)) it$label <- as.expression(it$label)
  w <- if (is.expression(it$label)) { pushViewport(viewport(gp = gpar(fontfamily = FONT, fontsize = 7)))
    v <- convertWidth(stringWidth(it$label), "in", valueOnly = TRUE); popViewport(); v } else str_in(it$label, 7)
  w + switch(it$type, ln = .30, dash = .17, txt = 0, .15) }
keyrow <- function(items, cx, y, gap = .18, head = NULL) {   # y: inches from the page top (row center)
  hw <- if (is.null(head)) 0 else str_in(head, 7) + .12
  ww <- vapply(items, kw, 0); x <- cx - (hw + sum(ww) + gap * (length(ww) - 1)) / 2; yy <- unit(H - y, "in")
  if (!is.null(head)) { grid.text(head, unit(x, "in"), yy, just = c("left", "centre"), gp = gpar(fontfamily = FONT, fontsize = 7)); x <- x + hw }
  for (i in seq_along(items)) { it <- items[[i]]
    if (it$type == "pt") grid.points(unit(x + .05, "in"), yy, pch = 21, size = unit(it$size %||% 1.6, "mm"),
                                     gp = gpar(col = "black", fill = it$fill, lwd = EDGE * .stroke / 2))
    if (it$type == "sq") grid.rect(unit(x, "in"), yy, unit(.095, "in"), unit(.095, "in"), just = "left", gp = gpar(fill = it$fill, col = NA))
    if (it$type == "ln") grid.lines(unit(c(x, x + .22), "in"), unit.c(yy, yy), gp = gpar(col = "black", lwd = it$lwd %||% 1, lty = it$lty %||% "solid", lineend = "butt"))
    if (it$type == "dash") grid.lines(unit(c(x + .005, x + .095), "in"), unit.c(yy, yy), gp = gpar(col = "black", lwd = it$lwd %||% 1, lineend = "butt"))
    if (it$type == "tick") grid.lines(unit(c(x + .05, x + .05), "in"), unit.c(yy - unit(.055, "in"), yy + unit(.055, "in")), gp = gpar(col = "black", lwd = .6 * .pt, lineend = "butt"))
    off <- switch(it$type, ln = .30, dash = .17, txt = 0, .15)
    grid.text(if (is.call(it$label)) as.expression(it$label) else it$label, unit(x + off, "in"), yy, just = c("left", "centre"), gp = gpar(fontfamily = FONT, fontsize = 7))
    x <- x + ww[i] + gap }
}
gt_height <- function(g) convertHeight(sum(g$heights), "in", valueOnly = TRUE)
fix_panels <- function(g, h) { pr <- unique(g$layout$t[grepl("^panel", g$layout$name)]); g$heights[pr] <- unit(h, "in"); g }   # each panel row h in high
# a gtable with its top-left corner at (x, y) and width w; returns its height; drawn only when draw is TRUE
place <- function(g, x, y, w, draw = TRUE) { h <- gt_height(g)
  if (draw) { pushViewport(viewport(x = unit(x, "in"), y = unit(H - y, "in"), width = unit(w, "in"), height = unit(h, "in"), just = c("left", "top")))
    grid.draw(g); popViewport() }
  h }

# ---- page --------------------------------------------------------------------------------------------------------------------
LET_X <- .10; GAP_PANEL <- .20
# column header: lines bottom-aligned at ybot (in from the page top); faces per line
hdr <- function(lines, faces, x, ybot, just = "left") { n <- length(lines)
  for (i in seq_len(n)) grid.text(lines[i], unit(x, "in"), unit(H - ybot + (n - i) * .115, "in"), just = c(just, "bottom"),
                                  gp = gpar(fontfamily = FONT, fontsize = 7, fontface = faces[i])) }
XB_C <- 5.20; WC <- W - XB_C - .06                          # C to the right of B
# Lay out the page and draw the panels named in `only` ("A" to "D"); returns the page positions.
render_s9 <- function(only) {
  POS <- list(); on <- function(k) k %in% only
  # A: headers, then the plot whose top is the panel top
  yA <- .10; topA <- yA + .58; hb <- topA - .05             # panel top; header baseline band
  xs <- A_X0 + cumsum(c(0, A_W))                             # left edge of every column
  if (on("A")) {
    heading("A", LET_X, yA)
    hdr(c("Elmentaite", "2020 + 2021", "mean of 3 parts"), c("bold", "bold", "plain"), xs[3], hb)
    hdr("Egozi 2021", "bold", xs[5], hb); hdr("Hickey 2023", "bold", xs[7], hb)
    hdr(c("McCauley 2023", "HIO graft EC"), c("bold", "plain"), xs[9], hb)
    hdr(c("Regulons with", "intestinal", "direction"), rep("plain", 3), xs[11] + A_W[11] / 2, hb, "centre")
    xr <- xs[12]
    grid.text("Not used for selection", unit((X_HELD + xr) / 2, "in"), unit(H - yA - .06, "in"), just = c("centre", "centre"), gp = gpar(fontfamily = FONT, fontsize = 7, fontface = "bold"))
    grid.lines(unit(c(X_HELD, xr), "in"), unit(rep(H - yA - .14, 2), "in"), gp = gpar(col = "black", lwd = .5))
  }
  hA <- place(patchworkGrob(pA), A_X0, topA, sum(A_W), on("A"))
  yAx <- topA + hA + .015
  k1 <- yAx + .24
  if (on("A")) {
    grid.text("EC minus other enteroendocrine, regulon target score", unit(A_XT_C, "in"), unit(H - yAx, "in"), just = c("centre", "top"), gp = gpar(fontfamily = FONT, fontsize = 7))
    keyrow(list(list(type = "pt", label = "higher in SC-SIC", fill = DOTCOL[["pass_ec"]]), list(type = "pt", label = "lower in SC-SIC", fill = DOTCOL[["pass_beta"]]),
                list(type = "pt", label = "not reproducible", fill = DOTCOL[["fail"]])), A_X0 + sum(A_W) / 2, k1, head = "Reproducible in SC-islets and intestine:")
    keyrow(list(list(type = "txt", label = "(+) activated targets"), list(type = "txt", label = "(−) repressed targets")), A_X0 + sum(A_W) / 2, k1 + .15, gap = .3)
  }
  POS$A <- k1 + .15 + .08
  # B (left) and C (right), tops aligned
  yB <- POS$A + GAP_PANEL
  if (on("B")) { heading("B", LET_X, yB)
    grid.text("Targets higher in SC-SIC", unit(A_X0, "in"), unit(H - yB - .25, "in"), just = c("left", "centre"), gp = gpar(fontfamily = FONT, fontsize = 7, fontface = "bold")) }
  if (on("C")) heading("C", XB_C - .02, yB)
  gB <- ggplotGrob(pB + theme(plot.background = element_blank())); pr <- unique(gB$layout$t[gB$layout$name == "panel"]); gB$heights[pr] <- unit(nrow(B) * PITCH_B, "in")
  pc <- unique(gB$layout$l[gB$layout$name == "panel"]); gB$widths[pc] <- unit(B_PW, "in")
  wB <- convertWidth(sum(gB$widths), "in", valueOnly = TRUE)
  xB0 <- A_X0 + .10
  hB <- place(gB, xB0, yB + .36, wB, on("B"))
  ptop <- yB + .36 + convertHeight(sum(gB$heights[seq_len(pr - 1)]), "in", valueOnly = TRUE)   # top of B's bar area
  xP <- xB0 + wB - convertWidth(sum(gB$widths[(pc + 1):ncol(gB)]), "in", valueOnly = TRUE)       # right end of the bars
  xT <- xP + c(B_TW[1], B_TW[1] + B_TW[2] + B_TW[3])
  kB <- yB + .36 + hB + .14
  if (on("B")) {
    for (k in 1:2) { grid.text(c("G2, G3 of total", "Lean, q")[k], unit(xT[k], "in"), unit(H - ptop + .07, "in"), just = c("left", "bottom"),
                               gp = gpar(fontfamily = FONT, fontsize = 7, fontface = "bold"))
      grid.text(B[[c("txt", "txt2")[k]]], unit(xT[k], "in"), unit(H - ptop - (nrow(B) - B$y + .5) * PITCH_B, "in"), just = c("left", "centre"),
                gp = gpar(fontfamily = FONT, fontsize = 7)) }
    keyrow(list(list(type = "sq", label = GL[2], fill = G2COL), list(type = "sq", label = GL[3], fill = G3COL)), (xB0 + xT[2] + .6) / 2, kB)
    keyrow(list(list(type = "ln", label = "G2 share of all G2 + G3 genes", lty = "22", lwd = .45 * .pt)), (xB0 + xT[2] + .6) / 2, kB + .15)
  }
  gC <- fix_panels(ggplotGrob(pC + theme(plot.background = element_blank())), 1.75); hC <- place(gC, XB_C, ptop - .05, WC, on("C"))
  kC <- ptop - .05 + hC + .06
  if (on("C")) keyrow(list(list(type = "pt", label = "Reanalysis", fill = "grey35", size = 1.4), list(type = "dash", label = "Published shares")), XB_C + WC / 2, kC, gap = .14)
  POS$B <- kB + .15 + .08; POS$C <- kC + .08
  # D: full width
  yD <- max(POS$B, POS$C) + GAP_PANEL
  if (on("D")) heading("D", LET_X, yD)
  gD <- fix_panels(ggplotGrob(pD + theme(plot.background = element_blank())), .42); hD <- place(gD, A_X0, yD + .02, 6.6, on("D"))
  kD <- yD + .02 + hD + .10
  if (on("D")) keyrow(list(list(type = "pt", label = quote(italic("CDX2") * "(+) regulon"), fill = SC_SIC, size = 2.1), list(type = "tick", label = "Median of the row"),
                           list(type = "pt", label = "Other reproducible regulon of Fig. 5B", fill = "grey72", size = 1.4)), A_X0 + 3.3, kD)
  POS$D <- kD + .08; POS$yB <- yB; POS$yD <- yD
  invisible(POS)
}
