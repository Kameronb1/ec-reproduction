# Figure 2C: four-marker classifier call rate within each primary study label (bars: called / cells with that
# label) and, right of the bars, the classifier's four markers (TPH1, SLC18A1, FEV, LMX1A) as dots: size = share of
# the label's cells with at least one read of the marker, color = the mean log expression of the marker among those
# cells. All five studies are stacked in one plot; the gene names are written once, under the last study; two key
# lines below.

XTITLE <- "Called EC (% of cells with that label)"
TOP_C <- 3.254                                # cap top of the first study title (inches from the page top)

labels_C <- read_table("C_primary_labels.tsv")[c("family", "author_label_original", "author_label_class", "n", "called", "rate", "label")]
# McCauley 2023 HIO grafts: Endocrine_1 (EC) and Endocrine_2 (other terminal EEC), HIO sample groups pooled
hio_C <- aggregate(cbind(n, rule_called) ~ author_label, read_table("hio_BC_classifier_calls_by_label.tsv"), sum)
hio_C <- hio_C[order(hio_C$author_label), ]
labels_C <- rbind(labels_C, data.frame(family = "hio", author_label_original = hio_C$author_label,
                                       author_label_class = c("EXPLICIT_EC", "TERMINAL_NON_EC_EEC"), n = hio_C$n,
                                       called = hio_C$rule_called, rate = 100 * hio_C$rule_called / hio_C$n,
                                       label = paste0(hio_C$rule_called, "/", hio_C$n)))
CLASS_COL <- c(EXPLICIT_EC = "#43A5B6", TERMINAL_NON_EC_EEC = "#8C8C8C", BROAD_EEC_PARENT_UNRESOLVED = "#D3D6D9",
               ENDOCRINE_PROGENITOR = "#444444")
EC_COL <- c(gse238274 = "#CCBB44", hio = "#6943AC")   # EC bars in the study's own color (tissue EC: teal)

# The markers per label: cells with >= 1 raw read, their share, and their mean log1p CPM (NA where no cell has a
# read), from workflow/scripts/figure2/c_marker_detection.R. Its cell counts and EC counts must equal the bars.
GENES_C <- c("TPH1", "SLC18A1", "FEV", "LMX1A")
mm_C <- read_table("C_marker_detection.tsv")
k <- match(paste(mm_C$family, mm_C$author_label_original), paste(labels_C$family, labels_C$author_label_original))
stopifnot(!anyNA(k), nrow(mm_C) == 4 * nrow(labels_C), mm_C$n == labels_C$n[k], mm_C$called == labels_C$called[k],
          !duplicated(paste(mm_C$family, mm_C$author_label_original, mm_C$gene)), mm_C$gene %in% GENES_C,
          abs(mm_C$share_detecting - mm_C$cells_detecting / mm_C$n) < 1e-12,
          is.na(mm_C$mean_log_detecting) == (mm_C$cells_detecting == 0))

# ---- vertical layout ----
# y is in inches below the plot panel's top (drawn negative), so each spacing below is a distance on the page.
fh <- function(s, size, face = "plain")       # ascent of a string in inches on the open device
  convertHeight(grobAscent(textGrob(s, gp = gpar(fontfamily = "Arial", fontsize = size, fontface = face))), "in", valueOnly = TRUE)
CAP8 <- fh("E", 8, "bold")                    # cap height of the 8 pt bold study titles
P_IN <- .101                                  # row pitch (>= 0.101 in: neighboring labels' parentheses and underscores clear)
T2R_IN <- .0745                               # a study title's baseline to its first row's center (title descenders
                                              # such as '|' and 'g' clear the first bar by about 0.015 in)
R2S_IN <- .060                                # a study's last row center to its bottom rule ...
R2S_DESC <- .068                              # ... or, when the last called/total has a comma, to clear the comma
GAP_IN <- .0365                               # a bottom rule's center to the next study title's cap top (0.030 in from
                                              # the 0.35-linewidth rule's lower edge)
TOP_IN <- CAP8 + .004                         # the panel top to the first title's baseline
yy <- TOP_IN; blocks_C <- list()
for (f in FAMS) {
  d0 <- labels_C[labels_C$family == f, ]
  ry <- yy + T2R_IN + (seq_len(nrow(d0)) - 1) * P_IN
  r2s <- if (grepl(",", d0$label[nrow(d0)], fixed = TRUE)) R2S_DESC else R2S_IN
  blocks_C[[f]] <- list(d = d0, title = yy, rows = ry, sep = max(ry) + r2s)
  yy <- max(ry) + r2s + GAP_IN + CAP8
}
TOTAL_C <- blocks_C[[length(blocks_C)]]$sep
rows_C <- do.call(rbind, lapply(FAMS, function(f) {
  b <- blocks_C[[f]]; d <- b$d
  d$fc <- unname(CLASS_COL[d$author_label_class])
  if (f %in% names(EC_COL)) d$fc[d$author_label_class == "EXPLICIT_EC"] <- EC_COL[[f]]
  d$y <- -b$rows; d
}))

# ---- horizontal layout ----
WC <- 4.90                                    # plot width: x 0.10 to 5.00 in
str_it <- function(s, size = PT) {
  pushViewport(viewport(gp = gpar(fontfamily = "Arial", fontsize = size, fontface = "italic")))
  v <- convertWidth(stringWidth(s), "in", valueOnly = TRUE); popViewport(); v
}
base_C <- function(xmax) ggplot() +
  scale_x_continuous(limits = c(0, xmax), breaks = c(0, 50, 100), labels = function(x) paste0(x, "%"), expand = c(0, 0)) +
  scale_y_continuous(limits = c(-TOTAL_C, 0), breaks = rows_C$y, labels = rows_C$author_label_original, expand = c(0, 0)) +
  labs(x = NULL, y = NULL) + coord_cartesian(clip = "off") + theme_fig2() +
  theme(panel.grid = element_blank(), panel.grid.major = element_blank(), axis.line.x = element_blank())
# Panel width in inches (set by the label column), then the matrix: four columns COL_P apart, right-aligned in the
# panel; the bars' 100% is as long as the longest bar plus its called/total still ends MGAP before the first dot column.
PSPAN_C <- panel_span(ggplotGrob(base_C(100) + geom_blank(data = rows_C, aes(50, y))), WC); WP <- diff(PSPAN_C)
C_MAX <- 2.7                                                  # ggplot point size of a 100% dot
dot_in <- (C_MAX * .pt + .15 * .stroke / 2) * .75 / 72        # its drawn diameter in inches (outline included)
stopifnot(dot_in < P_IN - .012)                               # the largest dot fits inside one row
COL_P <- .18; MGAP <- .16; ANG <- 45
gw <- str_it(GENES_C); gh <- fh("E", PT)
stopifnot(COL_P * sinpi(ANG / 180) > gh + .03)                # neighboring angled names stay apart
cc <- (WP - dot_in / 2 - .012) - (3:0) * COL_P               # column centers, inches from the panel's left
bar_zone <- cc[1] - dot_in / 2 - MGAP                         # bars and numbers end here
tw <- vapply(rows_C$label, str_in, 0); TXT_GAP <- .04
BW <- min((bar_zone - TXT_GAP - tw[rows_C$rate > 0]) / (rows_C$rate[rows_C$rate > 0] / 100))   # inches for 100%
stopifnot(all(rows_C$rate / 100 * BW + TXT_GAP + tw <= bar_zone + 1e-9), BW >= 2.0)
to_x <- function(inches) 100 * inches / BW
XMAX_C <- to_x(WP)

# dot color: mean log expression among the cells with a read, 3 to 8 (above 8 drawn at 8)
FILL_LIM <- c(3, 8)
FILL_RAMP <- grDevices::colorRampPalette(c("#FEE5D9", "#FB6A4A", "#99000D"))(240)
fill_of <- function(v) FILL_RAMP[pmax(1, pmin(240, 1 + round(239 * (pmin(pmax(v, FILL_LIM[1]), FILL_LIM[2]) - FILL_LIM[1]) / diff(FILL_LIM))))]
mm_C$x <- to_x(cc[match(mm_C$gene, GENES_C)])
mm_C$y <- rows_C$y[match(paste(mm_C$family, mm_C$author_label_original), paste(rows_C$family, rows_C$author_label_original))]
mm_C$size <- C_MAX * sqrt(mm_C$share_detecting)               # dot area proportional to the share, as in D
dots_C <- mm_C[mm_C$cells_detecting > 0, ]                    # no dot where no cell has a read
dots_C$fill <- fill_of(dots_C$mean_log_detecting)

grid_C <- do.call(rbind, lapply(blocks_C, function(b) data.frame(x = c(0, 50, 100), y0 = -(min(b$rows) - P_IN / 2), y1 = -b$sep)))
sep_C <- data.frame(y = -vapply(blocks_C, `[[`, 0, "sep"))
title_C <- data.frame(y = -vapply(blocks_C, `[[`, 0, "title"), label = paste(STUDY[FAMS], CONTEXT[FAMS], sep = " | "))
plot_C <- base_C(XMAX_C) +
  geom_segment(data = grid_C, aes(x = x, xend = x, y = y0, yend = y1), color = "#E4E4E4", linewidth = .2) +
  geom_rect(data = rows_C, aes(xmin = 0, xmax = rate, ymin = y - .34 * P_IN, ymax = y + .34 * P_IN), fill = rows_C$fc) +
  geom_text(data = rows_C, aes(x = rate + to_x(TXT_GAP), y = y, label = label), hjust = 0, size = TM, color = "black", family = "Arial") +
  geom_segment(data = sep_C, aes(x = 0, xend = XMAX_C, y = y, yend = y), color = "black", linewidth = .35) +
  geom_text(data = title_C, aes(x = 0, y = y, label = label), hjust = 0, vjust = 0, size = 8 / .pt, fontface = "bold",
            color = "black", family = "Arial") +
  geom_point(data = dots_C, aes(x = x, y = y), shape = 21, color = "black", stroke = .15, size = dots_C$size, fill = dots_C$fill)
grob_C <- ggplotGrob(plot_C)
stopifnot(abs(diff(panel_span(grob_C, WC)) - WP) < 1e-6)
hts <- convertHeight(grob_C$heights, "in", valueOnly = TRUE)   # the panel's own (null) height converts to 0
p_t <- grob_C$layout$t[grob_C$layout$name == "panel"]; p_b <- grob_C$layout$b[grob_C$layout$name == "panel"]
TOP_FIX <- sum(hts[seq_len(p_t - 1)]); BOT_FIX <- sum(hts[-seq_len(p_b)])
HC <- TOP_FIX + TOTAL_C + BOT_FIX

# ---- placement on the page ----
Y0_C <- TOP_C + CAP8 - TOP_IN - TOP_FIX                       # top of the plot
RULE_C <- Y0_C + TOP_FIX + TOTAL_C                            # the last study's bottom rule
LAB_DY <- .055                                                # rule to the gene names' end points
G_X <- .10 + PSPAN_C[1] + cc                                  # gene names' end points (upper right), x
XT_BASE <- Y0_C + HC + .075                                   # x-axis title baseline, 0.075 in under the plot
BAR_MID <- .10 + PSPAN_C[1] + BW / 2                          # center of the bars (0 to 100%) and of the x-axis title
# Key lines: line 1 0.12 in under the x-axis title's baseline, line 2 0.145 in under line 1.
KEY1 <- XT_BASE + .12; KEY2 <- KEY1 + .145
KEY1_LABS <- c("EC (primary study label)", "Other terminal EEC", "Unresolved EEC"); KEY1_GAP <- .12
k1w <- vapply(KEY1_LABS, str_in, 0) + .14; k1w[1] <- k1w[1] + .22
w1 <- sum(k1w) + KEY1_GAP * 2                                 # line 1's width
K2_T1 <- "Cells with a read:"; K2_T2 <- "Their mean log expression:"; K2_SZ <- c(.25, .5, 1); K2_LB <- c("25%", "50%", "100%")
K2_RW <- .45                                                  # color ramp width
k2w_size <- str_in(K2_T1) + .08 + sum(vapply(K2_LB, str_in, 0) + .22) - .07
w2 <- k2w_size + .25 + str_in(K2_T2) + .08 + str_in("3") + .04 + K2_RW + .04 + str_in("8+")   # line 2's width
# The angled gene names, as sampled points of their 45-degree rotated boxes; line 1's box (top of the text to the
# bottom of the color squares) at center c; the nearest distance between them.
name_pts <- do.call(rbind, lapply(seq_along(GENES_C), function(i) {
  uv <- expand.grid(u = seq(-gw[i], 0, length.out = 120), v = seq(-gh / 2, gh / 2, length.out = 9))
  data.frame(x = G_X[i] + uv$u * cospi(ANG / 180) - uv$v * sinpi(ANG / 180),
             y = RULE_C + LAB_DY - uv$u * sinpi(ANG / 180) - uv$v * cospi(ANG / 180))
}))
clear_1 <- function(c) {
  dx <- pmax(c - w1 / 2 - name_pts$x, 0, name_pts$x - (c + w1 / 2)); dy <- pmax(KEY1 - .005 - name_pts$y, 0, name_pts$y - (KEY1 + .095))
  min(sqrt(dx^2 + dy^2))
}
# Panel D's size key starts at D_KEY_L (centered in x 5.15-7.39 in), on about line 2's line.
D_KEY_L <- 5.15 + (2.24 - (str_in("EC cells detecting it:") + .08 + sum(vapply(K2_LB, str_in, 0) + .22) - .07)) / 2
# Both key lines share one center: the largest one, up to the bars' center, that keeps line 1 at least 0.03 in from
# every gene name and line 2 at least 0.45 in left of D's size key; line 2 must still start inside C (x >= 0.10).
cands <- seq(BAR_MID, .10 + w2 / 2, by = -.001)
KEY_MID <- cands[which(vapply(cands, clear_1, 0) >= .03 & cands + w2 / 2 <= D_KEY_L - .45)[1]]
stopifnot(!is.na(KEY_MID))
K1 <- KEY_MID + c(-1, 1) * w1 / 2; K2 <- KEY_MID + c(-1, 1) * w2 / 2   # the lines' left and right ends
# Line 2 lies wholly below the gene names, and both lines inside C's width.
stopifnot(KEY2 > max(name_pts$y) + .02, K1[1] >= .10, K2[1] >= .10, K1[2] <= 5.00)

draw_panel <- function() {
  draw_at(grob_C, .10, Y0_C, WC, HC)
  # the gene names once, under the last study: italic, 45 degrees, each name ending under its column
  grid.text(GENES_C, x = unit(G_X, "in"), y = unit(H - RULE_C - LAB_DY, "in"), rot = ANG, hjust = 1, vjust = .5,
            gp = gpar(fontfamily = "Arial", fontsize = PT, fontface = "italic", col = "black"))
  draw_xtitle(XTITLE, BAR_MID, XT_BASE)
  # key line 1: the label classes; the three EC colors share one label
  x <- K1[1]; y <- KEY1
  ec <- c(CLASS_COL[["EXPLICIT_EC"]], EC_COL[["gse238274"]], EC_COL[["hio"]])
  for (j in 1:3) grid.rect(x = unit(x + (j - 1) * .11, "in"), y = unit(H - y - .045, "in"), width = unit(.095, "in"),
                           height = unit(.095, "in"), just = "left", gp = gpar(fill = ec[j], col = NA))
  text_at(KEY1_LABS[1], x + .36, y - .005); x <- x + k1w[1] + KEY1_GAP
  for (i in 2:3) { key(KEY1_LABS[i], CLASS_COL[i], x, y); x <- x + k1w[i] + KEY1_GAP }
  # key line 2: dot size (dark neutral dots) and dot color (the ramp, as a 240-pixel raster stretched across the bar)
  x <- K2[1]; y <- KEY2
  text_at(K2_T1, x, y); x <- x + str_in(K2_T1) + .08
  for (i in 1:3) {
    grid.points(unit(x + .06, "in"), unit(H - y - .05, "in"), pch = 21,
                gp = gpar(col = "black", fill = "#333333", lwd = .15 * .stroke / 2, fontsize = C_MAX * sqrt(K2_SZ[i]) * .pt + .15 * .stroke / 2))
    text_at(K2_LB[i], x + .15, y); x <- x + str_in(K2_LB[i]) + .22
  }
  x <- x - .07 + .25
  text_at(K2_T2, x, y); x <- x + str_in(K2_T2) + .08
  text_at("3", x, y); x <- x + str_in("3") + .04
  grid.raster(matrix(FILL_RAMP, nrow = 1), x = unit(x, "in"), y = unit(H - y - .01, "in"), width = unit(K2_RW, "in"),
              height = unit(.08, "in"), just = c("left", "top"), interpolate = FALSE)
  text_at("8+", x + K2_RW + .04, y)
}
