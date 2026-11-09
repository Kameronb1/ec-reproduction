# Figure S14A: mature endocrine composition of every SC-islet sample (primary study labels): one bar per sample, stacked
# SC-SIC, SC-beta, SC-alpha, SC-delta, other hormone cells (share of the sample's mature endocrine cells), with its SC-SIC
# share and mature endocrine cell count at the right. Studies with grafts first (graft rows below a dotted rule under
# their in vitro rows), then three studies without grafts; samples with fewer than 50 mature endocrine cells pale.

PH <- 708 / 72                                                 # layout height (in); the strip shows its lower 5792 px
S1 <- read_table("s1_composition.tsv")
# The drawn rows in order: study block, block name and line, and row label (drawing-only, s1_rows.tsv); each sample's
# mature endocrine cells, arm and SC-SIC share from s1_composition.tsv (no SC-SIC row: 0%).
D <- read_drawing("s1_rows.tsv")
D$key <- paste(D$study, D$sample); S1$key <- paste(S1$study, S1$sample)
k <- match(D$key, S1$key); stopifnot(!anyNA(k))
D$mature <- S1$mature[k]; D$graft <- S1$arm[k] == "graft"
sic <- S1[S1$category == "SC-SIC", ]; D$sic_pct <- sic$pct[match(D$key, sic$key)]; D$sic_pct[is.na(D$sic_pct)] <- 0
D$sic_label <- as.character(round(D$sic_pct)); D$n_label <- formatC(D$mature, big.mark = ",", format = "d")
D$pale <- D$mature < 50

X_AXIS <- "Share of mature endocrine cells (%, primary study labels)"
STACK <- c("SC-SIC", "SC-beta", "SC-alpha", "SC-delta", "other hormone")                    # bar and key order
COL <- c(`SC-SIC` = "#AA3377", `SC-beta` = "#4477AA", `SC-alpha` = "#117733", `SC-delta` = "#AE8C44", `other hormone` = "gray70")
PALE_COL <- local({ m <- col2rgb(COL); setNames(rgb(t(round(0.4 * m + 0.6 * 255)), maxColorValue = 255), names(COL)) })   # 40% color on white
KEY_LAB <- c(`SC-SIC` = "SC-SIC", `SC-beta` = "SC-beta", `SC-alpha` = "SC-alpha", `SC-delta` = "SC-delta", `other hormone` = "other hormone (epsilon, gamma)")
PALE_KEY <- "< 50 mature endocrine cells (pale bar)"
BAND_FILL <- "#E8E8E8"

# ---- geometry (inches from the layout's top-left) ----
RH <- 0.103; GAP <- 0.045; GAP_VERES <- 0.06; BAR_FRAC <- 0.86      # row pitch, gap between studies (more before Veres 2019), bar / row
# graft rows sit SEPG lower than their study's in vitro rows, with the dotted rule RULE_OFF below the top of that gap;
# row text sits TS above its centered baseline
SEPG <- 0.012; RULE_OFF <- 0.0092; TS <- 0.005
Y0 <- 0.88 - 5 * SEPG                                              # top of the first row (five studies have graft rows)
XB0 <- 3.46; XB1 <- 6.46; XU <- (XB1 - XB0) / 100                  # bars: 0-100% at 3.46-6.46 in
X_STUDY <- 1.40; X_ROW <- 3.41; X_SIC <- 6.84; X_N <- 7.40; X_ROW_L <- 1.46; X_STUDY_L <- 0.10; BOX_PAD <- 0.015
HEAD_BASE <- c(0.62, 0.73, 0.84) - (0.88 - Y0)
CHIP_W <- 0.13; CHIP_H <- 0.09; CHIP_GAP <- 0.05; ENTRY_GAP <- 0.16; PALE_W <- 0.20; KEY_MID <- c(9.60, 9.745)

# rows, block by block (in vitro rows first, then the graft rows SEPG lower)
blocks <- unique(D$block); D$y_top <- NA_real_; y <- Y0
for (i in seq_along(blocks)) {
  ii <- which(D$block == blocks[i])
  if (i > 1) y <- y + GAP + (if (blocks[i] == "Veres2019") GAP_VERES else 0)
  sep <- any(D$graft[ii]) * SEPG
  D$y_top[ii] <- y + (seq_along(ii) - 1) * RH + D$graft[ii] * sep
  y <- y + length(ii) * RH + sep
}
Y_BARS_END <- y
D$y_mid <- D$y_top + RH / 2; D$base_y <- D$y_mid + CAP7 / 2 - TS
# study blocks: bold name and plain line on the baselines of the first two rows (a one-row block centers them on its row);
# the dotted rule between in vitro and graft rows
BLK <- do.call(rbind, lapply(blocks, function(b) { d <- D[D$block == b, ]; n <- nrow(d); g <- which(d$graft)[1]
  data.frame(block = b, name = d$block_name[1], line = d$block_line[1], top = min(d$y_top), bottom = max(d$y_top) + RH,
             first_base = if (n > 1) d$base_y[1] else d$base_y[1] - RH / 2, second_base = if (n > 1) d$base_y[2] else d$base_y[1] + RH / 2,
             rule_y = if (is.na(g)) NA_real_ else d$y_top[g] - SEPG + RULE_OFF) }))
s <- BLK$block == "Siehler2024"                                    # its centered label sits 15 px lower, clear of Zhu 2023's band
BLK$first_base[s] <- BLK$first_base[s] + 15 / 600; BLK$second_base[s] <- BLK$second_base[s] + 15 / 600
# light grey bands behind study blocks 1, 3, 5 and 7, whole pixels, 0.012 in beyond the block's rows
BAND <- BLK[seq(1, nrow(BLK), by = 2), ]
BAND$y0 <- round((BAND$top - 0.012) * 600) / 600; BAND$y1 <- round((BAND$bottom + 0.012) * 600) / 600

# bar segments: share of mature endocrine cells, stacked left to right in STACK order
SEG <- S1[S1$key %in% D$key, c("key", "category", "cells")]
SEG <- SEG[order(match(SEG$key, D$key), match(SEG$category, STACK)), ]
SEG$pct <- 100 * SEG$cells / D$mature[match(SEG$key, D$key)]; SEG$pale <- D$pale[match(SEG$key, D$key)]
cum <- ave(SEG$pct, SEG$key, FUN = cumsum)
SEG$x0 <- XB0 + XU * (cum - SEG$pct); SEG$x1 <- XB0 + XU * cum; SEG$y_mid <- D$y_mid[match(SEG$key, D$key)]

# key: two rows left-aligned with the bars; each entry is a chip and its label, placed by the label's measured width
key <- list()
for (r in 1:2) { x <- XB0
  for (k in list(STACK[1:4], c(STACK[5], "pale"))[[r]]) {
    cw <- if (k == "pale") PALE_W else CHIP_W
    key[[length(key) + 1]] <- list(k = k, row = r, x = x, cw = cw, tx = x + cw + CHIP_GAP)
    x <- x + cw + CHIP_GAP + text_w(if (k == "pale") PALE_KEY else KEY_LAB[[k]], 7) + ENTRY_GAP } }

draw_panel <- function() in_strip(0, 5792, function() {
  XI <- function(x) unit(x, "in"); YI <- function(y) unit(PH - y, "in")
  txt <- function(s, x, y_base, size, bold = FALSE, hjust = 0)
    grid.text(s, x = XI(x), y = YI(y_base), hjust = hjust, vjust = 0, gp = gpar(fontfamily = FONT, fontsize = size, fontface = if (bold) "bold" else "plain", col = "#000000"))
  rect <- function(x0, y0, x1, y1, fill, col = NA, lwd = 1)
    grid.rect(x = XI(x0), y = YI(y1), width = unit(x1 - x0, "in"), height = unit(y1 - y0, "in"), just = c("left", "bottom"),
              gp = gpar(fill = fill, col = col, lwd = lwd, alpha = 1, lty = "solid"))
  rect(X_STUDY_L, BAND$y0, X_N, BAND$y1, BAND_FILL)
  # column headers
  txt(c("SC-SIC", "(%)"), X_SIC, HEAD_BASE[2:3], 7, TRUE, hjust = 1)
  txt(c("Mature", "endocrine", "cells"), X_N, HEAD_BASE, 7, TRUE, hjust = 1)
  # dotted in vitro / graft rules, then bars (full color, then the pale < 50-cell samples) and the study boxes
  r <- BLK[!is.na(BLK$rule_y), ]
  grid.segments(XI(X_ROW_L), YI(r$rule_y), XI(X_N), YI(r$rule_y), gp = gpar(col = "gray60", lwd = 0.85, lty = "dotted", lineend = "butt"))
  hb <- RH * BAR_FRAC / 2
  for (pl in c(FALSE, TRUE)) { g <- SEG[SEG$pale == pl, ]
    rect(g$x0, g$y_mid - hb, g$x1, g$y_mid + hb, fill = if (pl) PALE_COL[g$category] else COL[g$category], col = "white", lwd = 0.3) }
  rect(XB0 - BOX_PAD, BLK$top, XB1 + BOX_PAD, BLK$bottom, fill = NA, col = "gray75", lwd = 0.8)
  # study labels, row labels, SC-SIC (%) and mature endocrine cell counts
  txt(BLK$name, X_STUDY, BLK$first_base, 7, TRUE, hjust = 1); txt(BLK$line, X_STUDY, BLK$second_base, 7, hjust = 1)
  txt(D$drawn_label, X_ROW, D$base_y, 7, hjust = 1)
  txt(D$sic_label, X_SIC, D$base_y, 7, hjust = 1)
  txt(D$n_label, X_N, D$base_y, 7, hjust = 1)
  # x axis
  tk <- c(0, 25, 50, 75, 100)
  grid.segments(XI(XB0 + XU * tk), YI(Y_BARS_END), XI(XB0 + XU * tk), YI(Y_BARS_END + 0.025), gp = gpar(col = "gray40", lwd = 0.8))
  txt(as.character(tk), XB0 + XU * tk, Y_BARS_END + 0.025 + 0.025 + CAP7, 7, hjust = 0.5)
  txt(X_AXIS, (XB0 + XB1) / 2, 9.43, 7, hjust = 0.5)
  # key
  for (it in key) { ym <- KEY_MID[it$row]
    if (it$k == "pale") {
      pk <- c("SC-SIC", "SC-beta", "SC-alpha"); xs <- it$x + it$cw * (0:3) / 3
      rect(xs[1:3], ym - CHIP_H / 2, xs[2:4], ym + CHIP_H / 2, fill = PALE_COL[pk], col = "white", lwd = 0.3)
      txt(PALE_KEY, it$tx, ym + CAP7 / 2, 7)
    } else {
      rect(it$x, ym - CHIP_H / 2, it$x + it$cw, ym + CHIP_H / 2, fill = COL[[it$k]])
      txt(KEY_LAB[[it$k]], it$tx, ym + CAP7 / 2, 7)
    } }
})
