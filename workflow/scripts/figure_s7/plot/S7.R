# Supplementary Figure S7: mature endocrine composition of every sequenced SC-islet sample, in vitro and graft
# (published labels). One stacked bar per sample (SC-SIC, SC-beta, SC-alpha, SC-delta, other hormone cells; pale when
# the sample has fewer than 50 mature endocrine cells), grouped by study and cell line, with the SC-SIC share and the
# number of mature endocrine cells at the right. Values from s1_composition.tsv (Figure 4 unit); row order and row
# labels from s7_rows.tsv in this folder.

S <- read_table("s1_composition.tsv", "figure4")             # one row per sample and cell type
OLD <- read.delim(file.path(dirname(PLOT_SCRIPT), "s7_rows.tsv"), stringsAsFactors = FALSE, check.names = FALSE, encoding = "UTF-8")
STACK <- c("SC-SIC", "SC-beta", "SC-alpha", "SC-delta", "other hormone")
stopifnot(all(S$category %in% STACK))
S$key <- paste(S$study, S$sample); OLD$key <- paste(OLD$study, OLD$sample)
stopifnot(setequal(unique(S$key), OLD$key), !anyDuplicated(OLD$key))

# study order: the four graft studies in Figure 4A order, then the in vitro only studies (Zhu 2023, Siehler 2024), then
# Veres 2019
BLOCKS <- c("A20_HUES8", "A20_WS4corr", "Balboa2022", "Du2022", "A23", "Zhu2023", "Siehler2024", "Veres2019")
stopifnot(setequal(BLOCKS, unique(OLD$block)))
D <- OLD[order(match(OLD$block, BLOCKS), OLD$row), c("key", "study", "sample", "context", "graft", "block", "block_name", "block_line", "drawn_label")]
D$graft <- as.logical(D$graft)
# within a block, in vitro rows first, then grafts
D <- D[order(match(D$block, BLOCKS), D$graft, seq_len(nrow(D))), ]
# the three extended-culture end points of Figure 4 are named as such
EXT <- c("Augsornworawat2023 s6_month6" = "stage 6 mo 6, extended culture", "Augsornworawat2023 s6_month12" = "stage 6 mo 12, extended culture",
         "Balboa2022 SRR13805898" = "stage 7 wk 6 #1, extended culture", "Balboa2022 SRR13805899" = "stage 7 wk 6 #2, extended culture")
stopifnot(all(names(EXT) %in% D$key))
D$label <- D$drawn_label; D$label[match(names(EXT), D$key)] <- unname(EXT)

# values from the composition table: mature endocrine cells, cells per type, share
D$mature <- sapply(D$key, function(k) { m <- unique(S$mature[S$key == k]); stopifnot(length(m) == 1); m })
D$cells_sum <- sapply(D$key, function(k) sum(S$cells[S$key == k]))
stopifnot(all(D$cells_sum == D$mature))
for (ct in STACK) D[[ct]] <- sapply(D$key, function(k) sum(S$cells[S$key == k & S$category == ct]))
D$sic_pct <- 100 * D$`SC-SIC` / D$mature
D$pale <- D$mature < 50
D$sic_label <- sprintf("%.0f", D$sic_pct); D$n_label <- formatC(D$mature, format = "d", big.mark = ",")

# ---- geometry (inches from the page top-left) ---------------------------------------------------------------------------
RH <- 0.102; GAP <- 0.038; GAP_SECTION <- 0.05; BAR_FRAC <- 0.80; SEPG <- 0.012; RULE_OFF <- 0.0092
CAP7 <- 0.716 * PT / 72
Y0 <- 0.46; HEAD_BASE <- c(0.19, 0.30, 0.41)
XB0 <- 3.46; XB1 <- 6.46; XU <- (XB1 - XB0) / 100
X_STUDY <- 1.40; X_ROW <- 3.41; X_ROW_L <- 1.46; X_STUDY_L <- 0.10; X_SIC <- 6.84; X_N <- 7.40; BOX_PAD <- 0.015
D$y_top <- NA_real_; y <- Y0; BL <- list()
for (i in seq_along(BLOCKS)) {
  ii <- which(D$block == BLOCKS[i])
  if (i > 1) y <- y + GAP + (if (BLOCKS[i] == "Zhu2023") GAP_SECTION else 0)
  sep <- any(D$graft[ii]) * SEPG
  slots <- max(length(ii), 2)                     # a one-sample block keeps two lines for its study name and line
  off <- (slots - length(ii)) * RH / 2            # its bar sits in the middle
  D$y_top[ii] <- y + off + (seq_along(ii) - 1) * RH + D$graft[ii] * sep
  top <- y; y <- y + slots * RH + sep
  BL[[i]] <- data.frame(block = BLOCKS[i], top = top, bottom = y)
}
BLK <- do.call(rbind, BL); Y_BARS_END <- y
D$y_mid <- D$y_top + RH / 2; D$base_y <- D$y_mid + CAP7 / 2 - 0.005
BLK$name <- D$block_name[match(BLK$block, D$block)]; BLK$line <- D$block_line[match(BLK$block, D$block)]
BLK$first_base <- BLK$top + RH / 2 + CAP7 / 2 - 0.005; BLK$second_base <- BLK$first_base + RH
BLK$rule_y <- sapply(BLK$block, function(b) { d <- D[D$block == b, ]; g <- which(d$graft)[1]; if (is.na(g)) NA_real_ else d$y_top[g] - SEPG + RULE_OFF })
BAND <- BLK[seq(1, nrow(BLK), by = 2), ]
BAND$y0 <- round((BAND$top - 0.012) * DPI) / DPI; BAND$y1 <- round((BAND$bottom + 0.012) * DPI) / DPI

SEG <- do.call(rbind, lapply(seq_len(nrow(D)), function(i) {
  p <- 100 * unlist(D[i, STACK]) / D$mature[i]; cum <- cumsum(p)
  data.frame(key = D$key[i], category = STACK, cells = unlist(D[i, STACK]), pct = p, x0 = XB0 + XU * (cum - p), x1 = XB0 + XU * cum,
             y_mid = D$y_mid[i], pale = D$pale[i], row.names = NULL) }))
SEG <- SEG[SEG$cells > 0, ]

COL <- c(`SC-SIC` = "#AA3377", `SC-beta` = "#4477AA", `SC-alpha` = "#117733", `SC-delta` = "#AE8C44", `other hormone` = "grey70")
PALE_COL <- local({ m <- col2rgb(COL); setNames(rgb(t(round(0.4 * m + 0.6 * 255)), maxColorValue = 255), names(COL)) })
KEY_LAB <- c(`SC-SIC` = "SC-SIC", `SC-beta` = "SC-beta", `SC-alpha` = "SC-alpha", `SC-delta` = "SC-delta", `other hormone` = "other hormone (epsilon, gamma)")
PALE_KEY <- "under 50 mature endocrine cells (pale bar)"
X_AXIS <- "Share of mature endocrine cells (%, published labels)"
BAND_FILL <- "#E8E8E8"

# text width on a 600-dpi cairo device
text_w <- local({ dev <- NULL; function(s, face = "plain") { cur <- dev.cur()
  if (is.null(dev) || !(dev %in% dev.list())) { png(tempfile(fileext = ".png"), width = 2, height = 2, units = "in", res = DPI, type = "cairo"); dev <<- dev.cur() }
  dev.set(dev); w <- vapply(s, function(x) convertWidth(grobWidth(textGrob(x, gp = gpar(fontfamily = FONT, fontsize = PT, fontface = face))), "in", valueOnly = TRUE), 0)
  if (cur > 1L) dev.set(cur); unname(w) } })

# key: one line from the left edge of the bars if it fits inside the page, otherwise two lines
CHIP_W <- 0.13; CHIP_H <- 0.08; CHIP_GAP <- 0.05; ENTRY_GAP <- 0.16; PALE_W <- 0.20
ents <- c(STACK, "pale"); ew <- ifelse(ents == "pale", PALE_W, CHIP_W) + CHIP_GAP + text_w(ifelse(ents == "pale", PALE_KEY, KEY_LAB[ents]))
TICK_LAB_BASE <- Y_BARS_END + 0.025 + 0.025 + CAP7
AXIS_BASE <- TICK_LAB_BASE + 0.14
KEY_X0 <- X_ROW_L
one_line <- KEY_X0 + sum(ew) + ENTRY_GAP * (length(ents) - 1) <= X_N
KEY_ROWS <- if (one_line) list(ents) else list(STACK[1:4], c(STACK[5], "pale"))
KEY_MID <- AXIS_BASE + 0.15 + (seq_along(KEY_ROWS) - 1) * 0.15
key <- list(); for (r in seq_along(KEY_ROWS)) { x <- KEY_X0
  for (k in KEY_ROWS[[r]]) { cw <- if (k == "pale") PALE_W else CHIP_W
    key[[length(key) + 1]] <- list(k = k, row = r, x = x, cw = cw, tx = x + cw + CHIP_GAP)
    x <- x + cw + CHIP_GAP + text_w(if (k == "pale") PALE_KEY else KEY_LAB[[k]]) + ENTRY_GAP } }

# ---- drawing --------------------------------------------------------------------------------------------------------------
draw_panel <- function() {
  XI <- function(x) unit(x, "in"); YI <- function(y) unit(H - y, "in")
  txt <- function(s, x, y_base, bold = FALSE, hjust = 0)
    grid.text(s, x = XI(x), y = YI(y_base), hjust = hjust, vjust = 0, gp = gpar(fontfamily = FONT, fontsize = PT, fontface = if (bold) "bold" else "plain", col = "#000000"))
  rect <- function(x0, y0, x1, y1, fill, col = NA, lwd = 1)
    grid.rect(x = XI(x0), y = YI(y1), width = unit(x1 - x0, "in"), height = unit(y1 - y0, "in"), just = c("left", "bottom"), gp = gpar(fill = fill, col = col, lwd = lwd))
  grid.rect(gp = gpar(fill = "white", col = NA))
  rect(X_STUDY_L, BAND$y0, X_N, BAND$y1, BAND_FILL)
  txt(c("SC-SIC", "(%)"), X_SIC, HEAD_BASE[2:3], TRUE, hjust = 1)
  txt(c("Mature", "endocrine", "cells"), X_N, HEAD_BASE, TRUE, hjust = 1)
  r <- BLK[!is.na(BLK$rule_y), ]
  grid.segments(XI(X_ROW_L), YI(r$rule_y), XI(X_N), YI(r$rule_y), gp = gpar(col = "grey60", lwd = 0.85, lty = "dotted", lineend = "butt"))
  hb <- RH * BAR_FRAC / 2
  for (pl in c(FALSE, TRUE)) { g <- SEG[SEG$pale == pl, ]
    if (nrow(g)) rect(g$x0, g$y_mid - hb, g$x1, g$y_mid + hb, fill = if (pl) PALE_COL[g$category] else COL[g$category], col = "white", lwd = 0.3) }
  for (i in seq_len(nrow(BLK))) { d <- D[D$block == BLK$block[i], ]
    rect(XB0 - BOX_PAD, min(d$y_top), XB1 + BOX_PAD, max(d$y_top) + RH, fill = NA, col = "grey75", lwd = 0.8) }
  txt(BLK$name, X_STUDY, BLK$first_base, TRUE, hjust = 1); txt(BLK$line, X_STUDY, BLK$second_base, hjust = 1)
  txt(D$label, X_ROW, D$base_y, hjust = 1)
  txt(D$sic_label, X_SIC, D$base_y, hjust = 1)
  txt(D$n_label, X_N, D$base_y, hjust = 1)
  tk <- c(0, 25, 50, 75, 100)
  grid.segments(XI(XB0 + XU * tk), YI(Y_BARS_END), XI(XB0 + XU * tk), YI(Y_BARS_END + 0.025), gp = gpar(col = "grey40", lwd = 0.8))
  txt(as.character(tk), XB0 + XU * tk, TICK_LAB_BASE, hjust = 0.5)
  txt(X_AXIS, (XB0 + XB1) / 2, AXIS_BASE, hjust = 0.5)
  for (it in key) { ym <- KEY_MID[it$row]
    if (it$k == "pale") { pk <- c("SC-SIC", "SC-beta", "SC-alpha"); xs <- it$x + it$cw * (0:3) / 3
      rect(xs[1:3], ym - CHIP_H / 2, xs[2:4], ym + CHIP_H / 2, fill = PALE_COL[pk], col = "white", lwd = 0.3); txt(PALE_KEY, it$tx, ym + CAP7 / 2)
    } else { rect(it$x, ym - CHIP_H / 2, it$x + it$cw, ym + CHIP_H / 2, fill = COL[[it$k]]); txt(KEY_LAB[[it$k]], it$tx, ym + CAP7 / 2) } }
}
