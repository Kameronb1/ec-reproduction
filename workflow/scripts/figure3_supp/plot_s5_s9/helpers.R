# Drawing helpers shared by the Figure S9-S13 panels: text size, colors, the plot theme, placement on the page, the
# gene-group heatmaps of S10 and S11, the term panels of S13 and the groups of S12. Sourced by
# workflow/scripts/common/draw_panel.R, which defines the page size W x H (inches) and DPI. Places on the page are given
# in whole pixels from its top-left corner.

FONT <- "Arial"; PT <- 7.5
PW <- 7.3                                               # width of a full-width panel (in)
GRID_COL <- "#E4E4E4"
GROUP_COL <- c(G1 = "#43A5B6", G2 = "#A89A2E", G3 = "#B4553F")       # gene groups, as in Figure 3F and 3G

# The ggplot look shared by every plotted panel.
panel_theme <- function(grid = c("x", "y", "both", "none"), legend = "none") {
  grid <- match.arg(grid)
  gl <- element_line(color = GRID_COL, linewidth = .2)
  theme_minimal(base_family = FONT, base_size = PT) +
    theme(text = element_text(color = "black", size = PT),
          axis.text = element_text(color = "black", size = PT),
          axis.title = element_text(color = "black", size = PT),
          strip.text = element_text(color = "black", size = 8, face = "bold", hjust = 0),
          plot.title = element_text(color = "black", size = 8, face = "bold"),
          legend.text = element_text(color = "black", size = PT),
          legend.title = element_text(color = "black", size = PT),
          legend.position = legend, legend.key.size = unit(.12, "in"),
          panel.grid.minor = element_blank(),
          panel.grid.major.x = if (grid %in% c("x", "both")) gl else element_blank(),
          panel.grid.major.y = if (grid %in% c("y", "both")) gl else element_blank(),
          axis.line.x = element_line(color = "black", linewidth = .35),
          axis.ticks.x = element_line(color = "black", linewidth = .3),
          axis.ticks.length = unit(1.3, "pt"),
          panel.spacing = unit(.16, "in"),
          plot.margin = margin(2, 6, 2, 2))
}

# A plotted panel of nominal height h (inches): its canvas is h - 0.41 in tall (whole pixels, rounded down) with its
# top-left corner at page pixel (x, y), and the plot fills the bottom h - 0.47 in of it, w inches wide. plot is a ggplot
# or patchwork, or a function that draws.
plot_panel <- function(plot, x, y, w, h) {
  h_px <- floor((h - .41) * DPI + 1e-6)
  pushViewport(viewport(x = unit(x / DPI, "in"), y = unit((H * DPI - y - h_px) / DPI, "in"), width = unit(w, "in"),
                        height = unit(h - (.30 + .17), "in"), just = c("left", "bottom")))
  if (is.function(plot)) plot() else print(plot, newpage = FALSE)
  popViewport()
}

# Signed axis labels with a true minus sign.
mlab <- function(b) vapply(b, function(z) if (is.na(z)) NA_character_ else gsub("-", "−", format(z, trim = TRUE, drop0trailing = TRUE)), "")

# ---- Gene-group heatmaps (S10, S11) -------------------------------------------------------------------------------------
# Drawn with grid at fixed inch positions inside a canvas PW inches wide: x from its left edge, y down from its top edge.
# The canvas, with its top-left corner at page pixel (x, y) and h inches tall, is drawn on a device of its own and placed
# on the page pixel for pixel: its color ramps are images stretched with interpolation, and the interpolated pixels
# depend on where on the device the image falls.
heatmap_canvas <- function(x, y, h, draw) {
  w_px <- floor(PW * DPI + 1e-6); h_px <- floor(h * DPI + 1e-6)
  f <- tempfile(fileext = ".png")
  png(f, width = w_px, height = h_px, units = "px", res = DPI, type = "cairo", bg = "white")
  grid.newpage(); draw(); invisible(dev.off())
  grid.raster(png::readPNG(f, native = TRUE), x = unit(x / DPI, "in"), y = unit((H * DPI - y) / DPI, "in"),
              width = unit(w_px / DPI, "in"), height = unit(h_px / DPI, "in"), just = c("left", "top"), interpolate = FALSE)
  unlink(f)
}

# The two heatmap tables; S10 and S11 share one level color scale. Columns: block, row, gene, then 5 fold changes and
# 9 per-study levels, in drawn column order.
T6 <- function() read_table("s6_gene_groups/g1_g2_heatmap.tsv")
T7 <- function() read_table("s6_gene_groups/g3_heatmap.tsv")
FC <- 4:8; LV <- 9:17
LV_RANGE <- function() range(as.matrix(rbind(T6(), T7())[, LV]))
# drawn columns: 9 levels (3 intestinal, 6 SC-islet studies in Figure 3A order), then 5 fold changes
LV_LAB <- c("Elmentaite 2020 + 2021\n(donors pooled)", "Hickey 2023 (nuclei)", "Egozi 2021", "Augsornworawat 2020", "Balboa 2022",
            "Zhu 2023", "Siehler 2024", "Du 2022", "Augsornworawat 2023\n(nuclei)")
FC_LAB <- c("Elmentaite 2020 + 2021\n(mean of the 2 datasets)", "Hickey 2023", "Egozi 2021",
            "vs other SC endocrine\n(median of 6 studies)", "minus intestinal EC (Fig. 3D)")
COLS <- data.frame(part = c(rep("lv", 9), rep("fc", 5)),
                   grp = c(rep("Intestinal EC", 3), rep("SC-SIC", 6), rep("EC vs other EEC", 3), rep("SC-SIC", 2)),
                   label = c(LV_LAB, FC_LAB), kind = c(rep("lv", 9), "gut", "gut", "gut", "sic", "d"), stringsAsFactors = FALSE)

# colors: levels on a magma ramp; fold changes on white-to-color ramps by where the gene is higher, ending at 6
CAPV <- 6                                               # larger fold changes are drawn at 6 and marked ">"
MAGMA <- viridisLite::magma(256, begin = .08, end = .96); MAGMA_AT <- seq(0, 1, length.out = 256)^1.5
lv_col <- function(v, rng) scales::gradient_n_pal(MAGMA, MAGMA_AT, "Lab")(scales::rescale(v, from = rng))
RAMP <- c(gut = "#43A5B6", sic = "#AA3377", other = "grey35")
ramp_col <- function(m, end) { m <- pmin(abs(m), CAPV) / CAPV
  rgb(grDevices::colorRamp(c("#FFFFFF", end), space = "Lab")(m), maxColorValue = 255) }
fc_meaning <- function(kind, v) ifelse(kind == "gut", ifelse(v >= 0, "gut", "other"),
                               ifelse(kind == "sic", ifelse(v >= 0, "sic", "other"), ifelse(v >= 0, "sic", "gut")))

# one block of rows -> one row per cell (column by column): value, fill, "+" (2-fold up in an intestinal study), ">" (beyond 6)
cells <- function(b, plus, rng) {
  n <- nrow(b)
  do.call(rbind, lapply(seq_len(14), function(j) {
    v <- if (j <= 9) b[[LV[j]]] else b[[FC[j - 9]]]
    k <- COLS$kind[j]
    data.frame(row = seq_len(n), col = j, v = v,
               fill = if (k == "lv") lv_col(v, rng) else vapply(seq_len(n), function(i) ramp_col(v[i], RAMP[[fc_meaning(k, v[i])]]), ""),
               plus = plus & k == "gut" & v >= 1, clip = k != "lv" & abs(v) > CAPV, stringsAsFactors = FALSE)
  }))
}

PT_GENE <- 7; PT_COL <- 7; PT_STRIP <- 7.5; PT_KEY <- 7.5
TILE_LWD <- .2 * 72.27 / 25.4                           # tile border, 0.2 mm
YU <- function(y) unit(1, "npc") - unit(y, "in"); XU <- function(x) unit(x, "in")
gp <- function(pt, face = "plain", ...) gpar(fontfamily = FONT, fontsize = pt, fontface = face, col = "black", ...)
str_in <- function(s, pt, face = "plain") {             # widest line of each string, inches
  pushViewport(viewport(gp = gpar(fontfamily = FONT, fontsize = pt, fontface = face))); on.exit(popViewport())
  vapply(s, function(z) max(vapply(strsplit(z, "\n")[[1]], function(l) convertWidth(stringWidth(l), "in", valueOnly = TRUE), 0)), 0, USE.NAMES = FALSE)
}
LINE_IN <- function(pt) pt * 0.9 / 72                   # line advance at lineheight 0.9
col_x <- function(left, p, g_facet, g_part) {          # column centers: gaps between column groups and between the two parts
  x <- numeric(14); cur <- left
  for (i in 1:14) {
    if (i > 1) cur <- cur + if (COLS$part[i] != COLS$part[i - 1]) g_part else if (COLS$grp[i] != COLS$grp[i - 1]) g_facet else 0
    x[i] <- cur + p / 2; cur <- cur + p
  }
  x
}
COL_OFF <- .05
label_depth <- function(angle) {                        # how far the rotated column labels reach below the rows
  L <- str_in(COLS$label, PT_COL); th <- lengths(strsplit(COLS$label, "\n")) * LINE_IN(PT_COL); a <- angle * pi / 180
  max(COL_OFF + L * sin(a) + th / 2 * cos(a))
}
# group header: color dot and bold name, starting y in below the canvas top; HEADER_H is the height it takes
HEADER_H <- .06 + .145
block_heading <- function(y, chip, header) {
  grid.circle(unit(.04 + .0375, "in"), YU(y + .06 + .0373), r = unit(.0375, "in"), gp = gpar(fill = chip, col = NA))
  grid.text(header, unit(.04 + .115, "in"), YU(y + .06), just = c("left", "top"), gp = gp(PT, "bold"))
}
draw_strips <- function(y_line, xs, p, fc_label) {     # bold column-group headers over a rule
  grp_id <- paste(COLS$part, COLS$grp)
  for (g in unique(grp_id)) {
    i <- which(grp_id == g); lab <- COLS$grp[i[1]]
    if (lab == "EC vs other EEC") lab <- fc_label
    grid.text(lab, XU(mean(xs[i])), YU(y_line - .035), just = c("centre", "bottom"), gp = gp(PT_STRIP, "bold", lineheight = .9))
    grid.lines(XU(c(min(xs[i]) - p / 2 + .02, max(xs[i]) + p / 2 - .02)), YU(c(y_line, y_line)), gp = gpar(col = "black", lwd = .7))
  }
}
draw_block <- function(cc, genes, y_top, xs, p, R, gene_right) {   # tiles, a frame per column group, + and > marks, gene names
  n <- length(genes)
  grid.rect(XU(xs[cc$col]), YU(y_top + (cc$row - .5) * R), width = unit(p, "in"), height = unit(R, "in"),
            gp = gpar(fill = cc$fill, col = "white", lwd = TILE_LWD))
  grp_id <- paste(COLS$part, COLS$grp)
  for (g in unique(grp_id)) { i <- which(grp_id == g)
    grid.rect(XU(min(xs[i]) - p / 2), YU(y_top), width = unit(max(xs[i]) - min(xs[i]) + p, "in"), height = unit(n * R, "in"),
              just = c("left", "top"), gp = gpar(fill = NA, col = "grey75", lwd = .45)) }
  if (any(cc$plus)) { q <- cc[cc$plus, ]
    grid.points(XU(xs[q$col]), YU(y_top + (q$row - .5) * R), pch = 3, size = unit(4.4, "pt"), gp = gpar(col = "grey10", lwd = .9)) }
  if (any(cc$clip)) { q <- cc[cc$clip, ]
    grid.text(">", XU(xs[q$col]), YU(y_top + (q$row - .5) * R), gp = gp(PT_GENE, "bold")) }
  grid.text(genes, XU(gene_right), YU(y_top + (seq_len(n) - .5) * R), just = c("right", "centre"), gp = gp(PT_GENE, "italic"))
}
draw_col_labels <- function(y_bot, xs, angle)
  grid.text(COLS$label, XU(xs), YU(y_bot + COL_OFF), rot = angle, just = c("right", "centre"), gp = gp(PT_COL, lineheight = .9))
KEY_H <- .62; KEY_GAP <- .17; KEY_STEP <- .115
ramp_bar <- function(cols, x, y, w)                     # a color ramp 0.085 in tall
  grid.raster(matrix(cols, nrow = 1), XU(x), YU(y), width = unit(w, "in"), height = unit(.085, "in"), just = c("left", "top"), interpolate = TRUE)
draw_keys <- function(y, x_lv, x_fc, x_sym, plus, rng) {   # level key, three fold-change ramps, symbol entries
  grid.text(expression("Per-study expression, log"[2]*"(1 + CPM)"), XU(x_lv), YU(y), just = c("left", "top"), gp = gp(PT_KEY))
  bw <- 1.3; yb <- y + KEY_GAP
  ramp_bar(lv_col(seq(rng[1], rng[2], length.out = 256), rng), x_lv, yb, bw)
  br <- scales::extended_breaks()(rng); br <- br[br >= rng[1] & br <= rng[2]]
  bx <- x_lv + (br - rng[1]) / diff(rng) * bw
  grid.segments(XU(bx), YU(yb + .085), XU(bx), YU(yb + .085 - .025), gp = gpar(col = "white", lwd = .8))
  grid.text(scales::label_number(style_negative = "minus")(br), XU(bx), YU(yb + .085 + .035), just = c("centre", "top"), gp = gp(PT_KEY))
  grid.text(expression("log"[2]*" fold change"), XU(x_fc), YU(y), just = c("left", "top"), gp = gp(PT_KEY))
  fw <- .9; labs <- c(gut = "Higher in intestinal EC", sic = "Higher in SC-SIC", other = "Higher in the other endocrine cells")
  for (k in seq_along(labs)) {
    yk <- y + KEY_GAP + (k - 1) * KEY_STEP
    ramp_bar(ramp_col(seq(0, CAPV, length.out = 256), RAMP[[names(labs)[k]]]), x_fc, yk, fw)
    grid.rect(XU(x_fc), YU(yk), width = unit(fw, "in"), height = unit(.085, "in"), just = c("left", "top"), gp = gpar(fill = NA, col = "grey60", lwd = .4))
    grid.text(labs[k], XU(x_fc + fw + .07), YU(yk + .0425), just = c("left", "centre"), gp = gp(PT_KEY))
  }
  yt <- y + KEY_GAP + 2 * KEY_STEP + .085
  tx <- x_fc + seq(0, CAPV, 2) / CAPV * fw
  grid.segments(XU(tx), YU(yt), XU(tx), YU(yt + .025), gp = gpar(col = "black", lwd = .5))
  grid.text(seq(0, CAPV, 2), XU(tx), YU(yt + .035), just = c("centre", "top"), gp = gp(PT_KEY))
  ys <- y + KEY_GAP + .0425
  if (plus) {
    grid.points(XU(x_sym + .05), YU(ys), pch = 3, size = unit(4.4, "pt"), gp = gpar(col = "grey10", lwd = .9))
    grid.text("2-fold up in that intestinal study", XU(x_sym + .15), YU(ys), just = c("left", "centre"), gp = gp(PT_KEY))
    ys <- ys + KEY_STEP
  }
  grid.text(">", XU(x_sym + .05), YU(ys), gp = gp(PT_GENE, "bold"))
  grid.text("size beyond 6 (drawn at 6)", XU(x_sym + .15), YU(ys), just = c("left", "centre"), gp = gp(PT_KEY))
}

# S10: G1 (A) above G2 (B), sharing one set of columns and keys. Positions in inches: x from the canvas's left edge,
# y down from its top; the canvas is 7.3 in wide and H inches tall, with its top-left corner at page pixel (120, 48).
s6_layout <- function() {
  t6 <- T6(); rng <- LV_RANGE()
  g1 <- t6[t6$block == "G1", ]; g2h <- t6[t6$block == "G2 higher in SC-SIC, Fig. 3D FDR < 0.05", ]; g2n <- t6[t6$block == "G2 no significant Fig. 3D difference", ]
  gl <- .04 + max(str_in(c(g1$gene, g2h$gene, g2n$gene), PT_GENE, "italic")); left <- gl + .05   # gene names right-aligned at gl
  p <- (PW - .05 - left - 2 * .05 - .22) / 14                                                      # column pitch
  R <- .105                                                                                       # row pitch
  y <- 0
  Y <- list(h1 = y); y <- y + HEADER_H
  Y$strip <- y + .15; y <- Y$strip + .04
  Y$b1 <- y; y <- y + nrow(g1) * R + .10
  Y$h2 <- y; y <- y + HEADER_H
  Y$sb_h <- y; y <- y + .145
  Y$b2h <- y; y <- y + nrow(g2h) * R + .06
  Y$sb_n <- y; y <- y + .145
  Y$b2n <- y; y <- y + nrow(g2n) * R
  Y$bot <- y; y <- y + label_depth(60) + .08
  Y$key <- y; y <- y + KEY_H
  list(g1 = g1, g2h = g2h, g2n = g2n, rng = rng, gl = gl, left = left, p = p, xs = col_x(left, p, .05, .22), R = R, Y = Y,
       H = ceiling((y + .04) * 600) / 600)                                                        # whole pixels at 600 dpi
}

# ---- Terms of the gene groups alone and with their mixed genes (S13) ---------------------------------------------------
PAIR <- c(G1 = "G1_PLUS_MIXED_G1_LIKE", G2 = "G2_PLUS_MIXED_G2_LIKE", G3 = "G3_PLUS_MIXED_G3_IF_NOT_SHARED")   # each group + its mixed genes
BLOCKS <- c("G3", "G2", "G1")                           # Figure 3G order

# The terms one panel draws, in drawn order (s8_terms.tsv: for each group, the terms Figure 3G draws for it first, then
# the best other terms called with or without the mixed genes; the name each term is drawn with; whether Figure 3G draws
# it). From the enrichment table all (one row per gene list and term; key gives each row's term as s8_terms.tsv names
# it): the term's adjusted p for the group alone (NA: no group gene in the term) and with its mixed genes, and the mixed
# genes in the term.
s8_terms <- function(panel, all, key) {
  LT <- read_table("s8_mixed/pathways/lists.tsv")
  D <- read.delim(file.path(dirname(PLOT_SCRIPT), "s8_terms.tsv"))
  D <- D[D$panel == panel, ]
  ia <- match(paste(D$block, D$term), paste(all$list, key))
  ib <- match(paste(PAIR[D$block], D$term), paste(all$list, key))
  data.frame(block = D$block, term = D$name, p_without = all$p_adj[ia], p_with = all$p_adj[ib], in3g = D$in_3g,
             mixed = vapply(seq_along(ib), function(i) paste(sort(setdiff(strsplit(all$genes[ib[i]], "/", fixed = TRUE)[[1]],
                                                                          LT$gene[LT$list == D$block[i]])), collapse = ", "), ""))
}

# One panel: per group a header (group size and the mixed genes it gains), then per term a line from the group alone
# (open dot; circled cross at 0 when the group has no gene in the term) to the group with its mixed genes (filled dot) on
# -log10 adjusted p, the term name left of x = 0 (bold: also in Figure 3G) and the mixed genes in the term right of the
# plot. Text sits at fixed inch offsets from the plot, so A and B line up. gtable_h: the height used to measure the plot's
# fixed parts. Returns the plot and the panel height h for plot_panel().
term_panel <- function(D, gtable_h) {
  LT <- read_table("s8_mixed/pathways/lists.tsv")
  NG <- vapply(names(PAIR), function(g) sum(LT$list == g), 0L)
  NA_ <- vapply(names(PAIR), function(g) length(setdiff(LT$gene[LT$list == PAIR[[g]]], LT$gene[LT$list == g])), 0L)
  BL <- c(G3 = sprintf("G3 (%d genes) vs G3 + %d mixed genes higher in SC-SIC (Fig. 3D; %d of the %d G2-like)", NG[["G3"]], NA_[["G3"]], NA_[["G3"]], NA_[["G2"]]),
          G2 = sprintf("G2 (%d genes) vs G2 + %d G2-like mixed genes", NG[["G2"]], NA_[["G2"]]),
          G1 = sprintf("G1 (%d genes) vs G1 + %d G1-like mixed genes", NG[["G1"]], NA_[["G1"]]))
  NOTE_BOLD <- 'bold("Bold")~"names: also shown in Fig. 3G"'   # note over the term names (plotmath)
  MM <- PT / .pt
  TERM_W <- 2.32; GENE_W <- 2.80                        # left of x = 0 / right of x = XMAX
  LAB_GAP <- 0.12; GENE_GAP <- 0.14; CANVAS_PAD <- 0.04
  P <- 0.150; LINE2 <- 0.118; HDR <- 0.200; BGAP <- 0.05; COLHDR <- 0.17; BOT <- 0.05   # row pitch, extra line, header, gaps
  XMAX <- 10; BREAKS <- seq(0, XMAX, 2); THR <- -log10(0.05)
  CHIP <- 0.085; HDR_X <- CANVAS_PAD + CHIP + 0.05
  KINDS <- c("Group alone", "Group with its mixed genes", "No group gene in the term")
  meas <- function(s, face = "plain")
    vapply(s, function(z) convertWidth(grobWidth(textGrob(z, gp = gpar(fontfamily = FONT, fontsize = PT, fontface = face))), "in", valueOnly = TRUE), 0)
  wrap_genes <- function(s, width) {                    # greedy wrap after a comma, by measured italic width
    if (!nzchar(s)) return("no mixed gene in the term")
    z <- strsplit(s, ", ", fixed = TRUE)[[1]]; lines <- character(); cur <- z[1]
    for (k in seq_along(z)[-1]) { try <- paste0(cur, ", ", z[k])
      if (meas(if (k == length(z)) try else paste0(try, ","), "italic") <= width) cur <- try else { lines <- c(lines, paste0(cur, ",")); cur <- z[k] } }
    paste(c(lines, cur), collapse = "\n") }

  D$mixed_lab <- vapply(D$mixed, wrap_genes, "", width = GENE_W - GENE_GAP - 0.06)
  D$face <- ifelse(nzchar(D$mixed), "italic", "plain")
  D$nl <- lengths(strsplit(D$mixed_lab, "\n", fixed = TRUE))
  D$tface <- ifelse(D$in3g, "bold", "plain")
  D$block <- factor(D$block, BLOCKS)
  # rows top to bottom (y in inches below the plot top): column header line, then per block a header and its rows
  y <- -COLHDR; HB <- NULL; D$yc <- NA_real_
  for (b in BLOCKS) {
    if (b != BLOCKS[1]) y <- y - BGAP
    top <- y - HDR; HB <- rbind(HB, data.frame(block = b, hy = y - HDR / 2, top = top, bottom = NA_real_))
    y <- top
    for (i in which(D$block == b)) { h <- P + (D$nl[i] - 1) * LINE2; D$yc[i] <- y - h / 2; y <- y - h }
    HB$bottom[HB$block == b] <- y
  }
  YMIN <- y - BOT
  D$v0 <- ifelse(is.na(D$p_without), 0, -log10(D$p_without)); D$v1 <- -log10(D$p_with)
  PTS <- rbind(data.frame(block = D$block, y = D$yc, x = D$v0, kind = ifelse(is.na(D$p_without), KINDS[3], KINDS[1])),
               data.frame(block = D$block, y = D$yc, x = D$v1, kind = KINDS[2]))
  PTS$kind <- factor(PTS$kind, KINDS[KINDS %in% PTS$kind])
  build <- function(upi, left_in) {                     # upi: x data units per inch; text placed at fixed inch offsets
    xl <- function(inch) inch * upi
    x_canvas <- -left_in * upi
    GR <- merge(HB, data.frame(x = BREAKS))
    ggplot() +
      geom_segment(data = GR, aes(x = x, xend = x, y = top, yend = bottom), colour = GRID_COL, linewidth = .2) +
      geom_segment(data = HB, aes(x = THR, xend = THR, y = top, yend = bottom, linetype = "Adjusted p = 0.05"), colour = "grey45", linewidth = .3,
                   key_glyph = "vline") +
      geom_segment(data = D, aes(x = v0, xend = v1, y = yc, yend = yc), colour = "grey72", linewidth = .35) +
      geom_point(data = PTS, aes(x = x, y = y, colour = block, shape = kind), fill = "white", size = 1.8, stroke = .7) +
      geom_text(data = D, aes(x = -xl(LAB_GAP), y = yc, label = term, fontface = tface), hjust = 1, size = MM, family = FONT, colour = "black") +
      geom_text(data = D, aes(x = XMAX + xl(GENE_GAP), y = yc, label = mixed_lab, fontface = face), hjust = 0, vjust = .5, lineheight = .95,
                size = MM, family = FONT, colour = "black") +
      annotate("rect", xmin = x_canvas + xl(CANVAS_PAD), xmax = x_canvas + xl(CANVAS_PAD + CHIP), ymin = HB$hy - CHIP / 2, ymax = HB$hy + CHIP / 2,
               fill = GROUP_COL[HB$block], colour = NA) +
      annotate("text", x = x_canvas + xl(HDR_X), y = HB$hy, label = BL[HB$block], hjust = 0, size = MM, family = FONT, fontface = "bold", colour = "black") +
      annotate("text", x = XMAX + xl(GENE_GAP), y = -COLHDR / 2, label = "Mixed genes in the term", hjust = 0, size = MM, family = FONT, colour = "black") +
      annotate("text", x = -xl(LAB_GAP), y = -COLHDR / 2, label = NOTE_BOLD, parse = TRUE, hjust = 1, size = MM, family = FONT, colour = "black") +
      annotate("segment", x = 0, xend = XMAX, y = YMIN, yend = YMIN, colour = "black", linewidth = .35) +
      scale_colour_manual(values = GROUP_COL, guide = "none") +
      scale_shape_manual(values = c("Group alone" = 21, "Group with its mixed genes" = 16, "No group gene in the term" = 13), name = NULL) +
      scale_linetype_manual(values = c("Adjusted p = 0.05" = "22"), labels = expression("Adjusted"~italic(p)~"= 0.05"), name = NULL) +
      scale_x_continuous(breaks = BREAKS, labels = scales::label_number(accuracy = 1, style_negative = "minus")) +
      coord_cartesian(xlim = c(0, XMAX), ylim = c(YMIN, 0), expand = FALSE, clip = "off") +
      labs(x = expression(-log[10]~"adjusted"~italic(p)), y = NULL) +
      panel_theme("none", legend = "bottom") +
      guides(shape = guide_legend(order = 1, override.aes = list(colour = "grey30", size = 1.8, stroke = .7)),
             linetype = guide_legend(order = 2, override.aes = list(colour = "grey45", linewidth = .3),
                                     theme = theme(legend.text = element_text(margin = margin(t = 1.5, l = 5.5, unit = "pt"))))) +
      theme(axis.line.x = element_blank(), axis.text.y = element_blank(), axis.ticks.y = element_blank(),
            legend.key.width = unit(.16, "in"), legend.key.height = unit(.13, "in"), legend.spacing.x = unit(.04, "in"),
            legend.margin = margin(0, 0, 0, 0), legend.box.margin = margin(2, 0, 0, 0), legend.box.spacing = unit(.03, "in"),
            legend.justification = "center",
            plot.margin = unit(c(0, GENE_W, .02, TERM_W), "in"))
  }
  geom_in <- function(g) {                              # fixed widths/heights of the plot around its one flexible panel cell
    gt <- ggplotGrob(g); pushViewport(viewport(width = unit(PW, "in"), height = unit(gtable_h, "in")))
    wd <- convertWidth(gt$widths, "in", valueOnly = TRUE); ht <- convertHeight(gt$heights, "in", valueOnly = TRUE); popViewport()
    nw <- unitType(gt$widths) == "null"; nh <- unitType(gt$heights) == "null"; pl <- gt$layout[gt$layout$name == "panel", ]
    list(left = sum(wd[seq_len(pl$l - 1)]), fixed_w = sum(wd[!nw]), fixed_h = sum(ht[!nh])) }
  G0 <- geom_in(build(XMAX / 2, TERM_W))                # first pass measures the fixed parts, second places the text in inches
  g <- build(XMAX / (PW - G0$fixed_w), G0$left)
  list(plot = g, h = .47 + geom_in(g)$fixed_h - YMIN)
}

# ---- Groups of S12 -----------------------------------------------------------------------------------------------------
SIC <- "#AA3377"; BETA <- "#4477AA"; ALPHA <- "#117733"
SHP <- c("Whole cells" = 16, "Single nuclei" = 1)
NUC_STUDY <- c("GSE199636", "SF_HICKEY2023_ADULT_EEC", "Siletti 2023")   # single-nucleus studies: open symbols
NAMES <- c("Brain serotonin neurons" = "Brain serotonin neurons, adult", "Brain dopamine neurons" = "Brain dopamine neurons, adult",
           "Primary gut EC" = "Native intestinal EC", "HIO EC" = "HIO graft EC", "Induced EC" = "Induced EC", "SC-SIC" = "SC-SIC", "SC-beta" = "SC-beta",
           "SC-alpha" = "SC-alpha", "Adult islet beta" = "Adult islet beta")
GCOL <- setNames(c("#222222", "#888888", "#43A5B6", "#6943AC", "#CCBB44", SIC, BETA, ALPHA, "#1C3A63"), NAMES)
