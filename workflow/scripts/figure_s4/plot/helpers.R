# Drawing helpers of the Figure S4 panels (the same sizes, colors, theme and key shapes as Figure 3). Sourced by
# workflow/scripts/common/draw_panel.R, which defines the page size W x H (inches) and DPI.
suppressPackageStartupMessages({ library(grid); library(ggplot2); library(data.table) })

FONT <- "Arial"; PT <- 7.5; PT_GENE <- 7
GRID_COL <- "#E4E4E4"
GROUP_COL <- c(G1 = "#43A5B6", G2 = "#A89A2E", G3 = "#B4553F")                      # gene groups
POP_COL <- c("Primary intestinal EC" = "#43A5B6", "Other primary EEC" = "#8C8C8C", "HIO graft EC" = "#6943AC",
             "SC-SIC" = "#AA3377", "SC-beta" = "#4477AA", "SC-alpha" = "#117733", "SC-delta" = "#AE8C44",
             "Adult islet beta" = "#1C3A63", "Adult islet alpha" = "#0B4F26")

# A panel drawn on a canvas of its own, w x h inches (whole pixels), then placed with its top-left corner at page
# pixel (x, y). draw() uses inches from the canvas's top-left corner (XU, YU below).
panel_canvas <- function(x, y, w, h, draw) {
  w_px <- round(w * DPI); h_px <- round(h * DPI)
  f <- tempfile(fileext = ".png")
  png(f, width = w_px, height = h_px, units = "px", res = DPI, type = "cairo", bg = "white")
  grid.newpage(); reset_metric_cache(); draw(); invisible(dev.off())
  grid.raster(png::readPNG(f, native = TRUE), x = unit(x / DPI, "in"), y = unit((H * DPI - y) / DPI, "in"),
              width = unit(w_px / DPI, "in"), height = unit(h_px / DPI, "in"), just = c("left", "top"), interpolate = FALSE)
  unlink(f)
}

# R caches the height of "M" per device; measuring at an unused size clears it, so text lands the same on every device
reset_metric_cache <- function() {
  pushViewport(viewport(gp = gpar(fontfamily = FONT, fontsize = 1.234))); convertHeight(stringHeight("M"), "in"); popViewport()
}

panel_theme <- function() theme_minimal(base_family = FONT, base_size = PT) + theme(
  text = element_text(colour = "black", size = PT), axis.text = element_text(colour = "black", size = PT),
  axis.title = element_text(colour = "black", size = PT), strip.text = element_text(colour = "black", size = PT, face = "bold", hjust = 0),
  legend.position = "none", panel.grid.minor = element_blank(), panel.grid.major.y = element_blank(),
  panel.grid.major.x = element_line(colour = GRID_COL, linewidth = .2),
  axis.line.x = element_line(colour = "black", linewidth = .35), axis.ticks.x = element_line(colour = "black", linewidth = .3),
  axis.ticks.length = unit(1.3, "pt"), plot.margin = margin(2, 6, 2, 2))

# grid text and shapes with inches measured from the top-left of the current viewport (height .VH)
.VH <- NA_real_
YU <- function(y) unit(.VH - y, "in"); XU <- function(x) unit(x, "in")
gp <- function(pt = PT, face = "plain", ...) gpar(fontfamily = FONT, fontsize = pt, fontface = face, col = "black", ...)
str_in <- function(s, pt = PT, face = "plain") {
  pushViewport(viewport(gp = gpar(fontfamily = FONT, fontsize = pt, fontface = face))); on.exit(popViewport())
  vapply(s, function(z) convertWidth(stringWidth(z), "in", valueOnly = TRUE), 0, USE.NAMES = FALSE)
}
with_measure <- function(expr) { f <- tempfile(fileext = ".png"); png(f, 7.5, 10, units = "in", res = 600, type = "cairo")
  on.exit({ dev.off(); unlink(f) }); force(expr) }

# a color ramp 0.085 in tall
ramp_bar <- function(cols, x, y, w)
  grid.raster(matrix(cols, nrow = 1), XU(x), YU(y), width = unit(w, "in"), height = unit(.085, "in"), just = c("left", "top"), interpolate = TRUE)

# key items drawn as shapes (symbol from x, row top y, label after w_sym) and a centered row of them
key_item <- function(label, w_sym, draw) list(label = label, w_sym = w_sym, draw = draw)
key_row <- function(items, left, width, y, gap = .16) {
  widths <- vapply(items, function(it) it$w_sym + str_in(it$label), 0)
  x <- left + (width - (sum(widths) + gap * (length(items) - 1))) / 2
  if (x < left - 1e-6) stop("key row wider than its box: ", round(sum(widths) + gap * (length(items) - 1), 3), " > ", width)
  for (i in seq_along(items)) { items[[i]]$draw(x, y)
    grid.text(items[[i]]$label, XU(x + items[[i]]$w_sym), YU(y + .0425), just = c("left", "centre"), gp = gp())
    x <- x + widths[i] + gap }
}
sq_item <- function(label, fill, col = NA, lwd = .5) key_item(label, .14, function(x, y)
  grid.rect(XU(x), YU(y), width = unit(.085, "in"), height = unit(.085, "in"), just = c("left", "top"), gp = gpar(fill = fill, col = col, lwd = lwd)))
dot_item <- function(label, open = FALSE) key_item(label, .14, function(x, y)
  grid.points(XU(x + .0425), YU(y + .0425), pch = if (open) 21 else 16, size = unit(.08, "in"), gp = gpar(col = "grey25", fill = "white", lwd = .9)))
tick_item <- function(label) key_item(label, .11, function(x, y)
  grid.lines(XU(c(x + .04, x + .04)), YU(c(y - .005, y + .09)), gp = gpar(col = "grey15", lwd = 1.6, lineend = "butt")))

# a panel letter, 11 pt bold, top-left at (x, y) page inches
letter_at <- function(l, x, y) grid.text(l, x = unit(x, "in"), y = unit(H - y, "in"), just = c("left", "top"), gp = gp(11, "bold"))

# the gene groups as named in the figures
GLAB <- c(G1 = "G1: EC genes low in SC-SICs (28)", G2 = "G2: shared EC genes (34)", G3 = "G3: separate SC-SIC program (84)")

# ---- Figure S4 page: both panels share one page layout (inches from the page's top-left corner) ----
COLS <- data.table(study = c("Elmentaite 2020 + 2021", "Egozi 2021", "Hickey 2023", "McCauley 2023", "Augsornworawat 2020",
                             "Balboa 2022", "Du 2022", "Augsornworawat 2023", "Zhu 2023", "Siehler 2024"),
                   pop = c(rep("Primary intestinal EC", 3), "HIO graft EC", rep("SC-SIC", 6)))
COLS[, label := study][study %in% c("Hickey 2023", "Augsornworawat 2023"), label := paste(study, "(nuclei)")]
COLS[, grp := as.integer(factor(pop, unique(pop)))]
L <- as.data.table(read_table("g2_g3_pooled_levels.tsv")); O <- as.data.table(read_table("gene_order.tsv"))
# the color scale of Figure 3B: magma, stretched toward the low end, from 0 to the ceiling of the 3B values
B3 <- read_table("g1_pooled_levels.tsv", figure = "figure3"); LV_RANGE <- c(0, ceiling(max(as.matrix(B3[, COLS$study]))))
MAGMA <- viridisLite::magma(256, begin = .08, end = .96); MAGMA_AT <- seq(0, 1, length.out = 256)^1.5
lv_col <- function(v) scales::gradient_n_pal(MAGMA, MAGMA_AT, "Lab")(scales::rescale(v, from = LV_RANGE))
G2A <- O[group == "G2" & part == "higher in SC-SICs"][order(position)]$gene
G2B <- O[group == "G2" & part == "no significant difference"][order(position)]$gene
G3O <- O[group == "G3"][order(position)]$gene
stopifnot(length(G2A) == 9L, length(G2B) == 25L, length(G3O) == 84L)
lev <- function(genes) as.matrix(L[match(genes, gene), COLS$study, with = FALSE])

MARGIN <- .15; X0 <- .30; SGAP <- .15; SW <- (W - MARGIN - .02 - X0 - 2 * SGAP) / 3            # three slots per row
SX <- X0 + (0:2) * (SW + SGAP)
M <- with_measure(list(gene_w = max(str_in(L$gene, PT_GENE, "italic")), lab_w = max(str_in(COLS$label, PT_GENE)),
                       h1 = str_in("Higher in SC-SICs, FDR < 0.05 (9)", PT, "bold"), h2 = str_in("No significant difference (25)", PT, "bold")))
GG <- .04                                                          # gap between column groups
CW <- (SW - M$gene_w - .05 - 2 * GG) / nrow(COLS)
R <- .123                                                          # row pitch
if (CW < .14) stop("columns narrower than 0.14 in")
colx <- function(k) SX[k] + M$gene_w + .05 + (seq_len(nrow(COLS)) - .5) * CW + (COLS$grp - 1) * GG   # column centers in slot k
HDR_H <- .18; SECT_H <- .17
A_HDR <- .20; A_LAB0 <- A_HDR + HDR_H                             # group header, then rotated column labels
A_BAR <- A_LAB0 + M$lab_w + .04; A_DATA <- A_BAR + .06 + .04
A_ROW1 <- A_DATA + SECT_H                                          # first data row
A_BOT <- A_ROW1 + 9 * R + .06 + SECT_H + 8 * R
B_HDR <- A_BOT + .22; B_BAR <- B_HDR + HDR_H; B_ROW1 <- B_BAR + .06 + .04
B_ROWS_END <- B_ROW1 + 28 * R

.VH <- H
cells <- function(k, y0, genes) {                                  # one block: tiles, group outlines, gene names
  x <- colx(k); v <- lev(genes); n <- length(genes)
  for (j in seq_len(nrow(COLS))) grid.rect(XU(x[j]), YU(y0 + (seq_len(n) - .5) * R), width = unit(CW, "in"), height = unit(R, "in"),
                                           gp = gpar(fill = lv_col(v[, j]), col = "white", lwd = .2 * 72.27 / 25.4))
  for (q in unique(COLS$grp)) { i <- which(COLS$grp == q)
    grid.rect(XU(min(x[i]) - CW / 2), YU(y0), width = unit(length(i) * CW, "in"), height = unit(n * R, "in"), just = c("left", "top"),
              gp = gpar(fill = NA, col = "grey75", lwd = .45)) }
  grid.text(genes, XU(SX[k] + M$gene_w), YU(y0 + (seq_len(n) - .5) * R), just = c("right", "centre"), gp = gp(PT_GENE, "italic"))
}
bars <- function(k, y) { x <- colx(k); for (q in unique(COLS$grp)) { i <- which(COLS$grp == q)
  grid.rect(XU(min(x[i]) - CW / 2), YU(y), width = unit(length(i) * CW, "in"), height = unit(.06, "in"), just = c("left", "top"),
            gp = gpar(fill = POP_COL[[COLS$pop[i[1]]]], col = NA)) } }
labels_up <- function(k, y_base) grid.text(COLS$label, XU(colx(k)), YU(y_base), rot = 90, just = c("left", "centre"), gp = gp(PT_GENE))
labels_down <- function(k, y_top) grid.text(COLS$label, XU(colx(k)), YU(y_top), rot = 90, just = c("right", "centre"), gp = gp(PT_GENE))
group_head <- function(lab, col, y) { grid.points(XU(X0 + .10), YU(y + .06), pch = 16, size = unit(.09, "in"), gp = gpar(col = col))
  grid.text(lab, XU(X0 + .18), YU(y + .06), just = c("left", "centre"), gp = gp(PT, "bold")) }
