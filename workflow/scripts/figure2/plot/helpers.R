# Drawing helpers shared by the Figure 2 panels: type sizes, the panel theme, the page layout, placement on the page
# and key items. Sourced by workflow/scripts/common/draw_panel.R, which defines the page size W x H (inches) and DPI.
# Every position is in inches from the page's top-left corner.

FONT <- "Arial"
PT <- 7.5                         # body text (pt): axis text, labels, keys
MM <- PT / .pt                    # the same size in ggplot's text units (mm)
LINE <- 1.2 * PT / 72             # one line of body text (in)

# The page layout: each plot's and each key's place (in).
LAY <- read.table(header = TRUE, text = "
item       x      y      w      h
A_plot     0.08   0.171  2.31   2.02
B_plot     2.50   0.13   4.95   2.12
C_plot     0.08   2.40   2.27   3.45
C_key      0.08   5.95   2.27   NA
D_ytitle   2.68   2.42   NA     NA
D_plot     2.42   2.70   2.36   3.15
D_key      2.42   5.95   2.36   NA
E_plot     4.84   2.56   2.61   3.45
E_key      4.84   6.11   2.61   NA
")
at <- function(item) { r <- LAY[LAY$item == item, ]; c(x = r$x, y = r$y, w = r$w, h = r$h) }

# The theme every panel starts from.
theme_fig2 <- function() theme_minimal(base_family = FONT, base_size = PT) + theme(
  text = element_text(family = FONT, colour = "black", size = PT),
  axis.text = element_text(colour = "black", size = PT), axis.title = element_text(colour = "black", size = PT),
  strip.text = element_text(colour = "black", size = PT, face = "bold", hjust = 0, margin = margin(1, 0, 2, 0)),
  panel.grid.minor = element_blank(), panel.grid.major = element_line(colour = "#E4E4E4", linewidth = 0.2),
  axis.line.x = element_line(colour = "black", linewidth = 0.3), axis.ticks.x = element_line(colour = "black", linewidth = 0.3),
  axis.ticks.length = unit(1.5, "pt"), legend.position = "none", plot.margin = margin(1, 3, 1, 1))

# Colors of the three gene groups and the mixed genes (Figures 2D and 2E).
GROUP_COLS <- c(G1 = "#43A5B6", G2 = "#A89A2E", G3 = "#B4553F", Mixed = "grey68")

# A drawing-only table shipped next to the plot scripts (label and dot placements).
read_placement <- function(name) read.delim(file.path(dirname(PLOT_SCRIPT), name), stringsAsFactors = FALSE)

# A plot's panel box (inches from the plot's top-left) and data ranges when the plot is drawn w x h inches, measured
# on a device like the page's.
panel_box <- function(p, w, h) {
  png(tempfile(fileext = ".png"), width = w, height = h, units = "in", res = DPI, type = "cairo"); on.exit(dev.off())
  g <- ggplotGrob(p)
  pushViewport(viewport(width = unit(w, "in"), height = unit(h, "in")))
  wd <- convertWidth(g$widths, "in", valueOnly = TRUE); ht <- convertHeight(g$heights, "in", valueOnly = TRUE); popViewport()
  nw <- unitType(g$widths) == "null"; nh <- unitType(g$heights) == "null"
  v <- as.numeric(g$widths[nw]); wd[nw] <- (w - sum(wd[!nw])) * v / sum(v)
  v <- as.numeric(g$heights[nh]); ht[nh] <- (h - sum(ht[!nh])) * v / sum(v)
  pn <- g$layout[g$layout$name == "panel", ]; cw <- c(0, cumsum(wd)); ch <- c(0, cumsum(ht))
  pp <- ggplot_build(p)$layout$panel_params[[1]]
  list(l = cw[pn$l], r = cw[pn$r + 1], t = ch[pn$t], b = ch[pn$b + 1], xr = pp$x$continuous_range, yr = pp$y$continuous_range)
}

# Width of a string in inches on the current device.
str_w <- function(s, size = PT, face = "plain") {
  pushViewport(viewport(gp = gpar(fontfamily = FONT, fontsize = size, fontface = face)))
  w <- convertWidth(stringWidth(s), "in", valueOnly = TRUE); popViewport(); w
}

# Text with its top-left corner at (x, y).
text_at <- function(s, x, y, size = PT, face = "plain")
  grid.text(s, x = unit(x, "in"), y = unit(H - y, "in"), just = c("left", "top"),
            gp = gpar(fontfamily = FONT, fontsize = size, fontface = face, col = "black"))

# A ggplot or patchwork with its top-left corner at (x, y), drawn w x h inches.
draw_at <- function(p, x, y, w, h) {
  pushViewport(viewport(x = x / W, y = 1 - y / H, width = w / W, height = h / H, just = c("left", "top")))
  grid.draw(if (inherits(p, "patchwork")) patchworkGrob(p) else ggplotGrob(p)); popViewport()
}
plot_at <- function(p, item) { b <- at(item); draw_at(p, b[["x"]], b[["y"]], b[["w"]], b[["h"]]) }

# Key items: a symbol drawn from x (row top y), its label w_sym to the right of x.
dot_item <- function(label, col) list(label = label, w_sym = .14, draw = function(x, y)
  grid.points(x = unit(x + .04, "in"), y = unit(H - y - .04, "in"), pch = 16, size = unit(.085, "in"), gp = gpar(col = col, fill = col, lwd = .8)))
ring_item <- function(label, col, fontsize, lwd) list(label = label, w_sym = .14, draw = function(x, y)   # an open ring, as ggplot draws shape 1
  grid.points(x = unit(x + .04, "in"), y = unit(H - y - .04, "in"), pch = 1, size = unit(1, "char"), gp = gpar(col = col, fontsize = fontsize, lwd = lwd)))
text_item <- function(label) list(label = label, w_sym = 0, draw = function(x, y) invisible())          # a label alone

# A row of key items centered in [left, left + width], row top at y.
key_row <- function(items, left, width, y, gap = .13) {
  widths <- vapply(items, function(it) it$w_sym + str_w(it$label), 0)
  x <- left + (width - (sum(widths) + gap * (length(items) - 1))) / 2
  for (i in seq_along(items)) { items[[i]]$draw(x, y); text_at(items[[i]]$label, x + items[[i]]$w_sym, y - .005); x <- x + widths[i] + gap }
}

# Key items in a grid, row by row (ncol columns, each as wide as its widest item), centered in [left, left + width];
# first row top at y.
key_grid <- function(items, left, width, y, ncol, gap = .25, step = .16) {
  widths <- vapply(items, function(it) it$w_sym + str_w(it$label), 0)
  col <- (seq_along(items) - 1) %% ncol + 1; row <- (seq_along(items) - 1) %/% ncol + 1
  cw <- vapply(seq_len(ncol), function(k) max(widths[col == k]), 0)
  total <- sum(cw) + gap * (ncol - 1); if (total > width) stop("key grid wider than its box: ", round(total, 3), " > ", width)
  x0 <- left + (width - total) / 2 + c(0, cumsum(cw + gap))[col]
  for (i in seq_along(items)) { yy <- y + (row[i] - 1) * step; items[[i]]$draw(x0[i], yy); text_at(items[[i]]$label, x0[i] + items[[i]]$w_sym, yy - .005) }
}
