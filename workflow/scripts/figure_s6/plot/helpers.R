# Drawing helpers of the Figure S6 panels (the same sizes, colors, theme and key shapes as Figure 3). Sourced by
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
