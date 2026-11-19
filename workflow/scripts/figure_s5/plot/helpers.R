# Drawing helpers shared by the Figure S5 panels: type sizes, the panel theme, placement on the page and key rows.
# Sourced by workflow/scripts/common/draw_panel.R, which defines the page size W x H (inches) and DPI.
# Every position is in inches from the page's top-left corner.

FONT <- "Arial"
PT <- 7.5                         # body text (pt): axis text, labels, keys
MM <- PT / .pt                    # the same size in ggplot's text units (mm)

# The theme every panel starts from.
theme_s5 <- function() theme_minimal(base_family = FONT, base_size = PT) + theme(
  text = element_text(family = FONT, colour = "black", size = PT),
  axis.text = element_text(colour = "black", size = PT), axis.title = element_text(colour = "black", size = PT),
  strip.text = element_text(colour = "black", size = PT, face = "bold", hjust = 0, margin = margin(1, 0, 2, 0)),
  panel.grid.minor = element_blank(), panel.grid.major = element_line(colour = "#E4E4E4", linewidth = 0.2),
  axis.line.x = element_line(colour = "black", linewidth = 0.3), axis.ticks.x = element_line(colour = "black", linewidth = 0.3),
  axis.ticks.length = unit(1.5, "pt"), legend.position = "none", plot.margin = margin(1, 3, 1, 1))

# A drawing-only table shipped next to the plot scripts.
read_placement <- function(name) read.delim(file.path(dirname(PLOT_SCRIPT), name), stringsAsFactors = FALSE)

# Width of a string in inches on the current device.
str_w <- function(s, size = PT, face = "plain") {
  pushViewport(viewport(gp = gpar(fontfamily = FONT, fontsize = size, fontface = face)))
  w <- convertWidth(stringWidth(s), "in", valueOnly = TRUE); popViewport(); w
}

# Text with its top-left corner at (x, y).
text_at <- function(s, x, y, size = PT, face = "plain")
  grid.text(s, x = unit(x, "in"), y = unit(H - y, "in"), just = c("left", "top"),
            gp = gpar(fontfamily = FONT, fontsize = size, fontface = face, col = "black"))

# A ggplot with its top-left corner at (x, y), drawn w x h inches.
draw_at <- function(p, x, y, w, h) {
  pushViewport(viewport(x = x / W, y = 1 - y / H, width = w / W, height = h / H, just = c("left", "top")))
  grid.draw(ggplotGrob(p)); popViewport()
}

# A row of key items centered in [left, left + width], row top at y.
key_row <- function(items, left, width, y, gap = .13) {
  widths <- vapply(items, function(it) it$w_sym + str_w(it$label), 0)
  x <- left + (width - (sum(widths) + gap * (length(items) - 1))) / 2
  for (i in seq_along(items)) { items[[i]]$draw(x, y); text_at(items[[i]]$label, x + items[[i]]$w_sym, y - .005); x <- x + widths[i] + gap }
}
