# Drawing helpers shared by the Figure S14 and S16 panels: placement on the page, text measuring, table reading and
# colors. Sourced by workflow/scripts/common/draw_panel.R, which defines the page size W x H (inches) and DPI.
#
# Each panel is laid out in inches on a strip of the page: w_px pixels wide from page column x_px, h_px pixels high from
# page row y_px (S14A and S16A span the page's width; S16B and C share the left column beside S16D). The strip is
# drawn on a device of its own (w_px x h_px, 600 dpi, cairo, white) and placed on the page at whole pixels, as these
# pages are assembled: drawn through the page's offset instead, positions round differently and angled text can land a
# pixel off.

FONT <- "Arial"
in_strip <- function(y_px, h_px, draw, x_px = 0, w_px = W * DPI) {
  page <- dev.cur(); f <- tempfile(fileext = ".png")
  png(f, width = w_px, height = h_px, res = DPI, type = "cairo", bg = "white", family = FONT)
  grid.newpage(); draw(); invisible(dev.off()); dev.set(page)
  grid.raster(png::readPNG(f, native = TRUE), x = unit(x_px / DPI, "in"), y = unit(H - y_px / DPI, "in"), width = unit(w_px / DPI, "in"),
              height = unit(h_px / DPI, "in"), just = c("left", "top"), interpolate = FALSE)
  unlink(f) }
# width (in) of each string on the open device
text_w <- function(s, size, face = "plain")
  unname(vapply(s, function(x) convertWidth(grobWidth(textGrob(x, gp = gpar(fontfamily = FONT, fontsize = size, fontface = face))), "in", valueOnly = TRUE), 0))
# a results table with text_cols read as text (numbers drawn exactly as stored, e.g. "2.3")
read_text_cols <- function(name, text_cols)
  read.delim(result_path(name), check.names = FALSE, fileEncoding = "UTF-8", colClasses = setNames(rep("character", length(text_cols)), text_cols))
# a drawing-only table shipped in this folder
read_drawing <- function(name) read.delim(file.path(REPO, "workflow", "scripts", FIGURE, "plot_supp", name), check.names = FALSE, fileEncoding = "UTF-8")

POPS <- c("SC_SIC", "SC_BETA", "SC_ALPHA"); POP_LABEL <- c(SC_SIC = "SC-SIC", SC_BETA = "SC-beta", SC_ALPHA = "SC-alpha")
POP_COL <- c(SC_SIC = "#AA3377", SC_BETA = "#4477AA", SC_ALPHA = "#117733")
CAP7 <- 0.716 * 7 / 72                                         # Arial cap height at 7 pt (in)
