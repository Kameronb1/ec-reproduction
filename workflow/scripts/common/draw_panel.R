# Draw one figure panel and save it as a 600 dpi PNG on a white background.
#
#   Rscript draw_panel.R <plot script> <box> <results dir> <out.png> <exact|tolerance>
#
# <plot script> is workflow/scripts/<figure>/<plot folder>/<panel>.R. It reads its tables with read_table() and
# defines draw_panel(), which draws the panel at its place on the figure page. helpers.R in the same
# folder is sourced first. <box> is "x,y,w,h": the panel's rectangle on the page, in pixels from the
# page's top-left corner.
#
# The panel is drawn alone on the whole figure page (7.5 x 10 in at 600 dpi, cairo device, Arial), then
# cut out by its box. Drawing on the whole page keeps every line and letter on the same pixel grid as in
# the published figure. The page outside the box must stay white, so the cut loses nothing. Anything
# drawn outside the box stops the run with "exact"; with "tolerance" it is only a warning, because fonts
# that measure differently can push text a few pixels out.

suppressPackageStartupMessages({ library(ggplot2); library(patchwork); library(grid) })
if (.Platform$OS.type == "windows") stopifnot(nzchar(Sys.setlocale("LC_CTYPE", "English_United States.utf8")))

args <- commandArgs(TRUE)
stopifnot(length(args) == 5, args[5] %in% c("exact", "tolerance"))
PLOT_SCRIPT <- normalizePath(args[1], winslash = "/", mustWork = TRUE)
BOX <- as.integer(strsplit(args[2], ",", fixed = TRUE)[[1]])
RESULTS <- normalizePath(args[3], winslash = "/", mustWork = TRUE)
OUT <- args[4]
FIGURE <- basename(dirname(dirname(PLOT_SCRIPT)))
REPO <- normalizePath(file.path(dirname(PLOT_SCRIPT), "..", "..", "..", ".."), winslash = "/", mustWork = TRUE)

W <- 7.5; H <- 10; DPI <- 600         # figure page, inches
stopifnot(length(BOX) == 4, BOX[1] >= 0, BOX[2] >= 0, BOX[1] + BOX[3] <= W * DPI, BOX[2] + BOX[4] <= H * DPI)

# A table made by the analysis rules: results/<figure>/<name>.
read_table <- function(name, figure = FIGURE)
  read.delim(file.path(RESULTS, figure, name), check.names = FALSE, fileEncoding = "UTF-8")
result_path <- function(name, figure = FIGURE) file.path(RESULTS, figure, name)

# The plots are built while the device is open, because some of them measure their text.
page <- tempfile(fileext = ".png")
png(page, width = W, height = H, units = "in", res = DPI, type = "cairo", bg = "white")
source(file.path(dirname(PLOT_SCRIPT), "helpers.R"))
source(PLOT_SCRIPT)
grid.newpage()
draw_panel()
invisible(dev.off())

# One packed RGBA integer per pixel, stored row by row; as a matrix, img[y, x].
raw <- png::readPNG(page, native = TRUE)
img <- matrix(as.integer(raw), nrow = H * DPI, ncol = W * DPI, byrow = TRUE)
white <- -1L                                      # opaque white
rows <- BOX[2] + seq_len(BOX[4]); cols <- BOX[1] + seq_len(BOX[3])
outside <- img; outside[rows, cols] <- white
if (any(outside != white)) {
  msg <- sprintf("%s draws %d pixels outside its box %s", basename(PLOT_SCRIPT), sum(outside != white), args[2])
  if (args[5] == "exact") stop(msg, call. = FALSE) else warning(msg, call. = FALSE)
}
cut <- structure(as.vector(t(img[rows, cols])), dim = c(BOX[4], BOX[3]), class = "nativeRaster",
                 channels = attr(raw, "channels"))
dir.create(dirname(OUT), recursive = TRUE, showWarnings = FALSE)
png::writePNG(cut, OUT, dpi = DPI)
unlink(page)
cat("wrote", OUT, "\n")
