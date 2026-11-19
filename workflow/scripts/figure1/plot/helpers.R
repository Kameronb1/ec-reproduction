# Drawing helpers shared by the Figure 1 panels: fonts, colors, text and key helpers, the page rows and the preliminary
# 72 dpi draw. Sourced by workflow/scripts/common/draw_panel.R, which defines the page size W x H (inches).
# Each panel is drawn inside its row of the page: positions are in inches from the row's top-left corner, y grows
# downwards and PH is the row height.

PT <- 7.5; TAG <- 11                     # body text and panel letter sizes (pt), Arial throughout
SEC <- "grey35"                          # secondary text (strip labels, tick labels, key labels)
COL <- c(gut_ec = "#43A5B6", sc_sic = "#AA3377", beta = "#4477AA", alpha = "#117733", delta = "#AE8C44",
         other = "#C9D3DA", hio_ec = "#6943AC")

# Rows of the page (inches from the page top): the B/C row, D under it, E under D.
ROW_TOP_A <- .045
ROW_H <- c(A = 1.30, BC = 3.84, D = 2.068 - .104, E = 1.50 + .26)
ROW_TOP <- c(A = ROW_TOP_A, BC = ROW_TOP_A + ROW_H[["A"]] + .14)
ROW_TOP[["D"]] <- ROW_TOP[["BC"]] + ROW_H[["BC"]]
ROW_TOP[["E"]] <- ROW_TOP[["D"]] + ROW_H[["D"]] + .06

# Draw draw(PH) inside a row: a viewport the page's width wide, its top ROW_TOP[[row]] inches below the page top.
in_row <- function(row, draw) {
  pushViewport(viewport(x = 0, y = unit(H - ROW_TOP[[row]], "in"), width = unit(W, "in"), height = unit(ROW_H[[row]], "in"),
                        just = c("left", "top")))
  draw(ROW_H[[row]])
  popViewport()
}

# ggplot2 measures each font's descent once per device type and reuses it for the rest of the session. The published
# page is drawn at 72 dpi before 600 dpi, so its plots' text heights are the 72 dpi measurements: draw() runs once on a
# 72 dpi page first, then on the open page.
draw_after_72dpi <- function(draw) {
  page <- dev.cur(); f <- tempfile(fileext = ".png")
  png(f, width = W, height = H, units = "in", res = 72, type = "cairo")
  grid.newpage(); draw()
  dev.off(); dev.set(page); unlink(f)
  draw()
}

# text at (x, y) in the row (inches from the row's top-left)
tx <- function(s, x, y, PH, size = PT, bold = FALSE, just = c("left", "top"), col = "black", italic = FALSE)
  grid.text(s, unit(x, "in"), unit(PH - y, "in"), just = just,
            gp = gpar(fontfamily = "Arial", fontsize = size, col = col, lineheight = 1,
                      fontface = if (bold && italic) "bold.italic" else if (bold) "bold" else if (italic) "italic" else "plain"))
str_in <- function(s, size = PT, bold = FALSE, italic = FALSE) {        # width of a string (in) on the current device
  pushViewport(viewport(gp = gpar(fontfamily = "Arial", fontsize = size,
                                  fontface = if (bold && italic) "bold.italic" else if (bold) "bold" else if (italic) "italic" else "plain")))
  z <- convertWidth(stringWidth(s), "in", valueOnly = TRUE); popViewport(); z
}
cap_in <- function(size) {                               # cap height (in) of bold Arial at this size, from the glyph "H"
  pushViewport(viewport(gp = gpar(fontfamily = "Arial", fontsize = size, fontface = "bold")))
  z <- convertHeight(grobHeight(textGrob("H")), "in", valueOnly = TRUE); popViewport(); z
}
# the panel letter: its top at y, drop inches lower
panel_heading <- function(letter, x, y, PH, drop = 0)
  tx(letter, x, y + cap_in(TAG) + drop, PH, TAG, TRUE, just = c("left", "bottom"))

# key items in row coordinates: a symbol, then its label
square_item <- function(label, col) list(label = label, w_sym = .14, draw = function(x, y, PH)
  grid.rect(x = unit(x, "in"), y = unit(PH - y - .04, "in"), width = unit(.085, "in"), height = unit(.085, "in"), just = "left", gp = gpar(fill = col, col = NA)))
dash_item <- function(label) list(label = label, w_sym = .14, draw = function(x, y, PH)
  grid.lines(x = unit(x + .045, "in"), y = unit(PH - y - c(-.01, .09), "in"), gp = gpar(col = "grey20", lwd = 0.45 * .pt, lty = "22")))
dot_key <- function(label, col) list(label = label, w_sym = .10, draw = function(x, y, PH)
  grid.points(unit(x + .03, "in"), unit(PH - y - .045, "in"), pch = 16, size = unit(.075, "in"), gp = gpar(col = col)))
text_item <- function(label) list(label = label, w_sym = 0, draw = function(x, y, PH) invisible())
# a row of key items centred in [left, left + width], row top at y; stops if wider than its box
key_row <- function(items, left, width, y, PH, gap = .13, size = PT) {
  widths <- vapply(items, function(it) it$w_sym + str_in(it$label, size), 0)
  total <- sum(widths) + gap * (length(items) - 1)
  if (total > width + 1e-6) stop("key row wider than its box: ", round(total, 3), " > ", width)
  x <- left + (width - total) / 2
  for (i in seq_along(items)) { items[[i]]$draw(x, y, PH); tx(items[[i]]$label, x + items[[i]]$w_sym, y - .005, PH, size); x <- x + widths[i] + gap }
}
