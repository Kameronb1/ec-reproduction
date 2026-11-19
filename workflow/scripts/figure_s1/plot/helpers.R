# Drawing helpers shared by the Supplementary Figure S1 panels: fonts, colors, the plot theme, text and key helpers, the
# page layout and the preliminary 72 dpi draw. Sourced by workflow/scripts/common/draw_panel.R, which defines the page
# size W x H (inches). Positions are in inches from the page's top-left corner.
suppressPackageStartupMessages({ library(patchwork); library(ggrepel) })

FONT <- "Arial"
PT <- 7; LETTER <- 11; MM <- PT / .pt
EDGE <- .15                                                    # hairline black outline of points
COL <- c(gut_ec = "#43A5B6", sc_sic = "#AA3377", beta = "#4477AA", alpha = "#117733", delta = "#AE8C44", hio_ec = "#6943AC")
SC_ORDER <- c("Augsornworawat 2020", "Balboa 2022", "Du 2022", "Augsornworawat 2023", "Zhu 2023", "Siehler 2024")
GSE <- c(GSE151117 = "Augsornworawat 2020", GSE167880 = "Balboa 2022", GSE185036 = "Du 2022", GSE199636 = "Augsornworawat 2023",
         GSE202497 = "Zhu 2023", GSE226346 = "Siehler 2024")
NUCLEI <- c("Hickey 2023", "Augsornworawat 2023")             # single-nucleus studies (open symbols)
POPCOL <- c(GUT_EC = COL[["gut_ec"]], SC_SIC = COL[["sc_sic"]], SC_BETA = COL[["beta"]], SC_ALPHA = COL[["alpha"]], HIO_GRAFT_EC = COL[["hio_ec"]])
fmt <- function(n) format(n, big.mark = ",", trim = TRUE)
minus <- scales::label_number(style_negative = "minus")
theme_pg <- function(grid = "x") theme_minimal(base_family = FONT, base_size = PT) + theme(
  text = element_text(family = FONT, colour = "black", size = PT),
  axis.text = element_text(colour = "black", size = PT), axis.title = element_text(colour = "black", size = PT),
  strip.text = element_text(colour = "black", size = PT, face = "bold", hjust = 0, vjust = 0, lineheight = 1, margin = margin(1, 0, 3, 0)),
  panel.grid.minor = element_blank(),
  panel.grid.major.x = if (grid %in% c("x", "both")) element_line(colour = "#E4E4E4", linewidth = .2) else element_blank(),
  panel.grid.major.y = if (grid %in% c("y", "both")) element_line(colour = "#E4E4E4", linewidth = .2) else element_blank(),
  axis.line.x = element_line(colour = "black", linewidth = .3), axis.ticks.x = element_line(colour = "black", linewidth = .3),
  axis.ticks.length = unit(1.5, "pt"), legend.position = "none", plot.margin = margin(1, 3, 1, 1))

# ---- page layout (in): two panels side by side (A, B), then C and D full width, each with its key below -------------
LX <- .10; X0 <- .30; XR <- W - .15                       # letter column, panel left, right edge
A_W <- 4.20; B_X <- X0 + A_W + .30; B_W <- XR - B_X
A_H <- 2.55; B_H <- 2.85; C_H <- 2.30; D_H <- 2.30
ROW1 <- .08
KEY_A_Y <- ROW1 + A_H + .10
KEY_B_Y <- ROW1 + B_H + .04
KEY_B_LINE <- .12                                          # the triangle entry's second line
ROW2 <- max(KEY_A_Y + 2 * .15, KEY_B_Y + 4 * .15 + KEY_B_LINE) + .20
KEY_C_Y <- ROW2 + C_H + .04
ROW3 <- KEY_C_Y + 2 * .15 + .20
KEY_D_Y <- ROW3 + D_H + .04

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

# ---- grid helpers (inches from the page's top-left) -----------------------------------------------------------------
tx <- function(s, x, y, size = PT, face = "plain", just = c("left", "top"), col = "black")
  grid.text(s, unit(x, "in"), unit(H - y, "in"), just = just, gp = gpar(fontfamily = FONT, fontsize = size, fontface = face, col = col, lineheight = 1))
str_in <- function(s, size = PT, face = "plain") { pushViewport(viewport(gp = gpar(fontfamily = FONT, fontsize = size, fontface = face)))
  on.exit(popViewport()); max(vapply(strsplit(s, "\n")[[1]], function(z) convertWidth(stringWidth(z), "in", valueOnly = TRUE), 0)) }
letter <- function(L, x, y) tx(L, x, y, LETTER, "bold")
place <- function(p, x, y, w, h) { pushViewport(viewport(unit(x, "in"), unit(H - y, "in"), unit(w, "in"), unit(h, "in"), just = c("left", "top")))
  grid.draw(if (inherits(p, "patchwork")) patchworkGrob(p) else ggplotGrob(p)); popViewport() }
# key symbols, centred at (x, y)
sym_dot <- function(x, y, fill, size = .075, open = FALSE) {
  if (open) { grid.points(unit(x, "in"), unit(H - y, "in"), pch = 1, size = unit(size + .02, "in"), gp = gpar(col = "black", lwd = EDGE * .stroke / 2))
    grid.points(unit(x, "in"), unit(H - y, "in"), pch = 1, size = unit(size - .01, "in"), gp = gpar(col = fill, lwd = .5 * .stroke / 2)) }
  else grid.points(unit(x, "in"), unit(H - y, "in"), pch = 21, size = unit(size, "in"), gp = gpar(col = "black", fill = fill, lwd = EDGE * .stroke / 2)) }
sym_tri <- function(x, y, fill) grid.points(unit(x, "in"), unit(H - y, "in"), pch = 24, size = unit(.07, "in"), gp = gpar(col = "black", fill = fill, lwd = EDGE * .stroke / 2))
sym_line <- function(x, y, col, lty = "solid", lwd = .5, w = .2) grid.segments(unit(x - w / 2, "in"), unit(H - y, "in"), unit(x + w / 2, "in"), unit(H - y, "in"),
                                                                            gp = gpar(col = col, lty = lty, lwd = lwd * .pt, lineend = "butt"))
sym_box <- function(x, y, fill, w = .16, h = .1) grid.rect(unit(x, "in"), unit(H - y, "in"), unit(w, "in"), unit(h, "in"), gp = gpar(fill = fill, col = NA))
# one key row: items = list(list(draw = function(x, y), w_sym, label)); centred in [left, left + width], text top at y
key_row <- function(items, left, width, y, gap = .16) {
  wd <- vapply(items, function(it) it$w_sym + str_in(it$label), 0); tot <- sum(wd) + gap * (length(items) - 1)
  if (tot > width + 1e-6) stop("key row wider than its box: ", round(tot, 3), " > ", width)
  x <- left + (width - tot) / 2
  for (i in seq_along(items)) { items[[i]]$draw(x + items[[i]]$w_sym / 2 - .02, y + .05); tx(items[[i]]$label, x + items[[i]]$w_sym, y); x <- x + wd[i] + gap }
}
item <- function(label, draw, w_sym = .16) list(label = label, draw = draw, w_sym = w_sym)
