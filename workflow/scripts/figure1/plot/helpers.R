# Drawing helpers shared by the Figure 1 panels: colors, the page rows, the preliminary 72 dpi draw and panel D's table. Sourced by
# workflow/scripts/common/draw_panel.R, which defines the page size W x H (inches).
# Each panel is drawn inside its row of the page: positions are in inches from the row's top-left corner, y grows
# downwards and PH is the row height.

SEC <- "grey35"                          # secondary text (strip labels, tick labels, key labels)
COL <- c(gut_ec = "#43A5B6", induced_ec = "#CCBB44", sc_sic = "#AA3377", beta = "#4477AA", alpha = "#117733",
         delta = "#AE8C44", other = "#C9D3DA", hio_ec = "#6943AC")

# Rows of the page (inches), on whole pixels at 600 dpi: the B/C row, D under it, and E filling the page below D
# down to 0.06 in above the page foot.
ROW_BC <- c(top = 1391 / 600, h = 2394 / 600)
ROW_D <- c(top = 3785 / 600, h = 1241 / 600)
ROW_E <- c(top = 5062 / 600, h = H - 5062 / 600 - .06)

# Draw draw(PH) inside a row: a viewport the page's width wide, its top `top` inches below the page top.
in_row <- function(row, draw) {
  pushViewport(viewport(x = 0, y = unit(H - row[["top"]], "in"), width = unit(W, "in"), height = unit(row[["h"]], "in"),
                        just = c("left", "top")))
  draw(row[["h"]])
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

# Panel D's per-study rows (cell groups of at least 30 cells), also the reference for panel E's cells.
NUC <- c("Hickey 2023", "Augsornworawat 2023")             # single-nucleus studies
# The tables name the native intestinal groups "Gut ..."; the figure names them "Intestinal ...".
display_group <- function(g) sub("^Gut", "Intestinal", g)
serotonin_rows <- function() {
  d <- read_table("serotonin_codetection.tsv")
  d <- d[d$drawn, ]
  d$grp[d$grp == "HIO EC"] <- "HIO graft EC"                # the name used in Figures 3 and 5
  d$grp <- display_group(d$grp)
  d
}
