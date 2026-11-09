# Drawing helpers shared by the Figure S6-S8 and S15 panels: text size, theme, study names, string widths and placement on
# the page. Sourced by workflow/scripts/common/draw_panel.R, which defines the page size W x H (inches) and DPI.

PT <- 7.5; FONT <- "Arial"                # body text, points
REF <- "#777777"; SIC <- "#AA3377"        # reference lines; SC-SIC
SCC <- "#882255"; SCC_LIGHT <- "#D9A6BF"   # SC-SIC 20-gene classifier wine (after SC-SIC magenta)
# The six SC-islet studies, in Figure 3A's order, by first author and year
STUDY <- c(GSE151117 = "Augsornworawat 2020", GSE167880 = "Balboa 2022", GSE202497 = "Zhu 2023", GSE226346 = "Siehler 2024",
           GSE185036 = "Du 2022", GSE199636 = "Augsornworawat 2023 (nuclei)")
minus <- scales::label_number(style_negative = "minus")    # axis numbers with a true minus sign

theme_rk <- function(grid = c("x", "y", "both", "none"), legend = "none") {
  grid <- match.arg(grid)
  gl <- element_line(color = "#E4E4E4", linewidth = .2)
  theme_minimal(base_family = FONT, base_size = PT) +
    theme(text = element_text(color = "black", size = PT),
          axis.text = element_text(color = "black", size = PT),
          axis.title = element_text(color = "black", size = PT),
          strip.text = element_text(color = "black", size = 8, face = "bold", hjust = 0),
          plot.title = element_text(color = "black", size = 8, face = "bold"),
          plot.subtitle = element_text(color = "black", size = PT),
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

# Widths of strings in inches on the open device; str_w_at measures them on a small cairo device at `res` dpi, the
# resolution the layout of a panel is measured at (text widths differ slightly between resolutions).
str_w <- function(s, pt = PT, face = "plain") {
  pushViewport(viewport(gp = gpar(fontfamily = FONT, fontsize = pt, fontface = face))); on.exit(popViewport())
  vapply(s, function(z) convertWidth(stringWidth(z), "in", valueOnly = TRUE), 0)
}
str_w_at <- function(s, res, pt = PT, face = "plain") {
  png(tempfile(fileext = ".png"), width = 2, height = 2, units = "in", res = res, type = "cairo"); on.exit(dev.off())
  str_w(s, pt, face)
}

# Total size (inches) of a gtable's fixed widths or heights: everything but its stretchable plot panels.
fixed_part <- function(u, conv) sum(vapply(which(unitType(u) != "null"), function(i) conv(u[i], "in", valueOnly = TRUE), 0))

# Draw a panel as the figure page is built: on its own 600 dpi cairo device on white, as.integer(w * DPI) x ph pixels,
# with draw() filling a viewport w x h inches anchored at the device's bottom-left corner; that image is then placed on
# the page with its top-left corner at page pixel (x, y). Drawn straight onto the page, raster images (color bars,
# tints) would be sampled on a different sub-pixel grid. The page's panels lie on whole pixels: rows start 0.08 in below
# the page top and 0.10 in below the row above, and every panel follows a 0.20 in letter column. Each plot height h is
# written with the page's own arithmetic (a panel height less the band above its plot), so it is the same double.
in_place <- function(x, y, w, h, ph, draw) {
  f <- tempfile(fileext = ".png")
  png(f, width = as.integer(w * DPI), height = ph, units = "px", res = DPI, type = "cairo", bg = "white")
  grid.newpage()
  pushViewport(viewport(x = unit(0, "in"), y = unit(0, "in"), width = unit(w, "in"), height = unit(h, "in"), just = c("left", "bottom")))
  draw()
  invisible(dev.off())
  img <- png::readPNG(f); unlink(f)
  grid.raster(img, x = unit(x / DPI, "in"), y = unit(H - y / DPI, "in"), width = unit(dim(img)[2] / DPI, "in"),
              height = unit(dim(img)[1] / DPI, "in"), just = c("left", "top"), interpolate = TRUE)
}

# ---- Figure S15 A and B share one grid so their gene rows line up across the page: row pitch P, column width CW, a
# header band HEAD above the gene grid, row headings from XL, gene names right-aligned at XG, study columns from XT, and
# study names at ANG degrees, TH high (all inches).
HEAD <- .22; P <- .16; CW <- .43
XL <- .06; XG <- .58; XT <- .62; ANG <- 55; TH <- .11
NS <- length(STUDY)
LBLH <- max(str_w_at(STUDY, 100)) * sin(ANG * pi / 180) + TH * cos(ANG * pi / 180)   # height of the angled study names
tx <- function(...) geom_text(..., inherit.aes = FALSE, size = PT / .pt, family = FONT, colour = "black")
