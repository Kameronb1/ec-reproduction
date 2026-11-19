# Drawing helpers shared by the Supplementary Figure S8 panels (the same as those of Figure 4): text size, colors, theme,
# text and key glyphs placed in page inches (measured from the page's top-left corner), panel positions of a drawn
# ggplot, and the arrow marks. Sourced by workflow/scripts/common/draw_panel.R before the panel's plot script.

W <- 7.5; H <- 10; DPI <- 600; PT <- 7.5                      # page in inches, PNG resolution, base text size (pt)
open_png <- function(file, w = W, h = H, res = DPI)            # the device the panels are measured on (as the page)
  png(file, width = round(w * res), height = round(h * res), units = "px", res = res, type = "cairo", bg = "white")
POP_COL <- c(SC_SIC = "#AA3377", SC_BETA = "#4477AA", SC_ALPHA = "#117733")
POP_LAB <- c(SC_SIC = "SC-SIC", SC_BETA = "SC-beta", SC_ALPHA = "SC-alpha")
GCOL <- c(`SC-SIC` = "#AA3377", `SC-alpha` = "#117733", `SC-beta` = "#4477AA")
ARW_LEN <- .042; ARW_ANGLE <- 25                               # every arrowhead: .042 in sides, 25 degrees, closed
ARW <- arrow(length = unit(ARW_LEN, "in"), type = "closed", angle = ARW_ANGLE)
HEAD_AX0 <- ARW_LEN * cos(ARW_ANGLE * pi / 180)                # the head's length along its axis (in)
HEAD_MIN_IN <- ARW_LEN + 2 / DPI                               # shorter changes are drawn as a bare arrowhead
STUB_IN <- 1 / (5 * DPI)                                       # a bare head sits on a 0.2 px stub
PT_SIZE <- .95; PT_STROKE <- .45                               # the open start circle
CIRCLE_PT <- PT_SIZE * ggplot2::.pt + PT_STROKE * ggplot2::.stroke / 2; CIRCLE_LWD <- PT_STROKE * ggplot2::.stroke / 2
CIRCLE_R <- .375 * CIRCLE_PT / 72 + CIRCLE_LWD / 96 / 2        # its outer radius (in)
BAND_GREY <- "#E8E8E8"                                         # row bands
BAND_ALPHA <- 1 - grDevices::col2rgb(BAND_GREY)[1, 1] / 255    # row bands: black at this alpha renders BAND_GREY

th <- function() theme_minimal(base_family = "Arial", base_size = PT) + theme(
  text = element_text(color = "black", size = PT), axis.text = element_text(color = "black", size = 7),
  axis.title = element_text(color = "black", size = 7), plot.title = element_text(color = "black", size = 7.5, face = "bold", margin = margin(0, 0, 2, 0)),
  strip.text = element_text(color = "black", size = 7.5, face = "bold", margin = margin(1, 0, 2, 0)),
  panel.grid.minor = element_blank(), panel.grid.major.y = element_blank(), panel.grid.major.x = element_line(color = "#E4E4E4", linewidth = .2),
  axis.line.x = element_line(color = "black", linewidth = .35), axis.ticks.x = element_line(color = "black", linewidth = .3),
  axis.ticks.length = unit(1.3, "pt"), plot.margin = margin(2, 4, 2, 2), legend.position = "none", plot.title.position = "panel")
nobg <- function(p) p + theme(plot.background = element_blank())   # the page's row bands show through the plots

# ---- text metrics on the open device (inches) -----------------------------------------------------------------
fm <- function(s, size, face, f) { pushViewport(viewport(gp = gpar(fontfamily = "Arial", fontsize = size, fontface = face))); v <- f(s); popViewport(); v }
str_in <- function(s, size = PT, face = "plain") fm(s, size, face, function(s) convertWidth(stringWidth(s), "in", valueOnly = TRUE))
asc_in <- function(s, size = PT, face = "plain") fm(s, size, face, function(s) convertHeight(unit(1, "strheight", s), "in", valueOnly = TRUE))
dsc_in <- function(s, size = PT, face = "plain") fm(s, size, face, function(s) convertHeight(unit(1, "strdescent", s), "in", valueOnly = TRUE))

# ---- drawing in page inches --------------------------------------------------------------------------------------
text_at <- function(s, x, y, size = PT, face = "plain", col = "black")      # text with its top at y
  grid.text(s, x = unit(x, "in"), y = unit(H - y, "in"), just = c("left", "top"), gp = gpar(fontfamily = "Arial", fontsize = size, fontface = face, col = col))
text_bl <- function(s, x, y, size = 7, face = "plain", just = "left")      # text on the baseline y
  grid.text(s, x = unit(x, "in"), y = unit(H - y, "in"), just = c(just, "bottom"), gp = gpar(fontfamily = "Arial", fontsize = size, fontface = face, col = "black"))
draw_at <- function(g, x, y, w, h) { pushViewport(viewport(x = x / W, y = 1 - y / H, width = w / W, height = h / H, just = c("left", "top")))
  grid.draw(g); popViewport() }

# ---- ggplot layout -----------------------------------------------------------------------------------------------
align_h <- function(ps) { gs <- lapply(ps, ggplotGrob); h <- do.call(unit.pmax, lapply(gs, function(g) g$heights)); lapply(gs, function(g) { g$heights <- h; g }) }
nonpanel_w <- function(g) { w <- g$widths; sum(convertWidth(w[unitType(w) != "null"], "in", valueOnly = TRUE)) }
nonpanel_h <- function(g) { h <- g$heights; sum(convertHeight(h[unitType(h) != "null"], "in", valueOnly = TRUE)) }
panels_in <- function(g, x, y, w, h) {                     # every panel of gtable g drawn at (x, y, w, h), in page inches
  L <- g$layout[grepl("^panel", g$layout$name), ]; L <- L[order(L$t, L$l), ]
  nw <- unitType(g$widths) == "null"; nh <- unitType(g$heights) == "null"
  free_w <- w - sum(convertWidth(g$widths[!nw], "in", valueOnly = TRUE)); free_h <- h - sum(convertHeight(g$heights[!nh], "in", valueOnly = TRUE))
  span <- function(u, idx, isnull, free, conv) { if (!length(idx)) return(0); uu <- u[idx]; nn <- isnull[idx]
    (if (any(!nn)) sum(conv(uu[!nn], "in", valueOnly = TRUE)) else 0) + (if (any(nn)) sum(as.numeric(uu[nn])) * free / sum(as.numeric(u[isnull])) else 0) }
  lapply(seq_len(nrow(L)), function(i) list(name = L$name[i], t = L$t[i], l = L$l[i],
    left_in = x + span(g$widths, seq_len(L$l[i] - 1L), nw, free_w, convertWidth), top_in = y + span(g$heights, seq_len(L$t[i] - 1L), nh, free_h, convertHeight),
    width_in = as.numeric(g$widths[L$l[i]]) * free_w / sum(as.numeric(g$widths[nw])), height_in = as.numeric(g$heights[L$t[i]]) * free_h / sum(as.numeric(g$heights[nh])))) }

# ---- marks -------------------------------------------------------------------------------------------------------
small_head <- function(z, col, ycol) {                    # a change too short for an arrow: a bare head, tip at the end value
  if (!nrow(z)) return(NULL)
  z$htip <- z$to; z$hbase <- z$to - z$stub * sign(z$to - z$from)
  geom_segment(data = z, aes(x = hbase, xend = htip, y = .data[[ycol]], yend = .data[[ycol]]), color = col, linewidth = .45, arrow = ARW) }
arrows_ef <- function(d, col) {                           # shaft, head on a stub, start circle on top; thr and stub in data units
  d$tip0 <- d$to - d$stub * sign(d$to - d$from)
  a <- d[abs(d$to - d$from) >= d$thr, ]; z <- d[abs(d$to - d$from) < d$thr, ]
  list(geom_blank(data = d, aes(x = from, y = y)), geom_blank(data = d, aes(x = to, y = y)),
       if (nrow(a)) geom_segment(data = a, aes(x = from, xend = to, y = y, yend = y), color = col, linewidth = .45),
       if (nrow(a)) geom_segment(data = a, aes(x = tip0, xend = to, y = y, yend = y), color = col, linewidth = .45, arrow = ARW),
       if (nrow(a)) geom_point(data = a, aes(x = from, y = y), shape = 21, fill = "white", color = col, size = PT_SIZE, stroke = PT_STROKE),
       small_head(z, col, "y")) }

# ---- key glyphs (y = the key line's top; glyphs centred .045 in below it, text top .005 in above it) ----------------
key_line <- function(x, y, col, lwd = 1.3, lty = 1, len = .2, butt = FALSE)
  grid.lines(unit(c(x, x + len), "in"), unit(rep(H - y - .045, 2), "in"), gp = if (butt) gpar(col = col, fill = col, lwd = lwd, lty = lty, lineend = "butt") else gpar(col = col, lwd = lwd, lty = lty))
START_GRAFT_KEY <- c(start = "Start", end = "Graft / extended culture")
start_graft_key <- function(x, y) {                       # open circle 'Start', arrow 'Graft / extended culture'
  grid.points(x = unit(x, "in"), y = unit(H - y - .045, "in"), pch = 21, size = unit(CIRCLE_PT, "pt"), gp = gpar(col = "#555555", fill = "white", lwd = CIRCLE_LWD))
  text_at(START_GRAFT_KEY[["start"]], x + .12, y - .005)
  x2 <- x + .12 + str_in(START_GRAFT_KEY[["start"]]) + .2
  grid.lines(x = unit(c(x2, x2 + .18), "in"), y = unit(rep(H - y - .045, 2), "in"), arrow = ARW, gp = gpar(col = "#555555", fill = "#555555", lwd = 1.3))
  text_at(START_GRAFT_KEY[["end"]], x2 + .23, y - .005); x2 + .23 + str_in(START_GRAFT_KEY[["end"]]) }
population_key <- function(x, y, pops = c("SC_SIC", "SC_BETA", "SC_ALPHA")) { for (k in pops) {   # coloured line and name per population
  key_line(x, y, POP_COL[[k]]); text_at(POP_LAB[[k]], x + .26, y - .005); x <- x + str_in(POP_LAB[[k]]) + .26 + .16 }; x - .16 }
