# Figures S2 and S3 helpers: type sizes, theme, colors, the panels' places on the two pages, and grid helpers that
# place text, plots and keys on the page. Page coordinates are inches from the page's top-left corner.

FONT <- "Arial"
PT <- 7                     # body text (pt)
MM <- PT / .pt              # the same size in ggplot text units
LINE <- 1.2 * PT / 72       # one line of body text (in)
EDGE <- 0.15                # hairline point outline (ggplot stroke)
L <- .25; R <- 7.25         # page margins (in)

# where the panels sit on the page (in): Figure S2 (A, B on top; C under A; D across the page) and Figure S3 (A to C
# stacked); bottoms of the panels above, as the page lays them out
S2_A_BOTTOM <- 3.7433333333333332
S2_B_BOTTOM <- 4.8624999999999998
S2_C_TOP <- 3.9433333333333334
S2_C_BOTTOM <- 5.2400000000000002
S2_D_TOP <- 5.46
S3_A_BOTTOM <- 2.2683333333333331
S3_B_TOP <- 2.5183333333333331
S3_B_BOTTOM <- 5.1161346340113472
S3_C_TOP <- 5.3661346340113472
stopifnot(S2_C_TOP == S2_A_BOTTOM + .2, S2_D_TOP == max(S2_B_BOTTOM, S2_C_BOTTOM) + .22,
          S3_B_TOP == S3_A_BOTTOM + .25, S3_C_TOP == S3_B_BOTTOM + .25)

# colors: EC by source, SC-SIC, other endocrine, the four-marker classes, fitting groups and clusters
EC_GUT <- "#43A5B6"; EC_HIO <- "#6943AC"
OTHER_EEC <- "#8C8C8C"; UNRES <- "#D3D6D9"; PROG <- "#444444"
SC_SIC <- "#AA3377"; OTHER_SC <- "#8C8C8C"
CLS_COL <- c(HIGH_CONFIDENCE_EC = "#8E3D29", LOW_CONFIDENCE_EC = "#D8A390", INSUFFICIENT_EC_EVIDENCE = "#D9D9D9")
CLS_LAB <- c(HIGH_CONFIDENCE_EC = "Classified as EC", LOW_CONFIDENCE_EC = "Above cutoff, not classified",
             INSUFFICIENT_EC_EVIDENCE = "Below cutoff")
MARK_COL <- c("0" = "#ECE8DF", "1" = "#D8C5AA", "2" = "#CB9B72", "3" = "#B86E46", "4" = "#8E3D29")
GROUP_COL <- c("E-MTAB-8901_fetal" = "#B9CC8F", "E-MTAB-8901_paediatric" = "#4F5D1C", "GSE185224" = "#7A7A7A",
               "GSE125970" = "#E48AAE", "E-MTAB-10187_intestinal" = "#2B3F6E")
GROUP_LAB <- c("E-MTAB-8901_fetal" = "Elmentaite 2020, fetal", "E-MTAB-8901_paediatric" = "Elmentaite 2020, child",
               "GSE185224" = "Burclaff 2022", "GSE125970" = "Wang 2020", "E-MTAB-10187_intestinal" = "Yu 2021a")
CLUSTER_PAL <- c("#9E82A4", "#F0A834", "#204F41", "#3E9B80", "#573D62", "#B29DFC", "#4A7473", "#734608", "#808722",
                 "#C7AFC7", "#0DCE99", "#7F6A4D", "#B745C6", "#D9886B", "#5B8DB8", "#A3A3A3")
# test-study headings (S2B, S3A)
TITLE_B <- list(gut_cell_atlas = c("Elmentaite 2021", "fetal + adult intestine"), egozi = c("Egozi 2021", "fetal + neonatal intestine"),
                hickey = c("Hickey 2023 (nuclei)", "adult intestine"), hio = c("McCauley 2023", "HIO graft"))

theme_p2 <- function() theme_minimal(base_family = FONT, base_size = PT) + theme(
  text = element_text(family = FONT, colour = "black", size = PT),
  axis.text = element_text(colour = "black", size = PT), axis.title = element_text(colour = "black", size = PT),
  strip.text = element_blank(), panel.grid.minor = element_blank(),
  panel.grid.major = element_line(colour = "#E4E4E4", linewidth = 0.2),
  axis.line.x = element_line(colour = "black", linewidth = 0.3), axis.ticks.x = element_line(colour = "black", linewidth = 0.3),
  axis.ticks.length = unit(1.5, "pt"), legend.position = "none", plot.margin = margin(0, 0, 0, 0))

str_w <- function(s, size = PT, face = "plain") {
  pushViewport(viewport(gp = gpar(fontfamily = FONT, fontsize = size, fontface = face)))
  w <- convertWidth(stringWidth(s), "in", valueOnly = TRUE); popViewport(); w
}
gp_t <- function(size = PT, face = "plain", col = "black") gpar(fontfamily = FONT, fontsize = size, fontface = face, col = col)
# text: (x, y) is the anchor in page inches; just as in grid
txt <- function(s, x, y, just = c("left", "top"), size = PT, face = "plain", col = "black", rot = 0)
  grid.text(s, x = unit(x, "in"), y = unit(H - y, "in"), just = just, rot = rot, gp = gp_t(size, face, col))
letter <- function(l, x, y) txt(l, x, y, size = 11, face = "bold")
rect_at <- function(x, y, w, h, fill, col = NA, lwd = .5)    # (x, y) top-left
  grid.rect(x = unit(x, "in"), y = unit(H - y, "in"), width = unit(w, "in"), height = unit(h, "in"), just = c("left", "top"),
            gp = gpar(fill = fill, col = col, lwd = lwd))
seg_at <- function(x0, y0, x1, y1, col = "black", lwd = .5, lty = 1)
  grid.segments(unit(x0, "in"), unit(H - y0, "in"), unit(x1, "in"), unit(H - y1, "in"), gp = gpar(col = col, lwd = lwd, lty = lty, lineend = "butt"))
# a ggplot whose whole plot area is drawn at (x, y) top-left, w x h
draw_at <- function(p, x, y, w, h) {
  pushViewport(viewport(x = unit(x, "in"), y = unit(H - y, "in"), width = unit(w, "in"), height = unit(h, "in"), just = c("left", "top")))
  grid.draw(if (inherits(p, "gg")) ggplotGrob(p) else p); popViewport()
}
# panel box of a ggplot drawn w x h: offsets (in) of the panel from the plot's left and top
panel_offsets <- function(p, w, h) {
  g <- ggplotGrob(p)
  pushViewport(viewport(width = unit(w, "in"), height = unit(h, "in")))
  wd <- convertWidth(g$widths, "in", valueOnly = TRUE); ht <- convertHeight(g$heights, "in", valueOnly = TRUE); popViewport()
  nw <- unitType(g$widths) == "null"; nh <- unitType(g$heights) == "null"
  wd[nw] <- (w - sum(wd[!nw])) * as.numeric(g$widths[nw]) / sum(as.numeric(g$widths[nw]))
  ht[nh] <- (h - sum(ht[!nh])) * as.numeric(g$heights[nh]) / sum(as.numeric(g$heights[nh]))
  pn <- g$layout[grepl("^panel", g$layout$name), ]
  list(l = cumsum(c(0, wd))[min(pn$l)], r = cumsum(c(0, wd))[max(pn$r) + 1], t = cumsum(c(0, ht))[min(pn$t)], b = cumsum(c(0, ht))[max(pn$b) + 1])
}
# bold study over plain type above each facet column
facet_heads <- function(p, x0, y0, w, h, heads, nfac, sp) {
  o <- panel_offsets(p, w, h); pw <- (o$r - o$l - (nfac - 1) * sp) / nfac
  for (i in seq_len(nfac)) { xl <- x0 + o$l + (i - 1) * (pw + sp)
    txt(heads[[i]][1], xl, y0, face = "bold"); txt(heads[[i]][2], xl, y0 + LINE) }
  list(o = o, pw = pw)
}

# key items: list(label, w_sym, draw(x, y)); y is the row's vertical centre
k_square <- function(label, fill, col = NA) list(label = label, w_sym = .15, draw = function(x, y)
  grid.rect(x = unit(x, "in"), y = unit(H - y, "in"), width = unit(.09, "in"), height = unit(.09, "in"), just = c("left", "centre"),
            gp = gpar(fill = fill, col = col, lwd = .6)))
k_line <- function(label, col, lty = 1, lwd = 1.1) list(label = label, w_sym = .24, draw = function(x, y)
  seg_at(x, y, x + .18, y, col = col, lwd = lwd, lty = lty))
k_vdash <- function(label) list(label = label, w_sym = .1, draw = function(x, y)
  seg_at(x + .03, y - .055, x + .03, y + .055, col = "black", lwd = .9, lty = "22"))
k_point <- function(label, fill, pch = 21, size = .075) list(label = label, w_sym = .14, draw = function(x, y)
  grid.points(unit(x + .04, "in"), unit(H - y, "in"), pch = pch, size = unit(size, "in"), gp = gpar(col = "black", fill = fill, lwd = .5)))
k_text <- function(label, face = "plain") list(label = label, w_sym = 0, face = face, draw = function(x, y) invisible())
k_multi <- function(label, fills) list(label = label, w_sym = .11 * length(fills) + .05, draw = function(x, y)
  for (j in seq_along(fills)) grid.rect(x = unit(x + (j - 1) * .11, "in"), y = unit(H - y, "in"), width = unit(.09, "in"),
                                        height = unit(.09, "in"), just = c("left", "centre"), gp = gpar(fill = fills[j], col = NA)))
item_w <- function(it) it$w_sym + str_w(it$label, face = if (is.null(it$face)) "plain" else it$face)
# one row of key items, from x (or centred in [x, x + width] when width is given), row centre at y
key_row <- function(items, x, y, width = NULL, gap = .14) {
  ws <- vapply(items, item_w, 0); tot <- sum(ws) + gap * (length(items) - 1)
  if (!is.null(width)) { if (tot > width + 1e-6) stop("key row wider than its box: ", round(tot, 3), " > ", width); x <- x + (width - tot) / 2 }
  for (i in seq_along(items)) {
    it <- items[[i]]; it$draw(x, y)
    txt(it$label, x + it$w_sym, y, just = c("left", "centre"), face = if (is.null(it$face)) "plain" else it$face)
    x <- x + ws[i] + gap
  }
  invisible(tot)
}
pct_lab <- function(v) paste0(v, "%")
