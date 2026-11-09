# Drawing helpers shared by the Figure 2 supplement panels: text size, colors, theme, study names, placement on
# the page, and the layout helpers that several panels of one page share. Sourced by
# workflow/scripts/common/draw_panel.R, which opens the page device and defines its size W x H (inches) and DPI.

PT <- 7.5; FONT <- "Arial"        # body text, points
tx <- PT / ggplot2::.pt          # body text in ggplot's text size units
MARKS <- c("0" = "#ECE8DF", "1" = "#D8C5AA", "2" = "#CB9B72", "3" = "#B86E46", "4" = "#8E3D29")   # markers detected, 0-4
RULE <- "#B86E46"; PANEL <- "#2A7F8E"; PANEL_LIGHT <- "#A9D3DA"   # four-marker classifier rust; intestinal 20-gene classifier dark teal (after intestinal EC teal)
TPH1_COL <- "#ABABAB"                                             # TPH1-only classifier grey
CLASS <- c(EXPLICIT_EC = "#43A5B6", TERMINAL_NON_EC_EEC = "#8C8C8C", BROAD_EEC_PARENT_UNRESOLVED = "#D3D6D9",
           ENDOCRINE_PROGENITOR = "#444444")                    # primary study label classes
CLASS_LABEL <- c(EXPLICIT_EC = "EC (primary study label)", TERMINAL_NON_EC_EEC = "Other terminal EEC",
                 BROAD_EEC_PARENT_UNRESOLVED = "Unresolved EEC", ENDOCRINE_PROGENITOR = "Progenitor")
EC_SRC <- c(gut = "#43A5B6", induced = "#CCBB44", hio = "#6943AC")   # EC labels by the cells' source
CALLS <- c(HIGH_CONFIDENCE_EC = "#8E3D29", LOW_CONFIDENCE_EC = "#D8A390", INSUFFICIENT_EC_EVIDENCE = "#D9D9D9")   # four-marker classes
GRID_COL <- "#E4E4E4"; REF <- "#777777"
# S2A-B: the 13 clusters of the fitting cells, and the names of the four-marker classes
CLUSTER_COL <- c("0" = "#9E82A4", "1" = "#F0A834", "2" = "#204F41", "3" = "#3E9B80", "4" = "#573D62", "5" = "#B29DFC",
                 "6" = "#4A7473", "7" = "#734608", "8" = "#808722", "9" = "#C7AFC7", "10" = "#0DCE99", "11" = "#7F6A4D", "12" = "#B745C6")
CLASS_KEY <- c(HIGH_CONFIDENCE_EC = "High-confidence EC (classified as EC)", LOW_CONFIDENCE_EC = "Low-confidence EC (not classified as EC)",
               INSUFFICIENT_EC_EVIDENCE = "Insufficient EC evidence")

theme_kit <- function(grid = c("x", "y", "none"), legend = "none") {
  grid <- match.arg(grid)
  gl <- element_line(color = GRID_COL, linewidth = .2)
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
          panel.grid.major.x = if (grid == "x") gl else element_blank(),
          panel.grid.major.y = if (grid == "y") gl else element_blank(),
          axis.line.x = element_line(color = "black", linewidth = .35),
          axis.ticks.x = element_line(color = "black", linewidth = .3),
          axis.ticks.length = unit(1.3, "pt"),
          panel.spacing = unit(.16, "in"),
          plot.margin = margin(2, 6, 2, 2))
}

# Width of the widest string in inches, on the open device.
str_in <- function(s, pt = PT, bold = FALSE, face = if (bold) "bold" else "plain") {
  pushViewport(viewport(gp = gpar(fontfamily = FONT, fontsize = pt, fontface = face)))
  on.exit(popViewport())
  max(vapply(s, function(z) convertWidth(stringWidth(z), "in", valueOnly = TRUE), 0))
}
in_w <- function(u) convertWidth(u, "in", valueOnly = TRUE)
in_h <- function(u) convertHeight(u, "in", valueOnly = TRUE)

# Study names (config/study_names.tsv): first author and year, short name, sample type. Fitting groups, test
# studies and the HIO groups are looked up by the group or study keys used in the result tables.
SN <- read.delim(file.path(REPO, "config", "study_names.tsv"), check.names = FALSE, colClasses = "character",
                 na.strings = character(0), fileEncoding = "UTF-8")
UNIT_OF <- c("E-MTAB-10187_intestinal" = "SF_YU2021_ENDODERM:fetal", "E-MTAB-8901_fetal" = "SF_ELMENTAITE2020_DEVGUT:fetal",
             "E-MTAB-8901_paediatric" = "SF_ELMENTAITE2020_DEVGUT:child", "GSE125970" = "SF_GSE125970_ADULT_INTESTINE",
             "GSE185224" = "SF_GSE185224_ADULT_REGIONAL_INTESTINE", "GSE214852_WT" = "SF_GSE214852_HIO_GRAFT",
             "E-MTAB-10268" = "SF_YU2021_ENDODERM:hio", gut_cell_atlas = "SF_ELMENTAITE2021_GUT_ATLAS",
             hickey = "SF_HICKEY2023_ADULT_EEC", egozi = "SF_EGOZI2021_DEVELOPMENTAL_EEC", gse238274 = "SF_INTESTINAL_GSE238274",
             gse214852_hio = "SF_GSE214852_HIO_GRAFT")
sn <- function(x, col) {
  v <- SN[[col]][match(ifelse(x %in% names(UNIT_OF), UNIT_OF[x], x), SN$unit_id)]
  if (anyNA(v)) stop("no study_names.tsv row for ", paste(x[is.na(v)], collapse = ", "))
  v
}
study <- function(x) sn(x, "author_year")           # "Hickey 2023"
study_short <- function(x) sn(x, "short_name")      # "Elmentaite 2020 fetal"
study_type <- function(x) sn(x, "type_words")       # "adult intestine"
fmt <- function(v) format(v, big.mark = ",", trim = TRUE)
# the four test studies and the McCauley 2023 HIO grafts scored by both classifiers, in drawing order
SCORED <- c("SF_ELMENTAITE2021_GUT_ATLAS", "SF_HICKEY2023_ADULT_EEC", "SF_EGOZI2021_DEVELOPMENTAL_EEC", "SF_INTESTINAL_GSE238274",
            "SF_GSE214852_HIO_GRAFT")

# ---- S3A-C: gtable column grid, column headers, hand-drawn keys and vertically stacked blocks ----
GRID_BC <- list(W = 7.3, G = 0.58, S = 0.12, R = 0.10)              # S3B-C: total, left gutter, gap between studies, right
GRID_BC$P <- (GRID_BC$W - GRID_BC$G - 4 * GRID_BC$S - GRID_BC$R) / 5  # panel width: S3B and S3C put each study at the same x
HEAD_BOLD_PT <- 8; HEAD_GAP_PT <- 10.4; HEAD_H <- 0.29; RULE_STOP <- 0.12
TOP <- .47                                                            # band above each panel's content (in)
PAD_B <- .03                                                          # blank strip at the foot of each panel (in)
panel_cols <- function(g) sort(unique(g$layout$l[grepl("^panel", g$layout$name)]))
panel_rows <- function(g) sort(unique(g$layout$t[grepl("^panel", g$layout$name)]))
total_h <- function(g) sum(vapply(seq_along(g$heights), function(i) in_h(g$heights[i]), 0))
other_h <- function(g) { pr <- panel_rows(g); sum(vapply(setdiff(seq_along(g$heights), pr), function(i) in_h(g$heights[i]), 0)) }
# every gtable column absolute: G left of the first panel, P per panel, S per gap, the rest on the right
fix_columns <- function(g, G, P, S, W) {
  pc <- panel_cols(g); n <- length(pc); nc <- length(g$widths)
  if (n > 1 && length(S) == 1) S <- rep(S, n - 1)
  wid <- vapply(seq_len(nc), function(i) if (i %in% pc) 0 else in_w(g$widths[i]), 0)
  left <- seq_len(pc[1] - 1); right <- (pc[n] + 1):nc; R <- W - G - n * P - sum(S[seq_len(n - 1)])
  wid[left[1]] <- wid[left[1]] + G - sum(wid[left])
  wid[right[length(right)]] <- wid[right[length(right)]] + R - sum(wid[right])
  wid[pc] <- P; if (n > 1) wid[pc[-n] + 1] <- S
  g$widths <- unit(wid, "in"); g
}
fix_rows <- function(g, panel_h) {                                    # panel rows get panel_h, other rows their natural height
  pr <- panel_rows(g)
  hh <- vapply(seq_along(g$heights), function(i) if (i %in% pr) 0 else in_h(g$heights[i]), 0)
  hh[pr] <- panel_h; g$heights <- unit(hh, "in"); g
}
# column header: bold study over plain sample type, from the panel's left edge; a spanning header gets a thin line
header_grob <- function(bold, plain, rule_to = NULL) {
  kids <- gList(textGrob(bold, x = unit(0, "npc"), y = unit(1, "npc"), just = c("left", "top"),
                         gp = gpar(fontfamily = FONT, fontsize = HEAD_BOLD_PT, fontface = "bold", col = "black")))
  pl <- if (is.list(plain)) plain else list(list(x = 0, s = plain))
  for (q in pl) kids <- gList(kids, textGrob(q$s, x = unit(q$x, "in"), y = unit(1, "npc") - unit(HEAD_GAP_PT, "pt"),
                                            just = c("left", "top"), gp = gpar(fontfamily = FONT, fontsize = PT, col = "black")))
  if (!is.null(rule_to))
    kids <- gList(kids, segmentsGrob(x0 = unit(str_in(bold, HEAD_BOLD_PT, TRUE) + .05, "in"), x1 = unit(rule_to, "in"),
                                     y0 = unit(1, "npc") - unit(3.3, "pt"), y1 = unit(1, "npc") - unit(3.3, "pt"), gp = gpar(col = REF, lwd = .6)))
  gTree(children = kids)
}
add_headers <- function(g, heads, above = 0) {                        # a header row above the first panel row (or `above` rows higher)
  pr <- panel_rows(g); pc <- panel_cols(g); pos <- pr[1] - 1 - above
  g <- gtable::gtable_add_rows(g, unit(HEAD_H, "in"), pos = pos)
  for (k in seq_along(heads)) {
    hd <- heads[[k]]
    span_w <- sum(in_w(g$widths[pc[hd$l]:pc[hd$r]]))
    g <- gtable::gtable_add_grob(g, header_grob(hd$bold, hd$plain, rule_to = if (isTRUE(hd$rule)) span_w - RULE_STOP else NULL),
                                 t = pos + 1, l = pc[hd$l], r = pc[hd$r], clip = "off", name = paste0("head-", k))
  }
  g
}
# hand-drawn key: rows of entries list(label, glyph, col), each row centered; grid = TRUE aligns columns across rows
draw_key <- function(rows, gap = .2, glyph_w = .16, grid = FALSE) {
  function() {
    vw <- convertWidth(unit(1, "npc"), "in", valueOnly = TRUE); nr <- length(rows)
    ww <- lapply(rows, function(r) vapply(r, function(e) str_in(e$label), 0) + glyph_w)
    if (grid) {
      nc <- max(lengths(rows)); cw <- vapply(seq_len(nc), function(j) max(vapply(ww, function(w) if (j <= length(w)) w[j] else 0, 0)), 0)
      ww <- lapply(rows, function(r) cw[seq_along(r)]); tot <- sum(cw) + gap * (nc - 1)
    }
    for (i in seq_len(nr)) {
      y <- unit((nr - i + .5) / nr, "npc")
      x <- (vw - (if (grid) tot else sum(ww[[i]]) + gap * (length(ww[[i]]) - 1))) / 2
      for (j in seq_along(rows[[i]])) {
        e <- rows[[i]][[j]]; gx <- x + .06
        ln <- function(lty = 1) grid.lines(x = unit(c(gx - .06, gx + .06), "in"), y = unit.c(y, y), gp = gpar(col = e$col, lwd = 1.4, lty = lty, lineend = "butt"))
        pt <- function(pch = 16, size = 4.6, fill = NA) grid.points(x = unit(gx, "in"), y = y, pch = pch, size = unit(size, "pt"), gp = gpar(col = e$col, fill = fill, lwd = 1))
        switch(e$glyph,
               point = pt(),
               line = ln(),
               linepoint = { ln(); pt(size = 3.6) },
               square = grid.rect(x = unit(gx, "in"), y = y, width = unit(.11, "in"), height = unit(.11, "in"), gp = gpar(fill = e$col, col = NA)),
               dashv = grid.lines(x = unit(c(gx, gx), "in"), y = unit.c(y - unit(.06, "in"), y + unit(.06, "in")),
                                  gp = gpar(col = e$col, lwd = 1, lty = "22", lineend = "butt")),
               diamond = pt(pch = 23, size = 4.2, fill = "white"))
        grid.text(e$label, x = unit(x + glyph_w, "in"), y = y, just = c("left", "centre"), gp = gpar(fontfamily = FONT, fontsize = PT, col = "black"))
        x <- x + ww[[i]][j] + gap
      }
    }
  }
}
# blocks of fixed height, stacked top to bottom: list(list(h = inches, draw = function, gtable or NULL for a gap))
stack_pieces <- function(pieces) {
  function() {
    y <- convertHeight(unit(1, "npc"), "in", valueOnly = TRUE)
    for (p in pieces) {
      pushViewport(viewport(x = unit(0, "in"), y = unit(y, "in"), width = unit(1, "npc"), height = unit(p$h, "in"), just = c("left", "top")))
      if (inherits(p$draw, "gtable")) grid.draw(p$draw) else if (is.function(p$draw)) p$draw()
      popViewport(); y <- y - p$h
    }
  }
}

# ---- S4 ----
# Height of a single-panel ggplot's plotting area, in inches, when drawn `content_h` inches tall.
panel_h <- function(p, content_h) {
  hs <- ggplotGrob(p)$heights
  content_h - sum(vapply(seq_along(hs), function(i) if (unitType(hs[i]) == "null") 0 else in_h(hs[i]), 0))
}
# S4C and S4D share their rows: one per scored study, bold name over the plain sample type, left of the plot; the dotted
# line at YDOT separates the studies whose primary study labels rank Fig. 2D's genes and train the classifier (above)
both <- function(gut, hio) { x <- read_table(gut); rbind(x, read_table(hio)[names(x)]) }   # test studies, then McCauley 2023
YR <- setNames(rev(seq_along(SCORED)), SCORED)
YDOT <- 2.5
rowlab <- data.frame(y = unname(YR), b = paste0(study(SCORED), ifelse(SCORED == "SF_HICKEY2023_ADULT_EEC", " (nuclei)", ""), "  "),
                     p = paste0(study_type(SCORED), "  "))
LW <- max(str_in(rowlab$b, PT, TRUE), str_in(rowlab$p, PT))          # row label width (in)

# Draw a panel at its place on the page. The panel's rectangle is x, y (top-left corner) and w_px x h_px, in page
# pixels; its content, w inches wide and h tall, is drawn `bottom` inches above the rectangle's foot, which leaves
# a 0.06 in margin above it. `draw` is a function drawing into the content's viewport, or a ggplot.
place <- function(draw, x, y, w_px, h_px, w, h, bottom = 0) {
  pushViewport(viewport(x = unit(x / DPI, "in"), y = unit((H * DPI - y - h_px) / DPI, "in"), width = unit(w_px / DPI, "in"),
                        height = unit(h_px / DPI, "in"), just = c("left", "bottom")))
  pushViewport(viewport(x = unit(0, "in"), y = unit(bottom, "in"), width = unit(w, "in"), height = unit(h, "in"),
                        just = c("left", "bottom")))
  if (is.function(draw)) draw() else print(draw, newpage = FALSE)
  popViewport(2)
}
