# Drawing helpers shared by the Figure 5 panels: text size, colors, theme, placement on the page, keys, the page
# rows, and the plots that are drawn in more than one panel image. Sourced by workflow/scripts/common/draw_panel.R,
# which defines the page size W x H (inches). Every position is in inches from the page's top-left corner.
#
# Plots side by side in a row share their plot-area top and bottom (align_h), so B, C and D are built together. Panel
# A's gene names reach below the top of the B-D row; the panel images divide the page there, so the parts of A below
# that line are drawn with B and C.

PT <- 7.5                                           # body text, points
TM <- 7 / ggplot2::.pt                              # 7 pt in ggplot's text size units
pct <- function(d) 100 * (2^d - 1)                  # log2 difference -> % more

# Gene groups G1-G3 and cell types, colored as in Figures 3 and 4
GCOL <- c(G1 = "#43A5B6", G2 = "#A89A2E", G3 = "#B4553F")
GCOL_LIGHT <- vapply(GCOL, function(c) grDevices::rgb(t(.4 * grDevices::col2rgb(c) + .6 * 255), maxColorValue = 255), "")   # 40% shade (5C bars)
GL <- c("G1: EC genes low in SC-SIC", "G2: EC genes SC-SIC have", "G3: SC-SIC's own program")
POP_COL <- c("SC-SIC" = "#AA3377", "SC-beta" = "#4477AA", "SC-alpha" = "#117733")
BAND <- "#D9D9D9"                                   # random-range bars and Figure 5E's row bands
AUG <- "Augsornworawat 2023"; ZHU <- "Zhu 2023"

# ggplot theme for panels A-F
th <- function(grid_x = TRUE, grid_y = FALSE) theme_minimal(base_family = "Arial", base_size = PT) + theme(
  text = element_text(color = "black", size = PT), axis.text = element_text(color = "black", size = 7),
  axis.title = element_text(color = "black", size = 7), plot.title = element_text(color = "black", size = 7.5, face = "bold", margin = margin(0, 0, 2, 0)),
  strip.text = element_text(color = "black", size = 7, face = "bold", margin = margin(1, 0, 2, 0)),
  panel.grid.minor = element_blank(),
  panel.grid.major.x = if (grid_x) element_line(color = "#E4E4E4", linewidth = .2) else element_blank(),
  panel.grid.major.y = if (grid_y) element_line(color = "#E4E4E4", linewidth = .2) else element_blank(),
  axis.line.x = element_line(color = "black", linewidth = .35), axis.ticks.x = element_line(color = "black", linewidth = .3),
  axis.ticks.length = unit(1.3, "pt"), plot.margin = margin(2, 4, 2, 2), legend.position = "none", plot.title.position = "panel")
NOXL <- theme(axis.line.x = element_blank(), axis.ticks.x = element_blank())
zero_h <- geom_hline(yintercept = 0, color = "#777777", linewidth = .3)
zero_v <- geom_vline(xintercept = 0, color = "#777777", linewidth = .3)
# The plots are drawn without a white background: a plot box can reach over the letter line beside its letter.
NOBG <- theme(plot.background = element_blank())

# Text with its top-left corner (or other `just` point) at (x, y); vp_h is the height of the current viewport.
text_at <- function(s, x, y, size = PT, bold = FALSE, col = "black", just = c("left", "top"), vp_h = H)
  grid.text(s, x = unit(x, "in"), y = unit(vp_h - y, "in"), just = just,
            gp = gpar(fontfamily = "Arial", fontsize = size, fontface = if (bold) "bold" else "plain", col = col))
# Width of a string in inches on the open device.
str_in <- function(s, size = PT) {
  pushViewport(viewport(gp = gpar(fontfamily = "Arial", fontsize = size)))
  v <- convertWidth(stringWidth(s), "in", valueOnly = TRUE)
  popViewport()
  v
}

# Draw a ggplot or gtable with its top-left corner at (x, y), w wide and h high.
draw_at <- function(p, x, y, w, h) {
  pushViewport(viewport(x = x / W, y = 1 - y / H, width = w / W, height = h / H, just = c("left", "top")))
  grid.draw(if (inherits(p, "gtable")) p else ggplotGrob(p))
  popViewport()
}
# ggplot grobs with equal row heights, so plots side by side share their plot-area top and bottom
align_h <- function(ps) {
  gs <- lapply(ps, ggplotGrob)
  h <- do.call(unit.pmax, lapply(gs, function(g) g$heights))
  lapply(gs, function(g) { g$heights <- h; g })
}
# Paint the page white outside the rectangle from (x0, y0) to (x1, y1), whose edges lie on whole pixels. The panel
# images cut the page into rectangles; A's gene names cross the cut between A and the B-D row. (Clipping cannot do this: the
# plots' axis and strip cells switch clipping off.)
white_outside <- function(x0, y0, x1, y1) {
  paint <- function(x, y, w, h) if (w > 0 && h > 0)
    grid.rect(unit(x, "in"), unit(H - y, "in"), unit(w, "in"), unit(h, "in"), just = c("left", "top"), gp = gpar(fill = "white", col = NA))
  paint(0, 0, W, y0); paint(0, y1, W, H - y1); paint(0, y0, x0, y1 - y0); paint(x1, y0, W - x1, y1 - y0)
}

# One row of key items (type sq, ln, pt, dash, vbar, bars3 or pb, each with a label), left-aligned or centered in
# [left, left + width]; pb is a gene group as 5B and 5C draw it, its light bar (it$bar) with its point (it$fill) on it.
# A shared key's dotted outline: centered on cx, w wide, from y1 down to y2 (inches from the page top).
key_frame <- function(cx, w, y1, y2)
  grid.rect(unit(cx, "in"), unit(H - (y1 + y2) / 2, "in"), width = unit(w, "in"), height = unit(y2 - y1, "in"),
            gp = gpar(fill = NA, col = "grey45", lwd = .6, lty = "dotted"))
row_keys <- function(items, left, width, y, gap = .16, center = TRUE) {
  ww <- vapply(items, function(it) str_in(it$label) + if (it$type %in% c("ln", "bars3")) .26 else .14, 0)
  total <- sum(ww) + gap * (length(ww) - 1)
  x <- if (center) left + (width - total) / 2 else left
  for (i in seq_along(items)) {
    it <- items[[i]]; yy <- H - y - .045
    if (it$type == "sq") grid.rect(unit(x, "in"), unit(yy, "in"), width = unit(.095, "in"), height = unit(.095, "in"),
                                   just = "left", gp = gpar(fill = it$col, col = NA))
    if (it$type == "ln") grid.lines(unit(c(x, x + .2), "in"), unit(rep(yy, 2), "in"),
                                    gp = gpar(col = it$col, lwd = it$lwd %||% 1.3, lineend = "butt"))
    if (it$type == "pt") grid.points(unit(x + .05, "in"), unit(yy, "in"), pch = it$pch %||% 16, size = unit(it$size %||% 1.5, "mm"),
                                     gp = gpar(col = it$col, fill = it$fill %||% it$col, lwd = it$plwd %||% .8))
    if (it$type == "dash") grid.lines(unit(c(x + .005, x + .095), "in"), unit(rep(yy, 2), "in"),
                                      gp = gpar(col = it$col, lwd = it$lwd %||% 1.3, lineend = "butt"))
    if (it$type == "vbar") grid.rect(unit(x + .05, "in"), unit(yy, "in"), width = unit(it$w %||% .065, "in"), height = unit(it$h %||% .12, "in"),
                                     gp = gpar(fill = it$col, col = NA))
    if (it$type == "bars3") for (k in 1:3)          # three short bars, one per group color (5C's bars)
      grid.rect(unit(x + .025 + (k - 1) * .07, "in"), unit(yy, "in"), width = unit(.05, "in"), height = unit(it$h %||% .12, "in"),
                gp = gpar(fill = it$col[k], col = NA))
    if (it$type == "pb") {
      grid.rect(unit(x + .05, "in"), unit(yy, "in"), width = unit(.05, "in"), height = unit(.13, "in"), gp = gpar(fill = it$bar, col = NA))
      grid.points(unit(x + .05, "in"), unit(yy, "in"), pch = 21, size = unit(1.5, "mm"), gp = gpar(col = "black", fill = it$fill, lwd = it$plwd))
    }
    text_at(it$label, x + if (it$type %in% c("ln", "bars3")) .26 else .14, y - .005)
    x <- x + ww[i] + gap
  }
  invisible(total)
}
PT_LWD <- .15 * ggplot2::.stroke / 2                # outline of the B-D points, as a grid line width (keys)

# ---- page rows ----------------------------------------------------------------------------------------------------------
# Positions are the published page's exact values; where a short decimal is not exactly that value, the sum that gives it
# is written instead.
A_Y <- .46 - .16 - 147 / 600; A_H <- 2.935          # A: top and height of its plot (the column titles level with the letter)
Y2 <- 2.86; C2 <- 3.22; H2 <- 1.48                  # B-D row: its top; B's plot top and the plots' height
GROW_CD <- 176 / 600                                # C and D reach this much higher, beside their letters
E_TOP <- 4.88 - (.16 - .24) - .03                   # E: top of its slot
FG_TOP <- E_TOP + 3.26 + .10                        # F row: plot top
FG_H <- H - FG_TOP - .06 - .16                      # F plot height: one key row and a margin below
ROW_BCD_TOP <- round(Y2 * DPI) / DPI                # the panel images divide A from the B-D row here (whole pixels)

# 5A's color ramp: magma as in the RNA-level pages of Figures S10 and S11, same ends and the same stretch toward
# the low end, 100 steps
A_PAL <- scales::gradient_n_pal(viridisLite::magma(256, begin = .08, end = .96), seq(0, 1, length.out = 256)^1.5, "Lab")(seq(0, 1, length.out = 100))

# ---- A: regulator RNA by cell type (heatmap) and SC-SIC minus SC-beta per SC-islet study (strip below it) -------------
plot_A <- function() {
  # the 34 regulators in their display groups (the heatmap's column blocks), in drawing order
  A_GROUPS <- read.delim(file.path(REPO, "workflow", "scripts", "figure5", "plot", "regulator_groups.tsv"),
                         comment.char = "#", check.names = FALSE, fileEncoding = "UTF-8")
  a <- read_table("regulator_rna/f5_5a_regulator_rna.tsv")
  TFS <- A_GROUPS$tf
  lev <- a[a$kind == "level", ]; cnt <- a[a$kind == "sc_sic_above_sc_beta_studies", ]
  NSC <- unique(lev$studies[lev$column == "SC-SIC"]); NAD <- unique(lev$studies[lev$column == "adult islet beta"])
  lev <- rbind(lev, read_table("regulator_rna/f5_5a_regulator_rna__hio_mccauley_rows.tsv")[, names(lev)])   # HIO graft EC row
  # rows top to bottom: table column, drawn label, row group
  ROWA <- data.frame(column = c("Hickey 2023 gut EC", "Elmentaite 2021 gut EC", "Elmentaite 2020 gut EC", "Egozi 2021 gut EC", "McCauley 2023 HIO graft EC",
                                "Singh 2024 induced EC", "SC-SIC", "SC-beta", "SC-alpha", "adult islet beta", "adult islet alpha"),
                     label = c("Hickey 2023 (nuclei)", "Elmentaite 2021", "Elmentaite 2020", "Egozi 2021", "McCauley 2023", "Singh 2024",
                               "SC-SIC", "SC-beta", "SC-alpha", "Beta", "Alpha"),
                     rg = c(rep("gut", 4), "hio", "ind", rep("sc", 3), rep("ad", 2)))
  RGL <- c(gut = "Native\nintestinal EC", hio = "HIO graft EC", ind = "Induced EC", sc = sprintf("SC-islet\n(%d studies)", NSC), ad = sprintf("Adult islet\n(%d studies)", NAD))
  SETL <- c(`EC / SC-SIC association` = "EC / SC-SIC\nassociation", `Endocrine association` = "Endocrine\nassociation",
            `Islet association` = "Islet\nassociation", `Intestinal association` = "Intestinal\nassociation")
  set_of <- function(tf) factor(SETL[A_GROUPS$group[match(tf, TFS)]], SETL)
  AL <- data.frame(tf = lev$tf, row = ROWA$label[match(lev$column, ROWA$column)], rg = ROWA$rg[match(lev$column, ROWA$column)], value = lev$value)
  AL$x <- factor(AL$tf, TFS); AL$set <- set_of(AL$tf)
  AL$rg <- factor(RGL[AL$rg], RGL); AL$y <- factor(AL$row, rev(ROWA$label))
  pAh <- ggplot(AL, aes(x, y)) +
    geom_tile(aes(fill = value), color = "white", linewidth = .3) +
    scale_fill_gradientn(colours = A_PAL, limits = range(lev$value)) +
    facet_grid(rg ~ set, scales = "free", space = "free", switch = "y") +
    scale_x_discrete(expand = c(0, 0)) + scale_y_discrete(expand = c(0, 0)) + labs(x = NULL, y = NULL) + th(FALSE) + NOXL +
    theme(axis.text.x = element_blank(),
          strip.placement = "outside", strip.text.y.left = element_text(angle = 0, hjust = 1, vjust = .5, face = "bold", size = 7, lineheight = .9, margin = margin(0, 4, 0, 0)),
          strip.text.x = element_text(size = 7, face = "bold", lineheight = .9, vjust = 0, margin = margin(0, 0, 2, 0)), strip.clip = "off",
          panel.spacing.x = unit(3, "pt"), panel.spacing.y = unit(2.5, "pt"), axis.text.y = element_text(size = 7, margin = margin(0, 1, 0, 0)),
          plot.margin = margin(2, 4, 0, 2), plot.background = element_blank())
  # Strip: SC-SIC minus SC-beta in each SC-islet study (dots; samples with at least 20 cells averaged within a study) and
  # their mean (bar), which is the heatmap's SC-SIC row minus its SC-beta row.
  sx <- read_table("regulator_rna/tf_expression_by_sample.tsv")
  sx <- sx[sx$cohort == "sc" & as.logical(sx$min20) & sx$tf %in% TFS, ]
  sm <- aggregate(mean_log2cpm ~ family + population + tf, sx, mean)
  AS <- merge(sm[sm$population == "SC_SIC", c("family", "tf", "mean_log2cpm")], sm[sm$population == "SC_BETA", c("family", "tf", "mean_log2cpm")],
              by = c("family", "tf"), suffixes = c("_sic", "_beta"))
  AS$d <- AS$mean_log2cpm_sic - AS$mean_log2cpm_beta
  AS$x <- factor(AS$tf, TFS); AS$set <- set_of(AS$tf)
  AM <- aggregate(d ~ tf + x + set, AS, mean)
  # the per-study differences give back the table's study counts (SC-SIC above SC-beta) and the heatmap's rows
  chk_n <- aggregate(d ~ tf, AS, function(v) c(n = length(v), up = sum(v > 0)))
  stopifnot(chk_n$d[, "n"] == NSC, chk_n$d[match(cnt$tf, chk_n$tf), "up"] == cnt$value)
  lv <- function(col) lev$value[lev$column == col][match(AM$tf, lev$tf[lev$column == col])]
  stopifnot(max(abs(AM$d - (lv("SC-SIC") - lv("SC-beta")))) <= 1e-9)
  AS$lab <- AM$lab <- " "    # empty row-group strip: keeps the axis beside the strip; its label is drawn on the page
  # dot color: which cell type has more RNA in that study; an exact 0 is a study where neither has any
  AS$side <- ifelse(AS$d > 0, "sic", ifelse(AS$d < 0, "beta", "none"))
  stopifnot(!any(AS$side == "none" & AS$mean_log2cpm_sic != 0))
  SETS0 <- split(A_GROUPS$tf, factor(A_GROUPS$group, unique(A_GROUPS$group)))
  pAs <- ggplot(AS, aes(x, d)) +
    geom_rect(data = data.frame(set = AS$set[match(TFS[c(TRUE, FALSE)], AS$tf)], k = ave(seq_along(TFS), rep(seq_along(SETS0), lengths(SETS0)), FUN = seq_along)[c(TRUE, FALSE)]),
              aes(xmin = k - .5, xmax = k + .5, ymin = -Inf, ymax = Inf), fill = "grey96", inherit.aes = FALSE) +   # alternating gene bands, to line names up with the heatmap
    zero_h +
    geom_crossbar(data = AM, aes(ymin = d, ymax = d), width = .62, linewidth = .2, colour = "black") +   # mean bar, behind the dots
    ggbeeswarm::geom_quasirandom(aes(colour = side), width = .22, size = .7, stroke = 0) +
    scale_colour_manual(values = c(sic = POP_COL[["SC-SIC"]], beta = POP_COL[["SC-beta"]], none = "grey60")) +
    facet_grid(lab ~ set, scales = "free_x", space = "free_x", switch = "y") +
    scale_x_discrete(expand = expansion(add = .5)) + scale_y_continuous(breaks = c(-3, 0, 3)) +
    labs(x = NULL, y = NULL) + th(FALSE, TRUE) + NOXL +
    theme(axis.text.x = element_text(angle = 90, hjust = 1, vjust = .5, face = "italic", size = 7), strip.text.x = element_blank(),
          strip.placement = "outside", strip.text.y.left = element_text(angle = 0, hjust = 1, vjust = .5, face = "bold", size = 7, lineheight = .9, margin = margin(0, 4, 0, 0)),
          axis.text.y = element_text(size = 7, hjust = 1, margin = margin(0, 1, 0, 0)),
          panel.spacing.x = unit(3, "pt"), plot.margin = margin(3, 4, 2, 2), plot.background = element_blank())
  # heatmap over the strip, which is .74 in high
  p <- (pAh / pAs + plot_layout(heights = unit(c(1, .45 + .29), c("null", "in")))) & NOBG
  list(plot = p + plot_annotation(theme = NOBG), lims = range(lev$value))
}
draw_A_plot <- function(p) draw_at(patchworkGrob(p), .04, A_Y, 7.40, A_H)

# ---- B, C, D: promoter chromatin of the gene groups ----------------------------------------------------------------------
# B and C are drawn at the same inches per x unit (one day of B as wide as one unit of C's culture axis) and D is 1.95 in
# wide; the three boxes are 0.10 in apart, from 0.148 in to 7.46 in. Their fixed widths (axis, labels, margins and facet
# spacing) are measured on the open device, so the plot areas share one scale. Returns the three aligned plots, their
# boxes (x and w, inches) and the lines where the panel images divide B from C and C from D: halfway across the gaps,
# on whole pixels.
row_BCD <- function() {
  # B: Zhu 2023 days 21, 32 and 39 within SC-SIC, SC-beta and SC-alpha, each group against the typical gene of the same
  # cell type, drawn as C: each group's median gene (line and dots) and the middle half of its genes (light bar), on C's
  # fold axis, deep enough that no bar is cut
  zb <- read_table("c5_relative_zhu_hickey_singh.tsv"); zb <- zb[zb$study == ZHU, ]
  ZX <- c(`21` = 1, `32` = 2.15, `39` = 3.3); FLOOR <- -3.2
  stopifnot(nrow(zb) == 27L, min(zb$q25) > FLOOR)
  zb$P <- factor(zb$population, c("SC-SIC", "SC-beta", "SC-alpha")); zb$x <- unname(ZX[as.character(zb$day)] + c(G1 = -.2, G2 = 0, G3 = .2)[zb$g])
  pB <- ggplot(zb, aes(x, rel, color = g)) + zero_h +
    geom_rect(aes(xmin = x - .09, xmax = x + .09, ymin = q25, ymax = q75), fill = unname(GCOL_LIGHT[zb$g]), colour = NA) +
    geom_line(aes(group = g), linewidth = .45) +
    geom_point(aes(fill = g), shape = 21, colour = "black", stroke = .15, size = 1.3) +
    facet_wrap(~P, nrow = 1) + scale_color_manual(values = GCOL) + scale_fill_manual(values = GCOL) +
    scale_x_continuous(breaks = unname(ZX), labels = names(ZX), limits = c(.62, 3.68)) +
    scale_y_continuous(breaks = log2(c(.125, .25, .5, 1, 2)), labels = c("×0.125", "×0.25", "×0.5", "×1", "×2")) + coord_cartesian(ylim = c(FLOOR, 1.75)) +
    labs(x = quote("Day of differentiation (" * bold("Zhu 2023") * ")"), y = "Openness relative to the\ntypical gene promoter of\nthe same cell type (fold)") +
    th(FALSE, TRUE) +
    theme(panel.spacing = unit(5, "pt"), axis.text.x = element_text(size = 7), axis.title.x = element_text(size = 7, margin = margin(t = 3.5, b = 9.5)))

  # C: culture weeks 2-4, months 6 and 12 and the median of the three grafts (Augsornworawat 2023), each population's
  # groups against its typical gene: each group's median gene (line and dots; the graft median a hollow diamond) and,
  # behind it, the middle half of the group's genes (a bar in a light shade of the group's color, 25th to 75th
  # percentile). Single grafts are in the tables but not drawn. The axis is log2-spaced and labeled in fold, as deep as
  # B's; C shares B's y title.
  POPL <- c("SC-SIC", "SC-beta", "SC-alpha"); TXP <- c("txp_1", "txp_2", "txp_3")
  ce <- read_table("c5_culture_extended.tsv"); cw <- read_table("c5_culture_within_pops_wk2to4.tsv")
  # weeks 2-4 are the weeks 2-4 table's values
  chk_w <- merge(ce[ce$sample %in% cw$sample, c("population", "sample", "g", "rel")], cw[, c("population", "sample", "g", "rel")], by = c("population", "sample", "g"))
  stopifnot(nrow(chk_w) == nrow(cw), max(abs(chk_w$rel.x - chk_w$rel.y)) == 0)
  cgx <- read_table("c5_culture_extended_per_gene.tsv")
  # weeks and months: the per-gene medians give back the drawn values; graft median: the median of the three grafts'
  # values, and each gene's graft-median value is the median of its three graft values
  cs <- ce[grepl("^s6_", ce$sample), ]
  chk_c <- merge(aggregate(rel ~ population + sample + g, cgx[grepl("^s6_", cgx$sample), ], median), cs[, c("population", "sample", "g", "rel")], by = c("population", "sample", "g"))
  stopifnot(nrow(chk_c) == nrow(cs), max(abs(chk_c$rel.x - chk_c$rel.y)) <= 1e-9)
  chk_g <- merge(aggregate(rel ~ population + g, ce[ce$sample %in% TXP, ], median), ce[ce$sample == "graft_median", c("population", "g", "rel")], by = c("population", "g"))
  stopifnot(nrow(chk_g) == 9L, max(abs(chk_g$rel.x - chk_g$rel.y)) <= 1e-12)
  chk_gg <- merge(aggregate(rel ~ population + g + gene, cgx[cgx$sample %in% TXP, ], median), cgx[cgx$sample == "graft_median", c("population", "g", "gene", "rel")], by = c("population", "g", "gene"))
  stopifnot(nrow(chk_gg) == sum(cgx$sample == "graft_median"), max(abs(chk_gg$rel.x - chk_gg$rel.y)) <= 1e-12)
  stopifnot(all(as.logical(ce$drawn) == !(ce$sample %in% TXP)))
  CX <- c(s6_week2 = 1, s6_week3 = 2.15, s6_week4 = 3.3, s6_month6 = 4.6, s6_month12 = 5.85, graft_median = 6.95); GDX <- c(G1 = -.2, G2 = 0, G3 = .2)
  cu <- ce[ce$sample %in% names(CX), ]; cu$P <- factor(cu$population, POPL); cu$x <- unname(CX[cu$sample] + GDX[cu$g])
  cg <- cgx[cgx$sample %in% names(CX), ]; cg$P <- factor(cg$population, POPL); cg$x <- unname(CX[cg$sample] + GDX[cg$g])
  cq <- do.call(rbind, lapply(split(cg, interaction(cg$g, cg$sample, cg$P, drop = TRUE)), function(d)
    data.frame(g = d$g[1], P = d$P[1], x = d$x[1], lo = unname(quantile(d$rel, .25)), hi = unname(quantile(d$rel, .75)))))
  pC <- ggplot(cu, aes(x, rel, color = g)) + zero_h +
    geom_vline(xintercept = c(3.95, 6.43), color = "#BBBBBB", linewidth = .2, linetype = "22") +
    geom_rect(data = cq, aes(xmin = x - .09, xmax = x + .09, ymin = lo, ymax = hi), fill = unname(GCOL_LIGHT[cq$g]), colour = NA, inherit.aes = FALSE) +
    geom_line(data = cu[cu$sample != "graft_median", ], aes(group = g), linewidth = .45) +
    geom_point(data = cu[cu$sample != "graft_median", ], aes(fill = g), shape = 21, colour = "black", stroke = .15, size = 1.3) +
    geom_point(data = cu[cu$sample == "graft_median", ], shape = 23, fill = "white", stroke = .6, size = 1.6) +
    # two-tier x labels: numbers on the ticks, "weeks" and "months" under short brackets below each group, at a fixed
    # distance below the axis floor
    geom_segment(data = data.frame(x = c(CX[[1]], CX[[4]]) - .25, xend = c(CX[[3]], CX[[5]]) + .25), aes(x = x, xend = xend, y = FLOOR - .92 * 4.95 / 4.45, yend = FLOOR - .92 * 4.95 / 4.45), colour = "black", linewidth = .2, inherit.aes = FALSE) +
    geom_text(data = data.frame(x = c(mean(CX[1:3]) - .3, mean(CX[4:5]) + .3), lab = c("weeks", "months")), aes(x = x, y = FLOOR - 1 * 4.95 / 4.45, label = lab), vjust = 1, size = 7 / ggplot2::.pt, family = "Arial", colour = "black", inherit.aes = FALSE) +
    facet_wrap(~P, nrow = 1) + scale_color_manual(values = GCOL) + scale_fill_manual(values = GCOL) +
    scale_x_continuous(breaks = unname(CX[1:5]), labels = c("2", "3", "4", "6", "12"), limits = c(.62, 7.3)) +
    scale_y_continuous(breaks = log2(c(.125, .25, .5, 1, 2)), labels = c("×0.125", "×0.25", "×0.5", "×1", "×2")) + coord_cartesian(ylim = c(FLOOR, 1.75), clip = "off") +
    labs(x = quote("Stage 6 culture (" * bold("Augsornworawat 2023") * ")"), y = NULL) + th(FALSE, TRUE) +
    theme(panel.spacing = unit(5, "pt"), axis.text.x = element_text(size = 7), axis.title.x = element_text(size = 7, margin = margin(t = 10.5, b = 2.5)))

  # D: grafts, each group's change from week 2 beyond the typical gene's, within SC-SIC, SC-beta and SC-alpha, against the
  # middle 95% of 10,000 random gene
  # sets; promoters where the promoter and distal-element calls agree in every graft and group, else distal elements
  WIN <- if (all(as.character(read_table("c5_graft_window_agreement.tsv")$same) %in% c("TRUE", "True", "true"))) "promoter" else "distal"
  gr <- read_table("c5_graft_relative_by_population.tsv"); gr <- gr[gr$window == WIN, ]
  gr$y <- pct(gr$group_rel); gr$lo <- pct(gr$null_q025_rel); gr$hi <- pct(gr$null_q975_rel)
  gr$u <- factor(c(txp_1 = "1", txp_2 = "2", txp_3 = "3")[gr$graft], c("1", "2", "3")); gr$win <- factor(gr$population, c("SC-SIC", "SC-beta", "SC-alpha"))
  pD <- ggplot(gr, aes(u, y, group = g)) + zero_h +
    geom_linerange(aes(ymin = lo, ymax = hi), position = position_dodge(.8), linewidth = .85, color = BAND) +
    geom_point(aes(fill = g), shape = 21, colour = "black", stroke = .15, size = 1.2, position = position_dodge(.8)) + facet_wrap(~win, nrow = 1) +
    scale_fill_manual(values = GCOL) + labs(x = "Graft (week 2 → graft)", y = "% more opening than the\ntypical gene's distal elements") +
    th(FALSE, TRUE) + NOXL + theme(panel.spacing = unit(5, "pt"), axis.text.x = element_text(size = 7))

  fixed_in <- function(p) { g <- ggplotGrob(p); nl <- unitType(g$widths) == "null"; convertWidth(sum(g$widths[!nl]), "in", valueOnly = TRUE) }
  K <- diff(c(.62, 7.3)) / diff(c(.62, 3.68))        # C's x range over B's
  fB <- fixed_in(pB); fC <- fixed_in(pC)
  X0 <- .148; X1 <- 7.46; WD <- 1.95; GAP <- .10
  WB <- (X1 - X0 - WD - 2 * GAP - fC + K * fB) / (1 + K); WC <- fC + K * (WB - fB)
  box <- list(B = c(x = X0, w = WB), C = c(x = X0 + WB + GAP, w = WC), D = c(x = X1 - WD, w = WD))
  cut <- c(BC = round((box$B[["x"]] + box$B[["w"]] + box$C[["x"]]) / 2 * DPI) / DPI,
           CD = round((box$C[["x"]] + box$C[["w"]] + box$D[["x"]]) / 2 * DPI) / DPI)
  c(setNames(align_h(list(pB + NOBG, pC + NOBG, pD + NOBG)), c("B", "C", "D")), list(box = box, cut = cut))
}
# Draw plot k ("B", "C" or "D") of row_BCD() in its box; each reaches up beside its letter.
draw_BCD <- function(g, k) draw_at(g[[k]], g$box[[k]][["x"]], C2 - GROW_CD, g$box[[k]][["w"]], H2 + GROW_CD)

# ---- F: the CDX2 knockout within each cell type ----------------------------------------------------------------------------
# F is the full page wide and uses the supplementary panels' theme at 7 pt (theme_fg), with its key drawn under it
# (no ggplot legends).
GRID_COL <- "#E4E4E4"; REF_COL <- "#777777"         # gridlines; zero lines
theme_fg <- function(grid = c("x", "y", "both", "none")) {
  grid <- match.arg(grid)
  gl <- element_line(color = GRID_COL, linewidth = .2)
  theme_minimal(base_family = "Arial", base_size = 7) +
    theme(text = element_text(color = "black", size = 7),
          axis.text = element_text(color = "black", size = 7),
          axis.title = element_text(color = "black", size = 7),
          strip.text = element_text(color = "black", size = 8, face = "bold", hjust = 0),
          legend.position = "none",
          panel.grid.minor = element_blank(),
          panel.grid.major.x = if (grid %in% c("x", "both")) gl else element_blank(),
          panel.grid.major.y = if (grid %in% c("y", "both")) gl else element_blank(),
          axis.line.x = element_line(color = "black", linewidth = .35),
          axis.ticks.x = element_line(color = "black", linewidth = .3),
          axis.ticks.length = unit(1.3, "pt"),
          panel.spacing = unit(.16, "in"),
          plot.margin = margin(2, 6, 2, 2),
          plot.background = element_blank())
}
FG_BAR <- "grey82"                                  # random-range bars
panel_F <- function() {
  # the knockout's effect within each cell type (no SC-beta subtraction), knockout minus wild type,
  # on G2's four EC markers, the other G2 genes and G3, in SC-SIC, SC-beta and SC-alpha,
  # pseudobulk (log2) and per-cell detection (percentage points of group genes), for the nuclei of the study's GEO deposit;
  # grey bars: middle 95% of 500 matched random gene sets; dots in each cell type's color, filled outside that range,
  # a ring inside it (both the same outer size).
  E <- read_table("cdx2_knockout/group_test/group_effects.tsv")
  SETD <- c("G2 serotonin core" = "G2: four\nEC markers", "G2 other" = "G2: other\ngenes", "G3" = "G3")
  E <- E[E$variant == "all transferred" & E$set %in% names(SETD) & E$population %in% c("SC_EC", "SC_beta", "SC_alpha"), ]
  stopifnot(nrow(E) == 18)
  k <- ifelse(grepl("per-cell", E$readout), 100, 1)                     # per-cell detection as percentage points
  E$v <- E$ko_minus_wt * k; E$lo <- E$random_q025 * k; E$hi <- E$random_q975 * k
  RD <- c("Pseudobulk\n(log2)", "Per-cell detection\n(points of group genes)")
  E$rd <- factor(ifelse(grepl("per-cell", E$readout), RD[2], RD[1]), RD)
  E$pop <- c(SC_EC = "SC-SIC", SC_beta = "SC-beta", SC_alpha = "SC-alpha")[E$population]
  # one facet row per gene set (bold label at left), SC-SIC, SC-beta, SC-alpha inside it
  E$grp <- factor(SETD[E$set], SETD)
  E$y <- c("SC-SIC" = 1, "SC-beta" = 2, "SC-alpha" = 3)[E$pop]
  E$out <- factor(ifelse(E$v < E$lo | E$v > E$hi, "Outside random range", "Inside"), c("Outside random range", "Inside"))
  pG <- ggplot(E) +
    geom_vline(xintercept = 0, colour = REF_COL, linewidth = .25) +
    geom_linerange(aes(y = y, xmin = lo, xmax = hi), colour = FG_BAR, linewidth = 1.6) +
    geom_point(data = E[E$out == "Inside", ], aes(y = y, x = v), shape = 21, fill = "white", colour = POP_COL[E$pop[E$out == "Inside"]], size = 1.44, stroke = .8) +
    geom_point(data = E[E$out != "Inside", ], aes(y = y, x = v), shape = 21, fill = POP_COL[E$pop[E$out != "Inside"]], colour = POP_COL[E$pop[E$out != "Inside"]], size = 2.1, stroke = .3) +
    facet_grid(grp ~ rd, scales = "free_x", switch = "y") +
    scale_y_reverse(breaks = 1:3, labels = c("SC-SIC", "SC-beta", "SC-alpha"), limits = c(3.5, .5), expand = c(0, 0)) +
    scale_x_continuous(expand = expansion(mult = .09)) +
    labs(x = expression(italic("CDX2")~"knockout minus wild type, within each cell type"), y = NULL) +
    theme_fg("x") +
    theme(strip.text.x = element_text(hjust = .5, size = 7, face = "bold"), strip.placement = "outside",
          strip.text.y.left = element_text(angle = 0, hjust = 1, size = 7, face = "bold"), panel.spacing.y = unit(.05, "in"),
          panel.border = element_rect(fill = NA, colour = "grey80", linewidth = .3),
          axis.text.y = element_text(colour = unname(POP_COL), face = "bold"),
          panel.grid.major.x = element_blank())

  ggplotGrob(pG)
}
FG_PT_LWD <- .3 * ggplot2::.stroke / 2              # outline of the F points, as a grid line width (keys)
