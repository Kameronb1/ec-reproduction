# Drawing helpers shared by the Figure 5 supplementary panels (pages S17-S19): text sizes and colors, the panel theme,
# placing a panel on its page, and the identity-factor rows of pages S18 and S19. Sourced by
# workflow/scripts/common/draw_panel.R, which defines the page size W x H (inches).
#
# Every panel spans 7.3 in from x = 0.20 in. Its content top (in whole pixels from the page top) puts its first ink
# 0.10 in under its letter's baseline, the panels following each other 0.12 in apart from the page top.

FONT <- "Arial"; PT <- 7.5; TXT <- PT / .pt         # body text (pt); TXT: ggplot text size of 7.5 pt
GRID_COL <- "#E4E4E4"; REF_COL <- "#777777"         # gridlines; zero lines
SC_SIC <- "#AA3377"; SC_BETA <- "#4477AA"
STROKE_THIN <- .12                                  # dot outlines: very thin black
PANEL_X <- .20; PANEL_W <- 7.3

# Panel theme: 7.5 pt black text, light major gridlines on the chosen axis, black x axis line, no background
theme_panel <- function(grid = c("x", "y", "both", "none"), legend = "none") {
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
          panel.grid.major.x = if (grid %in% c("x", "both")) gl else element_blank(),
          panel.grid.major.y = if (grid %in% c("y", "both")) gl else element_blank(),
          axis.line.x = element_line(color = "black", linewidth = .35),
          axis.ticks.x = element_line(color = "black", linewidth = .3),
          axis.ticks.length = unit(1.3, "pt"),
          panel.spacing = unit(.16, "in"),
          plot.margin = margin(2, 6, 2, 2),
          plot.background = element_blank())
}

# Widest of some strings, in inches, on the open device
str_in <- function(s, pt = PT, bold = FALSE) {
  pushViewport(viewport(gp = gpar(fontfamily = FONT, fontsize = pt, fontface = if (bold) "bold" else "plain")))
  on.exit(popViewport())
  max(vapply(s, function(z) convertWidth(stringWidth(z), "in", valueOnly = TRUE), 0))
}

# Draw a panel's content (a ggplot or patchwork, or a function drawing with grid) into its place: top_px pixels below
# the page top, px pixels high.
place <- function(content, top_px, px) {
  pushViewport(viewport(x = unit(PANEL_X, "in"), y = unit(1, "npc") - unit(top_px / DPI, "in"), width = unit(PANEL_W, "in"),
                        height = unit(px / DPI, "in"), just = c("left", "top")))
  if (is.function(content)) content() else print(content, newpage = FALSE)
  popViewport()
}

# ---- page S17: four rows of panels side by side ---------------------------------------------------------------------------
# A | B | C over their shared key; D | E; F | G; H | I. Each row's content is px pixels high from top pixels below the page
# top; its panels sit WG apart across PANEL_W at the widths ws (NA: the rest of the row). Positions inside a panel are
# inches from its bottom-left corner.
WG <- .15
S17_ROWS <- list(ABC = list(top = 171, px = 1200, ws = PANEL_W), DE = list(top = 1574, px = 1140, ws = c(2.30, NA)),
                 FG = list(top = 2892, px = 1260, ws = c(2.38, NA)), HI = list(top = 4369, px = 1560, ws = c(3.95, NA)))
in_vp <- function(x, y, w, h, f) {
  pushViewport(viewport(x = unit(x, "in"), y = unit(y, "in"), width = unit(w, "in"), height = unit(h, "in"), just = c("left", "bottom")))
  f(); popViewport()
}
# Draw f(w, vh) into slot k of a row of page S17, w wide and vh high (inches).
in_row <- function(row, k, f) {
  r <- S17_ROWS[[row]]; ws <- r$ws; ws[is.na(ws)] <- PANEL_W - sum(ws, na.rm = TRUE) - WG * (length(ws) - 1)
  xs <- c(0, cumsum(ws + WG))[seq_along(ws)]
  place(function() { vh <- convertHeight(unit(1, "npc"), "in", valueOnly = TRUE); in_vp(xs[k], 0, ws[k], vh, function() f(ws[k], vh)) }, r$top, r$px)
}
GCOL <- c(G1 = "#43A5B6", G2 = "#A89A2E", G3 = "#B4553F"); BANDC <- "#D9D9D9"   # gene-group colors as in Figure 5
GCOL_LIGHT <- vapply(GCOL, function(c) grDevices::rgb(t(.4 * grDevices::col2rgb(c) + .6 * 255), maxColorValue = 255), "")   # their 40% shades
# One row of key items centered on cx at height y (inches in the current viewport): gp (a point), pb (a gene group as A-C
# draw it, its point on its light bar), bar (the grey random-set bar), bars3 (three light group bars) or dash, each with
# its label.
ky2 <- function(items, cx, y, gap = .2) {
  wd <- function(it) str_in(it$label) + if (it$type %in% c("bars3", "dash")) .26 else .14
  ww <- vapply(items, wd, 0); x <- cx - (sum(ww) + gap * (length(ww) - 1)) / 2
  for (i in seq_along(items)) {
    it <- items[[i]]
    if (it$type == "pb") grid.rect(unit(x + .05, "in"), unit(y, "in"), unit(.05, "in"), unit(.13, "in"), gp = gpar(fill = it$bar, col = NA))
    if (it$type %in% c("gp", "pb")) grid.points(unit(x + .05, "in"), unit(y, "in"), pch = if (is.null(it$pch)) 21 else it$pch,
                                                size = unit(if (is.null(it$size)) 1.6 else it$size, "mm"), gp = gpar(col = "black", fill = it$fill, lwd = .4))
    if (it$type == "bar") grid.rect(unit(x + .05, "in"), unit(y, "in"), unit(.85 * .pt / 72, "in"), unit(.13, "in"), gp = gpar(fill = BANDC, col = NA))
    if (it$type == "bars3") for (j in 1:3)
      grid.rect(unit(x + .03 + (j - 1) * .07, "in"), unit(y, "in"), unit(.05, "in"), unit(.12, "in"), gp = gpar(fill = GCOL_LIGHT[[j]], col = NA))
    if (it$type == "dash") grid.lines(unit(c(x + .03, x + .2), "in"), unit(c(y, y), "in"), gp = gpar(col = "black", lwd = 1))
    grid.text(it$label, unit(x + if (it$type %in% c("bars3", "dash")) .26 else .14, "in"), unit(y, "in"), just = c("left", "centre"),
              gp = gpar(fontfamily = FONT, fontsize = PT, col = "black"))
    x <- x + ww[i] + gap
  }
  invisible(sum(ww) + gap * (length(ww) - 1))
}
# Figure 5's theme, colors and helpers (../plot/helpers.R) in their own environment, for A-C and F, which use Figure 5's
# text and mark sizes.
fig5_env <- function() {
  m <- new.env()
  sys.source(file.path(REPO, "workflow", "scripts", "figure5", "plot", "helpers.R"), envir = m)
  m
}

# ---- page S17, A-C: promoter accessibility of the gene groups --------------------------------------------------------------
# A: each group's median gene promoter in SC-SIC minus SC-beta and minus SC-alpha, per sample (Augsornworawat 2023
# circles, Zhu 2023 triangles). B and C: Figure 5C's measure (each group against the typical gene of the same cell type),
# drawn as Figure 5B draws Zhu 2023: each group's median gene (line and dots) on a light bar over the middle half of its
# genes, on a fold axis deep enough that no bar is cut. B: Hickey 2023 native EC and the other enteroendocrine cells by
# region, donors pooled; C: Singh 2024 induced EC, D, X and T4 cells (primary study labels) on days 5 and 6 of induction.
# The three are built together so that their plot areas line up; A is TOPW[1] and B TOPW[2] wide, C takes the rest of the
# row, and the key band TOPKEY high is under all three.
TOPW <- c(1.90, 2.25); TOPKEY <- .36
top_row <- function() {
  m <- fig5_env()
  m$FLOOR <- -3.65; m$FOLD_BR <- log2(c(.125, .25, .5, 1, 2)); m$FOLD_LAB <- c("×0.125", "×0.25", "×0.5", "×1", "×2")
  pA <- with(m, {
    pr <- read_table("c5_promoters_by_sibling.tsv"); pr <- pr[pr$sibling %in% c("Beta", "Alpha"), ]; pr$y <- pct(pr$median_diff)
    pr$sib <- factor(c(Beta = "vs SC-beta", Alpha = "vs SC-alpha")[pr$sibling], c("vs SC-beta", "vs SC-alpha")); pr$sn <- factor(pr$study, c(AUG, ZHU))
    pr$slot <- factor(paste(pr$sn, pr$sample), unique(paste(pr$sn, pr$sample)[order(pr$sn, pr$sample)]))
    ggplot(pr, aes(g, y, fill = g, shape = sn, group = slot)) + zero_h +
      geom_point(position = position_dodge(width = .84), size = 1.35, stroke = .15, colour = "black") +
      facet_wrap(~sib, nrow = 1) + scale_fill_manual(values = GCOL) + scale_shape_manual(values = c(21, 24)) +
      labs(x = quote(bold("Augsornworawat 2023") * " and " * bold("Zhu 2023")), y = "% more open in SC-SIC\n(group's median gene promoter)") +
      th(FALSE, TRUE) + NOXL + theme(panel.spacing = unit(5, "pt"), plot.background = element_blank())
  })
  pB <- with(m, {
    hk <- read_table("c5_relative_zhu_hickey_singh.tsv"); hk <- hk[hk$study == "Hickey 2023", ]; RX <- c(`small intestine` = 1, colon = 2.15)
    hk$P <- factor(c(EC = "Native EC", `other enteroendocrine` = "Other entero-\nendocrine")[hk$population], c("Native EC", "Other entero-\nendocrine"))
    hk$x <- unname(RX[hk$region] + c(G1 = -.2, G2 = 0, G3 = .2)[hk$g])
    stopifnot(nrow(hk) == 12, !anyNA(hk$P), !anyNA(hk$x), min(hk$q25) > FLOOR)
    ggplot(hk, aes(x, rel, color = g)) + zero_h +
      geom_rect(aes(xmin = x - .09, xmax = x + .09, ymin = q25, ymax = q75), fill = unname(GCOL_LIGHT[hk$g]), colour = NA) +
      geom_line(aes(group = g), linewidth = .45) + geom_point(aes(fill = g), shape = 21, colour = "black", stroke = .15, size = 1.3) +
      facet_wrap(~P, nrow = 1) + scale_color_manual(values = GCOL) + scale_fill_manual(values = GCOL) +
      scale_x_continuous(breaks = unname(RX), labels = c("Small\nintestine", "Colon"), limits = c(.5, 2.65)) +
      scale_y_continuous(breaks = FOLD_BR, labels = FOLD_LAB) + coord_cartesian(ylim = c(FLOOR, 1.75)) +
      labs(x = quote("Region (" * bold("Hickey 2023") * ")"), y = "Openness relative to the\ntypical gene promoter of\nthe same cell type (fold)") +
      th(FALSE, TRUE) + theme(panel.spacing = unit(5, "pt"), axis.text.x = element_text(size = 7, lineheight = .9), plot.background = element_blank())
  })
  pC <- with(m, {
    sg <- read_table("c5_relative_zhu_hickey_singh.tsv"); sg <- sg[sg$study == "Singh 2024", ]; DX <- c(`5` = 1, `6` = 2.15)
    sg$P <- factor(c(`induced EC` = "Induced EC", `D cells` = "D cells", `X cells` = "X cells", `T4 cells` = "T4 cells")[sg$population], c("Induced EC", "D cells", "X cells", "T4 cells"))
    sg$x <- unname(DX[as.character(sg$day)] + c(G1 = -.2, G2 = 0, G3 = .2)[sg$g])
    stopifnot(nrow(sg) == 24, !anyNA(sg$P), !anyNA(sg$x), min(sg$q25) > FLOOR)
    ggplot(sg, aes(x, rel, color = g)) + zero_h +
      geom_rect(aes(xmin = x - .09, xmax = x + .09, ymin = q25, ymax = q75), fill = unname(GCOL_LIGHT[sg$g]), colour = NA) +
      geom_line(aes(group = g), linewidth = .45) + geom_point(aes(fill = g), shape = 21, colour = "black", stroke = .15, size = 1.3) +
      facet_wrap(~P, nrow = 1) + scale_color_manual(values = GCOL) + scale_fill_manual(values = GCOL) +
      scale_x_continuous(breaks = unname(DX), labels = names(DX), limits = c(.5, 2.65)) +
      scale_y_continuous(breaks = FOLD_BR, labels = FOLD_LAB) + coord_cartesian(ylim = c(FLOOR, 1.75)) +
      labs(x = quote("Day of induction (" * bold("Singh 2024") * ")"), y = NULL) +
      th(FALSE, TRUE) + theme(panel.spacing = unit(5, "pt"), axis.text.x = element_text(size = 7), plot.background = element_blank())
  })
  ws <- c(TOPW, PANEL_W - sum(TOPW) - 2 * WG)
  list(g = m$align_h(list(pA, pB, pC)), x = c(0, cumsum(ws + WG))[1:3], w = ws)
}
# Draw plot k (1-3: A, B, C) of top_row() in its place above the key band.
draw_top <- function(t, k) in_row("ABC", 1, function(w, vh) in_vp(t$x[k], TOPKEY, t$w[k], vh - TOPKEY, function() grid.draw(t$g[[k]])))

# ---- page S17, F: the CDX2 knockout's cell-type shares -------------------------------------------------------------------
# Share of each transferred cell type, wild type -> knockout (dots: label transfer; dashes: Zhu 2023's published shares).
# Labels at the knockout end, pushed apart to at least 7.5 % units, with a leader where a label sits more than 2.5 % from
# its point. Drawn with Figure 5's helpers (fig5_env), at Figure 5's text and mark sizes.
shares_plot <- function() with(fig5_env(), {
  fa <- read_table("cdx2_knockout/figure5d_comparison.tsv"); fa <- fa[fa$arm == "all_stages", ]
  LABK <- c(SC_EC = "SC-SIC", SC_beta = "SC-beta", SC_alpha = "SC-alpha", SC_delta = "SC-delta", ENP3 = "ENP3")
  fa$lab <- factor(LABK[fa$label], LABK); fa$x <- ifelse(fa$genotype == "WT", 1, 2); fa$study <- ZHU
  wt <- fa[fa$genotype == "WT", ]; ko <- fa[fa$genotype == "CDX2-KO", ]; ko <- ko[match(wt$label, ko$label), ]
  ratio <- ko$ours[wt$label == "SC_EC"] / wt$ours[wt$label == "SC_EC"]
  KCOL <- c(`SC-SIC` = "#AA3377", `SC-beta` = "#4477AA", `SC-alpha` = "#117733", `SC-delta` = "#AE8C44", ENP3 = "grey60")
  KL <- data.frame(lab = as.character(ko$lab), y = 100 * ko$ours); KL <- KL[order(KL$y), ]
  KL$ly <- KL$y; for (i in seq_len(nrow(KL))[-1]) KL$ly[i] <- max(KL$y[i], KL$ly[i - 1] + 7.5)
  KL$text <- ifelse(KL$lab == "SC-SIC", sprintf("SC-SIC ×%.2f", ratio), KL$lab); KLL <- KL[(KL$ly - KL$y) > 2.5, ]
  ggplot(fa, aes(x, 100 * ours, color = lab, group = lab)) +
    annotate("segment", x = .8, xend = 2.2, y = seq(0, 60, 20), yend = seq(0, 60, 20), colour = "#E4E4E4", linewidth = .2) +
    geom_line(aes(linewidth = lab == "SC-SIC")) + geom_point(size = 1.1) +
    geom_segment(aes(x = x + ifelse(x == 1, -.12, .12) - .045, xend = x + ifelse(x == 1, -.12, .12) + .045, y = 100 * zhu_fig5d, yend = 100 * zhu_fig5d),
                 linewidth = .35, lineend = "butt") +
    geom_segment(data = KLL, aes(x = 2.20, xend = 2.29, y = y, yend = ly), inherit.aes = FALSE, linewidth = .2, colour = "black") +
    geom_text(data = KL, aes(x = 2.32, y = ly, label = text, fontface = ifelse(lab == "SC-SIC", "bold", "plain")), inherit.aes = FALSE,
              hjust = 0, size = TM, family = "Arial", color = "black") +
    facet_wrap(~study, labeller = as_labeller(function(x) rep("", length(x)))) +
    scale_linewidth_manual(values = c(`TRUE` = 1.0, `FALSE` = .45)) + scale_color_manual(values = KCOL) +
    scale_x_continuous(breaks = 1:2, labels = c("Wild\ntype", "CDX2\nknockout"), limits = c(.8, 3.25), expand = c(0, 0)) +
    scale_y_continuous(limits = c(0, 66), breaks = seq(0, 60, 20)) + coord_cartesian(clip = "off") +
    labs(x = NULL, y = "% of nuclei\n(five types = 100%)") + th(FALSE, TRUE) + NOXL + theme(panel.grid.major.y = element_blank())
})

# ---- page S17, H and I: the CDX2 knockout and the targets of regulons -------------------------------------------------------
# H: each CDX2(+) regulon's targets, knockout minus wild type, against 500 expression-matched random gene sets in SC-SIC,
# SC-beta and SC-alpha (pseudobulk over per-cell detection); I: the reproducible scorable Figure 5E regulons in SC-SIC
# (pseudobulk), by direction. regulon_panels() returns an environment holding drawH(w) and drawI(w) (each w inches wide).
regulon_panels <- function() {
  e <- new.env()
  local({
    K_E <- read_table("cdx2_knockout/regulon_targets/set_effects.tsv"); K_C <- read_table("cdx2_knockout/regulon_targets/context_regulons.tsv")
    K_PB <- "pseudobulk log2 (equal depth)"; K_PC <- "per-cell detection (depth-stratified)"
    K_NOT <- "grey50"; K_BAR <- "#D6D6D6"; K_MED <- "grey25"; K_OTHER <- "grey72"
    K_PT <- 2.0; K_ST <- STROKE_THIN                                  # CDX2 dots: size and black outline stroke
    # grid key items laid out left to right, centered on cx (inches in the current viewport); types pt, bar, tick
    key_items <- function(items, cx, y) {
      ww <- vapply(items, function(it) str_in(it$label) + .17, 0); x <- cx - (sum(ww) + .16 * (length(ww) - 1)) / 2
      for (i in seq_along(items)) { it <- items[[i]]; px <- unit(x + .06, "in"); py <- unit(y, "in")
        sz <- if (is.null(it$size)) K_PT else it$size; st <- if (is.null(it$stroke)) K_ST else it$stroke
        if (it$type == "pt") grid.points(px, py, pch = 21, gp = gpar(col = "black", fill = it$fill, fontsize = sz * .pt + st * .stroke / 2, lwd = st * .stroke / 2))
        if (it$type == "bar") grid.rect(px, py, unit(.13, "in"), unit(3.2 * .pt / 72 * .75, "in"), gp = gpar(fill = K_BAR, col = NA))
        if (it$type == "tick") grid.segments(px, py - unit(.055, "in"), px, py + unit(.055, "in"), gp = gpar(col = it$col, lwd = it$lwd * .pt))
        grid.text(it$label, unit(x + .17, "in"), py, just = c("left", "centre"), gp = gpar(fontfamily = FONT, fontsize = PT, col = "black"))
        x <- x + ww[i] + .16 }
    }
    minus <- function(x) sub("^-", "−", format(x, trim = TRUE, drop0trailing = TRUE))

    # H: the four CDX2(+) regulons (network and targets found in the expression matrix), three cell types, two readouts; the number right of each
    # row is the fraction of random sets as low or lower (bold below 0.05)
    K_SETS <- data.frame(set = c("GSE167880__CDX2(+)", "GSE185036__CDX2(+)", "HICKEY_2023_DRYAD__CDX2(+)", "EGOZI_2021_ZENODO__CDX2(+)"),
                         study = c("Balboa 2022", "Du 2022", "Hickey 2023", "Egozi 2021"),
                         rep = c(TRUE, TRUE, FALSE, FALSE), y = c(4, 3, 1.75, .75))
    # row label: study and the regulon's targets found in the expression matrix (in_matrix, one value per set)
    K_N <- tapply(K_E$in_matrix[K_E$scope == "all targets"], K_E$set[K_E$scope == "all targets"], unique)
    stopifnot(all(lengths(K_N[K_SETS$set]) == 1))
    K_SETS$lab <- sprintf("%s (%d)", K_SETS$study, as.integer(unlist(K_N[K_SETS$set])))
    K_POP <- c(SC_EC = "SC-SIC", SC_beta = "SC-beta", SC_alpha = "SC-alpha")
    KA <- merge(K_E[K_E$scope == "all targets" & K_E$set %in% K_SETS$set & K_E$population %in% names(K_POP), ], K_SETS, by = "set")
    stopifnot(nrow(KA) == 24)
    k <- KA$readout == K_PC; V <- c("ko_minus_wt", "random_median", "random_q025", "random_q975")
    KA[k, V] <- 100 * KA[k, V]                                    # per-cell detection in percentage points
    KA$pop <- factor(K_POP[KA$population], K_POP)
    KA$share <- ifelse(KA$random_share_at_or_below < .01, "<0.01", sprintf("%.2f", KA$random_share_at_or_below))
    KA$hit <- KA$random_share_at_or_below < .05
    k_block <- function(ro, xlab, top) {
      d <- KA[KA$readout == ro, ]
      r <- range(c(d$ko_minus_wt, d$random_q025, d$random_q975)); r <- r + c(-.05, .05) * diff(r)
      xt <- r[2] + .07 * diff(r); xr <- r[2] + .62 * diff(r)
      br <- scales::breaks_pretty(3)(r); br <- br[br >= r[1] & br <= r[2]]
      ggplot(d) +
        geom_vline(xintercept = 0, colour = REF_COL, linewidth = .3) +
        geom_segment(aes(x = random_q025, xend = random_q975, y = y, yend = y), colour = K_BAR, linewidth = 3.2 * .75, lineend = "butt") +
        geom_segment(aes(x = random_median, xend = random_median, y = y - .3, yend = y + .3), colour = K_MED, linewidth = .45) +
        geom_point(aes(ko_minus_wt, y, fill = rep), shape = 21, colour = "black", stroke = K_ST, size = K_PT) +
        geom_text(aes(x = xt, y = y, label = share, fontface = ifelse(hit, "bold", "plain")), hjust = 0, size = TXT, family = FONT) +
        (if (top) geom_text(data = data.frame(pop = factor("SC-SIC", K_POP)), aes(x = r[1] - .2 * diff(r), y = 5.45),
                            label = "CDX2(+) regulon:\nnetwork (targets)", hjust = 1, vjust = 0, size = TXT, family = FONT, fontface = "italic", lineheight = .95)) +
        facet_grid(. ~ pop) +
        scale_y_continuous(breaks = K_SETS$y, labels = K_SETS$lab, limits = c(.25, if (top) 5.15 else 4.5), expand = c(0, 0), oob = scales::oob_keep) +
        scale_x_continuous(limits = c(r[1], xr), breaks = br, labels = minus, expand = c(0, 0), oob = scales::oob_keep) +
        scale_fill_manual(values = c(`TRUE` = SC_SIC, `FALSE` = K_NOT), guide = "none") +
        coord_cartesian(clip = "off") + labs(x = xlab, y = NULL) + theme_panel("x") +
        theme(strip.text.x = if (top) element_text(face = "bold", size = PT, hjust = .3, margin = margin(0, 0, 2, 0)) else element_blank(),
              axis.text.y = element_text(hjust = 1, size = PT), panel.spacing.x = unit(.1, "in"), plot.margin = margin(if (top) 9 else 1, 2, 2, 2),
              axis.title.x = element_text(margin = margin(2, 0, 0, 0)))
    }
    pH <- k_block(K_PB, "CDX2 knockout minus wild type, pseudobulk (log2)", TRUE) /
          k_block(K_PC, "CDX2 knockout minus wild type, per-cell detection (points)", FALSE) + plot_layout(heights = c(1.12, 1))
    H_KEY <- .48                                                  # two key rows and the numbers' line under H
    drawH <- function(w) {
      vh <- convertHeight(unit(1, "npc"), "in", valueOnly = TRUE)
      pushViewport(viewport(x = 0, y = unit(H_KEY, "in"), width = unit(w, "in"), height = unit(vh - H_KEY, "in"), just = c("left", "bottom")))
      print(pH, newpage = FALSE); popViewport()
      cx <- w / 2 + .3
      key_items(list(list(type = "pt", label = "Reproducible CDX2(+) regulon", fill = SC_SIC), list(type = "pt", label = "Not reproducible", fill = K_NOT)), cx, H_KEY - .1)
      key_items(list(list(type = "bar", label = "Middle 95% of 500 random gene sets"), list(type = "tick", label = "Their median", col = K_MED, lwd = .45)), cx, H_KEY - .25)
      grid.text("Numbers: fraction of random sets as low or lower (bold, below 0.05)", unit(cx, "in"), unit(H_KEY - .40, "in"),
                gp = gpar(fontfamily = FONT, fontsize = PT, col = "black"))
    }

    # I: the 78 reproducible scorable Figure 5E regulons, SC-SIC, pseudobulk; each regulon's change minus the median of its
    # own random sets, by the direction of its factor's targets; black tick: the row's median
    KC <- K_C[K_C$scored %in% c("True", TRUE) & K_C$reproducible %in% c("True", TRUE), ]
    stopifnot(nrow(KC) == 78)
    KC$g <- factor(ifelse(KC$factor_direction == "SC-SIC higher", "h", "l"), c("h", "l"))
    K_HALF <- c(h = "Higher in\nSC-SIC", l = "Lower in\nSC-SIC"); K_HL <- unname(K_HALF)
    KC$gy <- 1; KC$half <- factor(K_HALF[as.character(KC$g)], K_HL)
    set.seed(7); KC$jy <- KC$gy + runif(nrow(KC), -.16, .16)   # vertical spread of the other regulons
    KC$cdx <- grepl("__CDX2\\(", KC$regulon)
    K_NM <- c(GSE167880 = "Balboa", GSE185036 = "Du", HICKEY_2023_DRYAD = "Hickey", EGOZI_2021_ZENODO = "Egozi")
    KC$lab <- ifelse(KC$cdx, K_NM[KC$network], "")
    KC$jy[KC$cdx] <- KC$gy[KC$cdx]                                 # CDX2 points on the row center
    stopifnot(sum(KC$cdx) == 2)
    K_MEDS <- aggregate(excess_over_random ~ g, KC, median); K_MEDS$gy <- 1; K_MEDS$half <- factor(K_HALF[as.character(K_MEDS$g)], K_HL)
    K_CNT <- aggregate(cbind(below = excess_over_random < 0, n = 1) ~ g, KC, sum); K_CNT$gy <- 1; K_CNT$half <- factor(K_HALF[as.character(K_CNT$g)], K_HL)
    K_CNT$txt <- sprintf("%d of %d", K_CNT$below, K_CNT$n)
    XB <- range(KC$excess_over_random); XB <- XB + c(-.03, .03) * diff(XB); XT <- XB[2] + .04 * diff(XB)
    pI <- ggplot(KC) +
      geom_vline(xintercept = 0, colour = REF_COL, linewidth = .3) +
      geom_point(data = KC[!KC$cdx, ], aes(excess_over_random, jy), shape = 21, fill = K_OTHER, colour = "black", stroke = STROKE_THIN * .7, size = 1.1) +
      geom_segment(data = K_MEDS, aes(x = excess_over_random, xend = excess_over_random, y = gy - .19, yend = gy + .19), colour = "black", linewidth = .6) +
      geom_point(data = KC[KC$cdx, ], aes(excess_over_random, jy), shape = 21, fill = SC_SIC, colour = "black", stroke = K_ST, size = K_PT + .3) +
      geom_text(data = KC[KC$cdx, ], aes(excess_over_random, jy + .27, label = lab), size = TXT, family = FONT, vjust = 0) +
      geom_text(data = K_CNT, aes(x = XT, y = gy, label = txt), hjust = 0, size = TXT, family = FONT) +
      geom_text(data = data.frame(half = factor(K_HL[1], K_HL)), aes(x = XT, y = 1.6, label = "Below 0"), hjust = 0, size = TXT, family = FONT, fontface = "italic") +
      facet_grid(half ~ ., switch = "y") +
      scale_y_continuous(breaks = NULL, limits = c(.72, 1.78), expand = c(0, 0)) +
      scale_x_continuous(limits = c(XB[1], XB[2] + .15 * diff(XB)), breaks = seq(-.2, .6, .2), labels = minus, expand = c(0, 0), oob = scales::oob_keep) +
      coord_cartesian(clip = "off") +
      labs(x = "CDX2 knockout minus wild type, relative to the\nregulon's own random gene sets (pseudobulk, log2)", y = NULL) +
      theme_panel("x") +
      theme(plot.margin = margin(1, 40, 2, 2),                 # right: room for the counts beside the rows
            strip.placement = "outside", strip.text.y.left = element_text(angle = 90, face = "bold", size = PT, lineheight = .95, hjust = .5),
            panel.background = element_rect(fill = "#F3F3F3", colour = NA), panel.spacing.y = unit(.06, "in"),
            axis.title.x = element_text(margin = margin(2, 0, 0, 0)))
    I_KEY <- .36                                                  # two key rows under I
    drawI <- function(w) {
      vh <- convertHeight(unit(1, "npc"), "in", valueOnly = TRUE)
      pushViewport(viewport(x = 0, y = unit(I_KEY, "in"), width = unit(w, "in"), height = unit(vh - I_KEY, "in"), just = c("left", "bottom")))
      print(pI, newpage = FALSE); popViewport()
      key_items(list(list(type = "pt", label = "CDX2(+) regulon", fill = SC_SIC, size = K_PT + .3), list(type = "tick", label = "Median of the row", col = "black", lwd = .6)),
                w / 2 + .2, I_KEY - .12)
      key_items(list(list(type = "pt", label = "Other reproducible Figure 5E regulon", fill = K_OTHER, size = 1.1, stroke = STROKE_THIN * .7)), w / 2 + .2, I_KEY - .27)
    }
  }, envir = e)
  e
}
# ---- identity factors (pages S18 and S19) -----------------------------------------------------------------------------
# The 14 Figure 5E factors in two halves by direction, then by symbol; italic symbol with a roman (+)/(−) suffix; y counts
# the rows from the bottom.
HALF <- c("Targets higher\nin SC-SIC", "Targets lower\nin SC-SIC")
factor_rows <- function() {
  FA <- unique(read_table("identity/gse214852_values_primary_labels.tsv")[, c("factor", "direction", "d")])
  FA$half <- factor(ifelse(FA$direction == "SC-SIC higher", HALF[1], HALF[2]), HALF); FA$tf <- sub("\\(.*$", "", FA$factor)
  FA <- FA[order(FA$half, FA$tf, FA$factor), ]
  FA$expr <- sprintf('italic("%s")~"(%s)"', FA$tf, ifelse(grepl("\\(\\+\\)$", FA$factor), "+", "−"))
  FA$y <- as.numeric(factor(FA$factor, rev(FA$factor)))
  FA
}
# Page S18: regulon dots colored as in Figure 5E, grey bands on every second factor row, a factor's regulons spread .11
# apart inside its row (at most +-.38)
REG_COL <- c(pass_ec = SC_SIC, pass_beta = SC_BETA, fail = "grey45")
REG_LAB <- c(pass_ec = "higher in SC-SIC", pass_beta = "lower in SC-SIC", fail = "not reproducible")
REG_KEY_TITLE <- "Reproducible in SC-islets and intestine:"
row_bands <- function(FA) transform(FA[seq_len(nrow(FA)) %% 2 == 0, c("factor", "half", "y")], ymin = y - .5, ymax = y + .5)
spread <- function(i) { n <- length(i); if (n == 1) 0 else seq(-1, 1, length.out = n) * min(.38, .11 * (n - 1)) }

# ---- page S19: where the factors' pooled targets fall -------------------------------------------------------------------
# A (factors with targets higher in SC-SIC) and B (lower) share the row pitch, the row-label column and the right-hand
# text column, so both are laid out from both tables. Inches; positions inside a panel are from its content's bottom.
PITCH <- .18; ML <- .04; DX <- .06                  # row pitch; left margin; plot's right end -> row text
GAP_T <- .15                                       # gap between the two text columns
LW_DASH <- .45; PT_SZ <- 2.3; PT_ST <- .5           # dashed line (mm); observed points (ggplot size, stroke)
TICK_H <- .10; LW_TICK <- .6                        # expected count: black vertical tick (height, line width in mm)
targets_layout <- function() {
  V <- read_table("identity/f5_s_gh_values_ALL_ELIGIBLE.tsv"); FA <- factor_rows()
  up <- FA[FA$direction == "SC-SIC higher", ]; lo <- FA[FA$direction == "SC-SIC lower", ]
  U <- V[V$part == "targets_G2_of_G2G3", ]; U <- U[match(up$factor, U$factor), ]
  A <- data.frame(factor = U$factor, expr = up$expr, g2 = as.integer(U$count), n = as.integer(U$total), e2 = U$expected)
  A$g3 <- A$n - A$g2; A$y <- rev(seq_len(nrow(A))); A$txt <- sprintf("%d, %d of %d", A$g2, A$g3, A$n)
  L <- V[V$part == "targets_beta_side", ]; L <- L[match(lo$factor, L$factor), ]
  B <- data.frame(factor = L$factor, expr = lo$expr, k = as.integer(L$count), e = L$expected); B$y <- rev(seq_len(nrow(B)))
  B$txt <- sprintf("%d (%.1f)", B$k, B$e)
  # second text column: each factor's test q (A: the lean toward the group it has more of than expected)
  TS <- read_table("identity/tests.tsv"); TS <- TS[TS$level == "factor", ]
  fq <- function(q) ifelse(q < .001, "<0.001", sprintf("%.2f", q))
  A$lean <- ifelse(A$g2 > A$e2, "G2", "G3")
  A$q <- TS$q[match(paste(A$factor, paste0("lean_", A$lean)), paste(TS$unit, TS$test))]; stopifnot(!anyNA(A$q))
  B$q <- TS$q[match(paste(B$factor, "beta_side"), paste(TS$unit, TS$test))]; stopifnot(!anyNA(B$q))
  A$txt2 <- paste(A$lean, fq(A$q)); B$txt2 <- fq(B$q)
  pushViewport(viewport(gp = gpar(fontfamily = FONT, fontsize = PT)))
  lab <- max(vapply(c(A$expr, B$expr), function(e) convertWidth(stringWidth(parse(text = e)), "in", valueOnly = TRUE), 0))
  popViewport()
  w_lab <- lab + 3 / 72 + .02                       # widest row label + its 3 pt gap + slack
  hdr <- list(A = c("G2, G3 of total", "Lean, q"), B = c("Observed (expected)", "q"))   # bold column headers
  w1 <- str_in(c(A$txt, B$txt, hdr$A[1], hdr$B[1])); w2 <- str_in(c(A$txt2, B$txt2, hdr$A[2], hdr$B[2]))
  mr <- DX + w1 + GAP_T + w2 + .12                  # right columns: the row texts plus >= 0.10 in of white
  list(A = A, B = B, HDR = hdr, X2 = DX + w1 + GAP_T, W_LAB = w_lab, MR = mr, PW = PANEL_W - ML - w_lab - mr)   # PW: plot width, the same in A and B
}
theme_targets <- theme_panel("none") + theme(axis.text.y = element_text(hjust = 1, margin = margin(r = 3)), plot.margin = margin(2, 0, 2, 0))
targets_y <- function(d) scale_y_continuous(breaks = d$y, labels = parse(text = d$expr), limits = c(.5, nrow(d) + .5), expand = c(0, 0), oob = scales::oob_keep)
# the built plot with fixed widths (left margin | row labels | plot | right text column) and plot height rows x pitch
fix_widths <- function(p, nrow, TL) {
  g <- ggplotGrob(p)
  pc <- unique(g$layout$l[g$layout$name == "panel"]); pr <- unique(g$layout$t[g$layout$name == "panel"])
  al <- unique(g$layout$l[g$layout$name == "axis-l"])
  wv <- rep(0, ncol(g)); wv[1] <- ML; wv[al] <- TL$W_LAB; wv[pc] <- TL$PW; wv[ncol(g)] <- TL$MR
  g$widths <- unit(wv, "in"); g$heights[pr] <- unit(nrow * PITCH, "in"); g
}
# One key row centered at cx, its center line y in from the content's bottom: color squares (sq), the dashed line (ln), the
# expected tick (tk) or an observed point (pt), each followed by its label.
KEY_LN <- 5 * 4 * LW_DASH * .pt / 96                # key line: three whole dashes
KEY_LX <- .33                                       # label offset after the key line
key_row <- function(items, cx, y) {
  ww <- vapply(items, function(it) str_in(it$label) + if (it$type == "ln") KEY_LX else .14, 0)
  x <- cx - (sum(ww) + .2 * (length(ww) - 1)) / 2
  for (i in seq_along(items)) {
    it <- items[[i]]
    if (it$type == "sq") grid.rect(unit(x, "in"), unit(y, "in"), unit(.095, "in"), unit(.095, "in"), just = "left", gp = gpar(fill = it$col, col = NA))
    if (it$type == "ln") grid.lines(unit(c(x, x + KEY_LN), "in"), unit(c(y, y), "in"),
                                    gp = gpar(col = "black", lwd = LW_DASH * .pt, lty = "dashed", lineend = "butt"))
    if (it$type == "tk") grid.lines(unit(c(x + .05, x + .05), "in"), unit(y + c(-1, 1) * TICK_H / 2, "in"), gp = gpar(col = "black", lwd = LW_TICK * .pt, lineend = "butt"))
    if (it$type == "pt") grid.points(unit(x + .05, "in"), unit(y, "in"), pch = 21, size = unit(1, "char"),
                                     gp = gpar(col = "grey25", fill = it$fill, fontsize = PT_SZ * .pt + PT_ST * .stroke / 2, lwd = PT_ST * .stroke / 2))
    grid.text(it$label, unit(x + if (it$type == "ln") KEY_LX else .14, "in"), unit(y, "in"), just = c("left", "centre"),
              gp = gpar(fontfamily = FONT, fontsize = PT, col = "black"))
    x <- x + ww[i] + .2
  }
}
# A panel's content: its bold row-group label, the plot g (gh in high, its top vh in above the content's bottom) and the
# key row centered under the plot.
draw_targets <- function(g, gh, vh, items, head, TL) {
  pushViewport(viewport(x = 0, y = unit(vh, "in"), width = unit(PANEL_W, "in"), height = unit(gh, "in"), just = c("left", "top")))
  grid.draw(g); popViewport()
  key_row(items, ML + TL$W_LAB + TL$PW / 2, vh - gh - .10)
  grid.text(head, x = unit(ML, "in"), y = unit(vh + .10 - .08, "in"), just = c("left", "bottom"),
            gp = gpar(fontfamily = FONT, fontsize = PT, fontface = "bold", col = "black"))
}
