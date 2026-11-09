# Figure S7B: intestinal EC's gene-by-gene changes followed in five populations. Each point is a gene judged in both
# intestinal EC and the population: its intestinal EC fold change against other EEC (x, median of three intestinal
# studies) against the population's fold change against its comparison cells (y). Intestinal EC-up genes (up at least
# 2-fold in 2 of 3 intestinal studies) in the population's color, the EEC and islet marker genes black and labeled, other
# genes grey; dashed line: equal fold change; colored line: least-squares fit without the marker genes. Each strip gives
# Spearman rho without the marker genes and the share of intestinal EC-up genes also up 2-fold in the population. Key in
# the sixth cell.

G <- read_table("s3_signature/gene_signatures.tsv.gz")
for (b in grep("judgeable|up_in_half|excluded18", names(G), value = TRUE)) G[[b]] <- G[[b]] == "True"
SC_POP <- c(SC_SIC = "SC-SIC vs other SC endocrine", SC_BETA = "SC-beta vs other SC endocrine\n(not SC-SIC)",
            SC_ALPHA = "SC-alpha vs other SC endocrine\n(not SC-SIC)", INDUCED_EEC = "Induced EEC (Singh 2024, day 6):\nEC vs other EEC",
            HIO_GRAFT_EC = "HIO graft (McCauley 2023):\nEC vs other EEC")
SC_COL <- c(SC_SIC = SIC, SC_BETA = "#4477AA", SC_ALPHA = "#117733", INDUCED_EEC = "#CCBB44", HIO_GRAFT_EC = "#6943AC")
SC_XL <- c(-8.5, 7); SC_YL <- c(-8, 9); SC_LAB_PT <- 7
# one row per gene judged in both intestinal EC and the population; strip = population, rho and the 2-fold share
L <- do.call(rbind, lapply(names(SC_POP), function(p) {
  ok <- G$gut_judgeable & G[[paste0(p, "_judgeable")]] & !is.na(G[[paste0(p, "_lfc")]]) & !is.na(G$gut_ec_lfc)
  d <- data.frame(gene = G$gene, x = G$gut_ec_lfc, y = G[[paste0(p, "_lfc")]], up = G$gut_ec_up_in_half, mk = G$excluded18, pop = p)[ok, ]
  k <- d[!d$mk, ]                                                 # marker genes left out of rho and the fit
  rho <- cor(rank(k$x), rank(k$y)); fit <- coef(lm(y ~ x, k)); n2 <- sum(k$y[k$up] >= 1)
  d$strip <- sub("= -", "= −", sprintf("%s\nρ = %.2f · %d/%d also up 2-fold", SC_POP[[p]], rho, n2, sum(k$up)), fixed = TRUE)
  d$a <- fit[[1]]; d$b <- fit[[2]]; d
}))
L$strip <- factor(L$strip, unique(L$strip))
L$kind <- ifelse(L$mk, "marker", ifelse(L$up, "up", "other"))
FITS <- unique(L[c("strip", "pop", "a", "b")])
# marker-gene labels: label center (lx, ly) and, where lead is TRUE, a leader from (sx0, sy0) to (sx1, sy1), in data
# units; placed once by a search that keeps labels off points, lines and each other (drawing-only file)
LABS <- read.delim(file.path(dirname(PLOT_SCRIPT), "scatter_labels.tsv"), stringsAsFactors = FALSE)
LABS$strip <- factor(as.character(L$strip[match(LABS$pop, L$pop)]), levels(L$strip))

plot_B <- ggplot(L, aes(x, y)) +
  geom_hline(yintercept = 0, colour = REF, linewidth = .25) + geom_vline(xintercept = 0, colour = REF, linewidth = .25) +
  geom_abline(slope = 1, intercept = 0, colour = "grey45", linewidth = .3, linetype = "dashed") +
  geom_point(data = L[L$kind == "other", ], colour = "grey78", size = .25, stroke = 0) +
  geom_point(data = L[L$kind == "up", ], aes(colour = pop), size = .8, stroke = 0) +
  geom_abline(data = FITS, aes(slope = b, intercept = a, colour = pop), linewidth = .5) +
  geom_point(data = L[L$mk, ], colour = "black", size = 1.1, stroke = 0) +
  geom_segment(data = LABS[LABS$lead, ], aes(x = sx0, y = sy0, xend = sx1, yend = sy1), colour = "grey25", linewidth = .2, inherit.aes = FALSE) +
  geom_text(data = LABS, aes(x = lx, y = ly, label = gene), size = SC_LAB_PT / .pt, family = FONT, fontface = "italic",
            colour = "black", hjust = .5, vjust = .5, inherit.aes = FALSE) +
  facet_wrap(~strip, ncol = 3) +
  scale_colour_manual(values = SC_COL, guide = "none") +
  scale_x_continuous(breaks = c(-8, -4, 0, 4), labels = minus) +
  scale_y_continuous(breaks = c(-5, 0, 5), labels = minus) +
  coord_cartesian(xlim = SC_XL, ylim = SC_YL, expand = FALSE) +
  labs(x = expression("Intestinal EC vs other EEC (" * log[2] * " fold change, median of 3 intestinal studies)"),
       y = expression("Population vs its comparison cells (" * log[2] * " fold change)")) +
  theme_rk("both") +
  theme(strip.text = element_text(size = 7.5, face = "bold", hjust = 0, vjust = 0, lineheight = .95, margin = margin(1, 0, 3, 0)),
        panel.spacing.x = unit(.15, "in"), panel.spacing.y = unit(.10, "in"), plot.margin = margin(2, 8, 2, 2))

# key in the empty sixth cell, glyphs as drawn
KEYB <- list(list(g = "dot", col = "grey72", size = 1.2, txt = "Other genes"),
             list(g = "dots5", size = 1.4, txt = "Intestinal EC-up genes (at least\n2-fold in at least 2 of 3 intestinal\nstudies), in the population's color"),
             list(g = "dot", col = "black", size = 1.6, txt = sprintf("EEC and islet marker genes\n(%d; left out of ρ)", sum(G$excluded18))),
             list(g = "dash", txt = "Equal fold change"),
             list(g = "line5", txt = "Fit (marker genes left out),\npopulation's color"))
GLW <- .34; TXT_XB <- .42; LINE_HB <- .128; ROW_GAPB <- .075
key_grob <- function(cell_w) {
  txt <- vapply(KEYB, `[[`, "", "txt")
  nl <- vapply(strsplit(txt, "\n", fixed = TRUE), length, 0L)
  hk <- sum(nl * LINE_HB) + (length(KEYB) - 1) * ROW_GAPB
  kids <- list(); ytop <- hk
  for (i in seq_along(KEYB)) {
    k <- KEYB[[i]]; hb <- nl[i] * LINE_HB; yc <- ytop - hb / 2
    kids[[length(kids) + 1]] <- textGrob(txt[i], x = unit(TXT_XB, "in"), y = unit(yc, "in"), just = c("left", "centre"),
                                         gp = gpar(fontfamily = FONT, fontsize = PT, col = "black", lineheight = .92))
    if (k$g == "dot") kids[[length(kids) + 1]] <- pointsGrob(unit(GLW / 2, "in"), unit(yc, "in"), pch = 16,
                                                             gp = gpar(col = k$col, fontsize = k$size * .pt))
    if (k$g == "dots5") kids[[length(kids) + 1]] <- pointsGrob(unit(seq(.03, GLW - .03, length.out = 5), "in"), unit(rep(yc, 5), "in"), pch = 16,
                                                               gp = gpar(col = unname(SC_COL), fontsize = k$size * .pt))
    if (k$g == "dash") { x0 <- .025 + (0:3) * .084                  # four whole dashes, the plotted line's pattern
      kids[[length(kids) + 1]] <- segmentsGrob(unit(x0, "in"), unit(rep(yc, 4), "in"), unit(x0 + .042, "in"), unit(rep(yc, 4), "in"),
                                               gp = gpar(col = "grey45", lwd = .3 * .pt, lineend = "butt")) }
    if (k$g == "line5") { xs <- seq(0, GLW, length.out = 6)
      kids[[length(kids) + 1]] <- segmentsGrob(unit(xs[-6], "in"), unit(rep(yc, 5), "in"), unit(xs[-1], "in"), unit(rep(yc, 5), "in"),
                                               gp = gpar(col = unname(SC_COL), lwd = .5 * .pt, lineend = "butt")) }
    ytop <- ytop - hb - ROW_GAPB
  }
  # the key's bottom sits on the bottom of the empty sixth cell, 0.10 in from its left
  gTree(children = do.call(gList, kids), vp = viewport(x = unit(.10, "in"), y = unit(0, "npc"), width = unit(cell_w - .10, "in"),
                                                        height = unit(hk, "in"), just = c("left", "bottom")))
}
draw_B <- function() {
  gt <- ggplotGrob(plot_B)                                        # built on the open device, which measures the text
  pl <- gt$layout[grepl("^panel", gt$layout$name), ]
  lR <- max(pl$l); tB <- max(pl$t)                                # the empty sixth cell
  cell_w <- (convertWidth(unit(1, "npc"), "in", valueOnly = TRUE) - fixed_part(gt$widths, convertWidth)) / 3
  gt <- gtable::gtable_add_grob(gt, key_grob(cell_w), t = tB, l = lR, b = tB, r = lR, clip = "off", name = "key-sixth-cell")
  grid.draw(gt)
}

draw_panel <- function() in_place(120, 1746, 7.3, 6.55 - (.30 + (.17 + .13)), 3570, draw_B)
