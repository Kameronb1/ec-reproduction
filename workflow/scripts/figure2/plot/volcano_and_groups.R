# Shared by 2C, 2C_key, 2D and 2D_key (sourced by each): the plots of Figure 2C (plot_C) and 2D (plot_D).
#
# 2C: SC-SIC against intestinal EC for every tested gene (donor-level limma-voom fit, study as a block): difference in
# log2 CPM against -log10 p. Genes with FDR >= 0.05 grey; open rings mark the four EC marker genes. Labels: ggrepel
# places eleven (seeded, with empty-label keep-out points along the FDR line, on the rings and where other labels sit);
# seven sit at fixed positions (positions_c_fixed_labels.tsv and positions_c_fixed_leaders.tsv, placed in data units for
# a 2.29 x 2.153 in plot and kept here at the same distance in inches from their genes). Every label has a thin leader.
# 2D: the three gene groups on 2C's x axis. Rows = in how many of the three intestinal studies the gene is 2-fold up in
# EC; G1 and G2 sit in the top row, G3 in the bottom row, the mixed genes in between. Fourteen genes are labeled.

dv <- read_table("donor_fit.tsv.gz")
dv$difference <- dv$logFC; dv$nlp <- -log10(dv$p)
fdr_p <- max(dv$p[dv$FDR < 0.05])                                    # p at the FDR 0.05 line
d_labels <- read_placement("positions_c_labels.tsv")
lab <- d_labels[d_labels$gene != "", ]; m <- match(lab$gene, dv$gene)                 # the labeled genes sit at their marks
stopifnot(!anyNA(m), max(abs(lab$difference - dv$difference[m]), abs(lab$nlp - dv$nlp[m])) < 1e-9)
d_fixed <- read_placement("positions_c_fixed_labels.tsv")
d_leaders <- read_placement("positions_c_fixed_leaders.tsv")
stopifnot(setequal(d_fixed$gene, d_leaders$gene), !anyDuplicated(d_leaders$gene))   # one leader per fixed label
x_rng <- range(dv$difference); y_rng <- range(dv$nlp)
D_XLIM <- x_rng + c(-1, 1) * (diff(x_rng) * 0.05)                                     # data range plus 5% each side,
D_YLIM <- c(y_rng[1] - 0.05 * diff(y_rng) - 0.75, y_rng[2] + 0.05 * diff(y_rng))   # and 0.75 more below for the ring labels
CD_XTITLE <- "SC-SIC minus intestinal EC\n(log2 CPM)"                                # 2C and 2D share the x limits and title

plot_C <- ggplot(dv, aes(difference, nlp)) +
  geom_hline(yintercept = -log10(fdr_p), linetype = "22", colour = "grey55", linewidth = 0.3) +
  geom_vline(xintercept = 0, colour = "grey75", linewidth = 0.3) +
  geom_point(data = dv[dv$call == "NO_CALL", ], colour = "grey78", size = 0.45, alpha = 0.6) +
  geom_point(data = dv[dv$call != "NO_CALL", ], aes(colour = call), size = 1.1, alpha = 0.85) +
  geom_point(data = dv[dv$label_gene, ], shape = 1, size = 2, colour = "grey10") +
  ggrepel::geom_text_repel(data = d_labels, aes(label = gene, fontface = face), size = MM, colour = "black", family = FONT,
                           segment.size = 0.2, segment.colour = "black", bg.colour = "white", bg.r = 0.12, seed = 20260925,
                           box.padding = unit(0.35, "lines"), min.segment.length = 0, force = 3, max.time = 5, max.iter = 1e5, max.overlaps = Inf) +
  annotate("text", x = max(dv$difference), y = -log10(fdr_p), label = "FDR 0.05", hjust = 1, vjust = -0.5, size = MM, colour = "black", family = FONT) +
  geom_segment(data = d_leaders, aes(x = x, y = y, xend = xend, yend = yend), inherit.aes = FALSE, colour = "black", linewidth = 0.2) +
  ggrepel::geom_text_repel(data = d_fixed, aes(x = x, y = y, label = gene, fontface = face), inherit.aes = FALSE, size = MM, family = FONT,
                           colour = "black", bg.colour = "white", bg.r = 0.12, force = 0, force_pull = 0, max.iter = 0,
                           min.segment.length = Inf, box.padding = 0, point.padding = 0, max.overlaps = Inf, seed = 1) +   # drawn where placed
  scale_colour_manual(values = c(SC_SIC_HIGHER = "#AA3377", NATIVE_EC_HIGHER = "#43A5B6"), guide = "none") +
  coord_cartesian(xlim = D_XLIM, ylim = D_YLIM, expand = FALSE) +
  labs(x = CD_XTITLE, y = expression("−log10" ~ italic(p))) +
  theme_fig2()

# The fixed labels keep their offset from their own gene's point in inches: the offset in data units times the ratio of
# the 2.29 x 2.153 in plot's panel size to this plot's, per axis.
local({
  b <- at("C_plot"); p2 <- panel_box(plot_C, 2.29, 2.153); p3 <- panel_box(plot_C, b[["w"]], b[["h"]])
  kx <- (p2$r - p2$l) / (p3$r - p3$l); ky <- (p2$b - p2$t) / (p3$b - p3$t)
  stopifnot(identical(plot_C$layers[[8]]$data, d_leaders), identical(plot_C$layers[[9]]$data, d_fixed))
  gx <- dv$difference[match(d_fixed$gene, dv$gene)]; gy <- dv$nlp[match(d_fixed$gene, dv$gene)]
  plot_C$layers[[9]]$data <<- transform(d_fixed, x = gx + (x - gx) * kx, y = gy + (y - gy) * ky)
  lx <- dv$difference[match(d_leaders$gene, dv$gene)]; ly <- dv$nlp[match(d_leaders$gene, dv$gene)]
  plot_C$layers[[8]]$data <<- transform(d_leaders, x = lx + (x - lx) * kx, y = ly + (y - ly) * ky,
                                        xend = lx + (xend - lx) * kx, yend = ly + (yend - ly) * ky)
  # keep-out points where the fixed labels now sit (center and both ends of each label), so ggrepel keeps its labels off them
  upi <- diff(p3$xr) / (p3$r - p3$l); f <- plot_C$layers[[9]]$data
  png(tempfile(fileext = ".png"), 7.5, 10, units = "in", res = DPI, type = "cairo")
  hw <- mapply(function(s, fc) str_w(s, face = fc), f$gene, f$face) / 2 * upi; invisible(dev.off())
  keep <- do.call(rbind, lapply(c(-1, 0, 1), function(s) data.frame(gene = "", difference = f$x + s * hw, nlp = f$y, face = "plain")))
  stopifnot(identical(names(keep), names(plot_C$layers[[6]]$data)))
  plot_C$layers[[6]]$data <<- rbind(plot_C$layers[[6]]$data, keep)
  # INS sits at the top edge of the panel: a fixed label left of its point (right end 0.05 in from it) with a leader; its
  # repelled label becomes an empty keep-out point at the same spot
  r <- which(plot_C$layers[[6]]$data$gene == "INS"); stopifnot(length(r) == 1)
  gx <- dv$difference[dv$gene == "INS"]; gy <- dv$nlp[dv$gene == "INS"]
  png(tempfile(fileext = ".png"), 7.5, 10, units = "in", res = DPI, type = "cairo"); w <- str_w("INS", face = "bold.italic") * upi; invisible(dev.off())
  cx <- gx - 0.05 * upi - w / 2
  plot_C$layers[[6]]$data$gene[r] <<- ""
  plot_C$layers[[9]]$data <<- rbind(plot_C$layers[[9]]$data, data.frame(gene = "INS", face = "bold.italic", x = cx, y = gy))
  plot_C$layers[[8]]$data <<- rbind(plot_C$layers[[8]]$data, data.frame(gene = "INS", x = cx + w / 2 + 0.01 * upi, y = gy, xend = gx - 0.012 * upi, yend = gy))
})

# ---- 2D ----
GKEY <- c(G1_EC_GENES_SC_SIC_LACK = "G1", G2_EC_GENES_SC_SIC_HAVE = "G2", G3_SC_SIC_OWN_PROGRAM = "G3")   # table values
F_ROWS <- c("0 or 1 of 3", "2 of 3 (mixed)", "All 3 intestinal studies")              # bottom to top
gr <- read_table("gene_groups.tsv.gz"); mx <- read_table("mixed_genes.tsv")
fd <- rbind(data.frame(gene = gr$gene, g = unname(GKEY[gr$group]), row = ifelse(gr$gut_up_of_3 == 3, F_ROWS[3], F_ROWS[1])),
            data.frame(gene = mx$gene, g = "Mixed", row = F_ROWS[2]))
fd$x <- dv$logFC[match(fd$gene, dv$gene)]
fd$row <- factor(fd$row, F_ROWS); fd$g <- factor(fd$g, names(GROUP_COLS))
set.seed(20261004); fd$yj <- as.integer(fd$row) + runif(nrow(fd), -0.24, 0.24)          # spread within the row
n <- table(fd$g)
F_YLAB <- c(sprintf("0 or 1 of 3\nG3: %d", n[["G3"]]), sprintf("2 of 3\n(mixed)\n%d, no group", n[["Mixed"]]),
            sprintf("All 3\nG1: %d\nG2: %d", n[["G1"]], n[["G2"]]))
# the labeled genes: the four EC markers and the genes the text names
D_LABEL_GENES <- list(G1 = c("REG4", "FOLH1", "CEACAM20", "HOXA10", "OR51E1"),
                      G2 = c("TPH1", "SLC18A1", "FEV", "LMX1A", "DDC", "ASCL1"),
                      G3 = c("NEUROD4", "GLP1R", "STAC"))
stopifnot(all(unlist(lapply(names(D_LABEL_GENES), function(g) D_LABEL_GENES[[g]] %in% fd$gene[fd$g == g]))))
# Label placement, a fixed rule computed for the drawn panel size (pw x ph in): labels sit in lanes that hold no points.
# Top row (G1, G2): lanes stacked above the row and under it (between the top and the middle row). In x order each label
# takes the lane with the lowest cost: its sideways shift from its gene (shifts over 0.12 in count five times), plus the
# equivalent of 0.5 in when its leader would pass through a label of a nearer lane on the same side, or when it would sit
# across the leader of a label in a farther lane. G3: the leftmost labeled gene in a column left of the G3 points, the
# rightmost in a column right of them if it fits inside the panel, the rest in the lane above the bottom row. Within a lane
# labels start centered on their gene and are pushed apart (merged clusters centered on their mean) until none overlap.
D_YTOP <- 3.70                                      # top of the y axis; leaves room above the top row
D_LANE <- c(g3_above = 1.53)                        # lane center (data units); the top-row lanes follow the text height
pack_1d <- function(want, w, gap, lo, hi) {           # centers closest to want with no overlap, inside [lo, hi]
  o <- order(want); want <- want[o]; w <- w[o]
  cl <- lapply(seq_along(want), function(i) i)
  repeat {
    pos <- unlist(lapply(cl, function(k) { tot <- sum(w[k]) + gap * (length(k) - 1); m <- mean(want[k])
      m <- min(max(m, lo + tot / 2), hi - tot / 2); m - tot / 2 + cumsum(w[k] + gap) - gap - w[k] / 2 }))
    j <- which(diff(pos) < (w[-1] + w[-length(w)]) / 2 + gap - 1e-9 & diff(rep(seq_along(cl), lengths(cl))) != 0)
    if (!length(j)) break
    a <- rep(seq_along(cl), lengths(cl))[j[1]]; cl[[a]] <- c(cl[[a]], cl[[a + 1]]); cl[[a + 1]] <- NULL
  }
  pos[order(o)]
}
d_label_layout <- function(pw, ph) {
  upx <- diff(D_XLIM) / pw; upy <- (D_YTOP - 0.55) / ph   # data units per inch
  g <- fd[fd$gene %in% unlist(D_LABEL_GENES), c("gene", "g", "x", "yj")]; g <- g[order(g$x), ]
  png(tempfile(fileext = ".png"), 3, 3, units = "in", res = DPI, type = "cairo")
  g$w <- vapply(g$gene, function(s) str_w(s, face = "italic"), 0) * upx; th <- PT / 72 * upy; invisible(dev.off())   # th: label height (one em)
  gap <- 0.05 * upx; g$cx <- NA_real_; g$cy <- NA_real_; step <- th * 1.15
  nA <- max(1, min(3, floor((D_YTOP - 0.01 - 3.26) / step))); nB <- max(1, min(3, floor((2.74 - 2.255) / step)))
  tl <- rbind(data.frame(lane = paste0("top_above", seq_len(nA)), side = 1, depth = seq_len(nA), y = 3.26 + step * (seq_len(nA) - 0.5)),
              data.frame(lane = paste0("top_below", seq_len(nB)), side = -1, depth = seq_len(nB), y = 2.74 - step * (seq_len(nB) - 0.5)))
  tl <- tl[order(tl$depth, -tl$side), ]                # nearest lanes first, above before below
  lanes <- c(setNames(tl$y, tl$lane), D_LANE)
  g$lane <- NA_character_; right <- setNames(rep(-Inf, nrow(tl)), tl$lane)
  placed <- data.frame(side = numeric(), depth = numeric(), l = numeric(), r = numeric(), px = numeric()); pad <- 0.03 * upx
  for (i in which(g$g %in% c("G1", "G2"))) {          # g is in x order
    shift <- pmax(g$x[i], right + gap + g$w[i] / 2) - g$x[i]
    bad <- vapply(seq_len(nrow(tl)), function(m) { c0 <- g$x[i] + shift[[m]]; l <- c0 - g$w[i] / 2 - pad; r <- c0 + g$w[i] / 2 + pad
      same <- placed$side == tl$side[m]
      any(same & placed$depth < tl$depth[m] & placed$l <= g$x[i] & placed$r >= g$x[i]) || any(same & placed$depth > tl$depth[m] & placed$px >= l & placed$px <= r) }, TRUE)
    over <- g$x[i] + shift + g$w[i] / 2 > D_XLIM[2] - gap
    cost <- shift + 4 * pmax(0, shift - 0.12 * upx) + bad * 0.5 * upx + over * 100 * upx
    m <- which.min(cost); ln <- tl$lane[m]; g$lane[i] <- ln; c0 <- g$x[i] + shift[[m]]; right[[ln]] <- c0 + g$w[i] / 2
    placed <- rbind(placed, data.frame(side = tl$side[m], depth = tl$depth[m], l = c0 - g$w[i] / 2 - pad, r = c0 + g$w[i] / 2 + pad, px = g$x[i])) }
  g3 <- which(g$g == "G3"); g3x <- range(fd$x[fd$g == "G3"])
  g$lane[head(g3, 1)] <- "g3_left"
  if (g3x[2] + 0.06 * upx + g$w[tail(g3, 1)] <= D_XLIM[2]) g$lane[tail(g3, 1)] <- "g3_right"
  g$lane[is.na(g$lane)] <- "g3_above"
  for (ln in setdiff(names(lanes), character())) { i <- which(g$lane == ln); if (!length(i)) next
    g$cx[i] <- pack_1d(g$x[i], g$w[i], gap, D_XLIM[1] + gap, D_XLIM[2] - gap); g$cy[i] <- lanes[[ln]] }
  vgap <- 0.03 * upy
  i <- which(g$lane == "g3_left"); g$cx[i] <- g3x[1] - 0.06 * upx - g$w[i] / 2; g$cy[i] <- pack_1d(g$yj[i], rep(th, length(i)), vgap, 0.6, 1.45)
  i <- which(g$lane == "g3_right"); g$cx[i] <- g3x[2] + 0.06 * upx + g$w[i] / 2; g$cy[i] <- pack_1d(g$yj[i], rep(th, length(i)), vgap, 0.6, 1.45)
  # leaders: from the label edge facing the point to the point (lanes: the top or bottom edge, at the point's x where the
  # label spans it; side columns: the inner edge)
  side <- g$lane %in% c("g3_left", "g3_right")
  sx <- ifelse(side, g$cx + ifelse(g$lane == "g3_left", 1, -1) * (g$w / 2 + 0.01 * upx), pmin(pmax(g$x, g$cx - g$w / 2 + 0.03 * upx), g$cx + g$w / 2 - 0.03 * upx))
  sy <- ifelse(side, g$cy, g$cy + ifelse(g$yj > g$cy, 1, -1) * (th / 2 + 0.01 * upy))
  list(labels = data.frame(gene = g$gene, x = g$cx, y = g$cy, w = g$w, h = th, lane = g$lane),
       leaders = data.frame(gene = g$gene, x = sx, y = sy, xend = g$x, yend = g$yj))
}
# Row labels (the number of studies in bold): drawn twice over the axis text, which stays (transparent) so the layout
# does not move; the first line in bold, the lines under it plain, each with the other lines blank. The axis text ends
# 3.36 pt left of the panel; the y range is 0.55 to D_YTOP.
f_row_label <- function(i, first) {
  ln <- strsplit(F_YLAB[i], "\n", fixed = TRUE)[[1]]; ln[if (first) -1 else 1] <- " "
  annotation_custom(textGrob(paste(ln, collapse = "\n"), x = unit(0, "npc") - unit(3.36, "pt"), y = unit((i - .55) / (D_YTOP - .55), "npc"),
                             hjust = 1, gp = gpar(fontfamily = FONT, fontsize = PT, col = "black", lineheight = .85,
                                                  fontface = if (first) "bold" else "plain")))
}

plot_D <- ggplot(fd, aes(x, yj, colour = g)) +
  geom_vline(xintercept = 0, colour = "grey55", linewidth = 0.25) +
  geom_point(size = 0.95, alpha = 0.9, stroke = 0) +
  scale_colour_manual(values = GROUP_COLS, guide = "none") +
  scale_y_continuous(breaks = 1:3, labels = F_YLAB, limits = c(0.55, D_YTOP), expand = c(0, 0)) +
  scale_x_continuous(breaks = c(-10, -5, 0, 5, 10)) +
  lapply(1:3, f_row_label, first = TRUE) + lapply(1:3, f_row_label, first = FALSE) + coord_cartesian(xlim = D_XLIM, expand = FALSE, clip = "off") +
  labs(x = CD_XTITLE, y = NULL) +
  theme_fig2() + theme(panel.grid.major.y = element_blank(), axis.text.y = element_text(lineheight = 0.85, colour = "transparent"))
# the gene labels and their leaders, placed for 2D's drawn panel size, drawn with a white halo
local({
  b <- at("D_plot"); pb <- panel_box(plot_D, b[["w"]], b[["h"]]); L <- d_label_layout(pb$r - pb$l, pb$b - pb$t)
  plot_D <<- plot_D + geom_segment(data = L$leaders, aes(x = x, y = y, xend = xend, yend = yend), inherit.aes = FALSE, colour = "grey40", linewidth = 0.2) +
    ggrepel::geom_text_repel(data = L$labels, aes(x = x, y = y, label = gene), inherit.aes = FALSE, size = MM, family = FONT, fontface = "italic",
                             colour = "black", bg.colour = "white", bg.r = 0.12, force = 0, force_pull = 0, max.iter = 0,
                             min.segment.length = Inf, box.padding = 0, point.padding = 0, max.overlaps = Inf, seed = 1)   # drawn where placed
})
