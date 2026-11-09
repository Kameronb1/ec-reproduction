# Figure S15A: the 20 genes that best separate SC-SIC from other SC endocrine cells (rows, in their limma rank with the
# six SC-islet studies weighted equally) in each study (columns), as a dot plot like Figure 2D: color = the study's mean
# log expression in SC-SIC minus other SC endocrine cells; size = share of the study's SC-SICs with >= 1 read (mean over
# its samples). Color bar and size key below.

eff <- read_table("s1_marker_limma/heatmap_effects.tsv")
D20 <- unique(eff$gene_symbol[order(eff$panel_rank)])          # the 20 genes, ranked

LO <- "#2F6DA8"; HI <- "#8E3D29"; LIM <- c(-2, 6)                 # blue below 0, brown above; beyond 6 drawn at 6
ny <- length(D20); KEYH <- .66; D_MAX <- 4.0                     # dot diameter (mm) for 100% of SC-SICs (Fig. 2D: 3.6)
WA <- XT + NS * CW + .10                                         # plot width and height (inches)
HA <- HEAD + ny * P + .06 + LBLH + .14 + KEYH + .04
y0 <- HA - HEAD; yb <- y0 - ny * P                               # top and bottom of the gene grid
# dot columns tighter than B's (as Fig. 2D), centered in A's width; 4.0 mm dots nearly fill a row (pitch 0.16 in)
CWA <- .24; OFF <- NS * (CW - CWA) / 2; XTA <- XT + OFF; XGA <- XG + OFF
eff$x <- XTA + (match(eff$family_id, names(STUDY)) - .5) * CWA
eff$y <- y0 - (match(eff$gene_symbol, D20) - .5) * P
gl <- data.frame(x = XGA, y = y0 - (seq_len(ny) - .5) * P, lab = D20)
sl <- data.frame(x = XTA + (seq_len(NS) - .5) * CWA, y = yb - .06, lab = unname(STUDY))
KT <- yb - .06 - LBLH - .14; BW <- 1.50; KX0 <- WA / 2 - BW / 2  # key: title, color bar, ticks
pal <- function(v) scales::div_gradient_pal(LO, "white", HI, "Lab")(scales::rescale_mid(scales::squish(v, LIM), to = c(0, 1), from = LIM, mid = 0))
vv <- seq(LIM[1], LIM[2], length.out = 241); vv <- (head(vv, -1) + tail(vv, -1)) / 2   # 240 color-bar steps
ticks <- data.frame(x = KX0 + BW * (c(-2, 0, 2, 4, 6) - LIM[1]) / diff(LIM), lab = c("−2", "0", "2", "4", "6"))
# size key under the color bar, centered: title, then dots for 25%, 50% and 100%
SK <- "SC-SICs detecting it:"; sz <- c(.25, .5, 1); sl_k <- c("25%", "50%", "100%")
sw <- str_w_at(c(SK, sl_k), 100); KY <- KT - .52
sx <- WA / 2 - (sw[1] + .08 + sum(sw[-1] + .22) - .07) / 2
kx <- sx + sw[1] + .08 + c(0, cumsum(sw[-1] + .22))[1:3]
skey <- data.frame(x = kx + .06, y = KY, s = D_MAX * sqrt(sz)); slab <- data.frame(x = kx + .15, y = KY, lab = sl_k)

plot_A <- ggplot() +
  geom_point(data = eff, aes(x, y, fill = effect, size = sic_detection), shape = 21, colour = "black", stroke = .15) +
  scale_fill_gradient2(low = LO, mid = "white", high = HI, midpoint = 0, limits = LIM, oob = scales::squish, guide = "none") +
  scale_size_area(max_size = D_MAX, limits = c(0, 1), guide = "none") +
  tx(data = gl, aes(x, y, label = lab), hjust = 1, vjust = .5, fontface = "italic") +
  tx(data = sl, aes(x, y, label = lab), angle = ANG, hjust = 1, vjust = .5) +
  tx(data = data.frame(x = XTA + NS * CWA / 2, y = y0 + .05, lab = "SC-islet study"), aes(x, y, label = lab), hjust = .5, vjust = 0) +
  tx(data = data.frame(x = WA / 2, y = KT, lab = "SC-SIC minus other SC endocrine (log expression)"), aes(x, y, label = lab), hjust = .5, vjust = 1) +
  annotation_raster(matrix(pal(vv), nrow = 1), xmin = KX0, xmax = KX0 + BW, ymin = KT - .15 - .09, ymax = KT - .15, interpolate = FALSE) +
  tx(data = ticks, aes(x, y = KT - .27, label = lab), hjust = .5, vjust = 1) +
  tx(data = data.frame(x = sx, y = KY, lab = SK), aes(x, y, label = lab), hjust = 0, vjust = .5) +
  geom_point(data = skey, aes(x, y), size = skey$s, shape = 21, fill = "#8C8C8C", colour = "black", stroke = .15, inherit.aes = FALSE) +
  tx(data = slab, aes(x, y, label = lab), hjust = 0, vjust = .5) +
  scale_x_continuous(limits = c(0, WA), expand = c(0, 0)) + scale_y_continuous(limits = c(0, HA), expand = c(0, 0)) +
  coord_cartesian(clip = "off") + theme_void(base_family = FONT) + theme(plot.margin = margin(0, 0, 0, 0), legend.position = "none")

draw_panel <- function() in_place(120, 48, WA, HA, 3322, function() print(plot_A, newpage = FALSE))
