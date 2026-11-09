# Figure S8A: every expressed gene's HIO graft EC minus native intestinal EC difference (x) against its SC-SIC minus
# native intestinal EC difference (y), log2 CPM, on equal axes, colored by its three-way pattern. The plane is tinted where the
# effect-size part of each pattern rule holds (the closest pair of the three groups and the 1 log2 threshold), with the
# wedge edges drawn; dashed line: equal difference. Wedge names give each pattern's gene count. The top-scoring genes of
# each pattern and direction, the genes of two GO terms and eight genes in the dense center are labeled.

suppressPackageStartupMessages(library(data.table))
main <- as.data.table(read_table("s4_hio/genes_all.tsv.gz"))
pc <- c(STEM_SHARED = SIC, GUT_SHARED = "#43A5B6", HIO_ONLY = "#6943AC", none = "#C8C8C8")
pl <- main[expressed == TRUE]; pl[, pattern := factor(pattern, names(pc))]       # every expressed gene is drawn
cnt <- table(pl$pattern)
fmt <- function(n) format(as.integer(n), big.mark = ",", trim = TRUE)
W_IN <- 7.3; H_IN <- 5.85
XLIM <- c(-13, 12); LIM <- c(-13.5, 16.5)                # x: data plus room for labels; y: data range plus 0.5
PH <- H_IN - .47 - .42                                   # panel height: the plot area minus the x axis and margins
UNIT_IN <- PH / diff(LIM); PW <- UNIT_IN * diff(XLIM)    # inches per log2 unit, the same on both axes
# tint: the effect-size part of each pattern rule (closest pair of groups, differences of at least 1 log2)
T_LOG2 <- 1
region <- function(x, y) { hs <- x - y
  cl <- c("HS", "HP", "SP")[max.col(-cbind(abs(hs), abs(x), abs(y)), ties.method = "first")]
  ifelse(cl == "HS" & sign(x) == sign(y) & abs(x) >= T_LOG2 & abs(y) >= T_LOG2, "STEM_SHARED",
  ifelse(cl == "HP" & sign(y) == -sign(hs) & abs(y) >= T_LOG2 & abs(hs) >= T_LOG2, "GUT_SHARED",
  ifelse(cl == "SP" & sign(x) == sign(hs) & abs(x) >= T_LOG2 & abs(hs) >= T_LOG2, "HIO_ONLY", "none"))) }
GR <- CJ(x = seq(XLIM[1], XLIM[2], by = .04), y = seq(LIM[1], LIM[2], by = .04)); GR[, r := region(x, y)]; GR <- GR[r != "none"]
tint <- c(STEM_SHARED = SIC, GUT_SHARED = "#43A5B6", HIO_ONLY = "#6943AC")
# wedge edges: six rays from the origin (y = 2x, y = x/2, y = -x), each drawn from where a tinted wedge starts
EDGE <- rbindlist(lapply(list(c(1, 2), c(-1, -2), c(2, 1), c(-2, -1), c(1, -1), c(-1, 1)), function(d) {
  u <- d / sqrt(sum(d^2)); n <- c(-u[2], u[1]); t <- seq(0, 40, by = .01)
  px <- t * u[1]; py <- t * u[2]
  tinted <- region(px + .02 * n[1], py + .02 * n[2]) != "none" | region(px - .02 * n[1], py - .02 * n[2]) != "none"
  inside <- px >= XLIM[1] & px <= XLIM[2] & py >= LIM[1] & py <= LIM[2]
  k <- which(tinted & inside)
  data.table(x = px[min(k)], y = py[min(k)], xend = px[max(k)], yend = py[max(k)]) }))
# wedge names with their gene counts
WL <- data.table(r = c("GUT_SHARED", "STEM_SHARED", "HIO_ONLY"),
  lab = c(sprintf("SC-SIC differs; HIO closer to\nnative intestinal EC (%s)", fmt(cnt[["GUT_SHARED"]])),
          sprintf("HIO and\nSC-SIC differ\ntogether (%s)", fmt(cnt[["STEM_SHARED"]])),
          sprintf("HIO differs;\nSC-SIC closer to native\nintestinal EC (%s)", fmt(cnt[["HIO_ONLY"]]))),
  x = c(XLIM[1] + .3, XLIM[2] - .3, XLIM[2] - .3), y = c(LIM[2] - .35, 13.6, -4.85), hj = c(0, 1, 1))

# ---- gene labels, with leaders where a label sits away from its point
# label-box centers (xc, yc, log2 units) of the labeled genes (drawing-only file)
PAD <- .016; INK_H <- .072                               # inches: halo + gap around the ink; cap height of 7 pt Arial
LAYOUT <- fread(file.path(dirname(PLOT_SCRIPT), "hio_label_places.tsv"))
lab <- pl[gene %in% LAYOUT$gene]
lab[, w := (str_w_at(gene, 100, 7, "italic") + 2 * PAD) / UNIT_IN]; BH <- (INK_H + 2 * PAD) / UNIT_IN
GAP_LEAD <- .055 / UNIT_IN                               # a label farther than this from its point gets a leader
LB <- merge(LAYOUT, lab[, .(gene, HP, SP, w)], by = "gene")
LB[, `:=`(xmin = xc - w / 2, xmax = xc + w / 2, ymin = yc - BH / 2, ymax = yc + BH / 2)]
LB[, `:=`(ex = pmin(pmax(HP, xmin), xmax), ey = pmin(pmax(SP, ymin), ymax))]           # nearest point of the label box
LB[, dist_in := sqrt((ex - HP)^2 + (ey - SP)^2) * UNIT_IN]
LB[, `:=`(lx0 = ifelse(dist_in / UNIT_IN > GAP_LEAD, HP, NA_real_), ly0 = ifelse(dist_in / UNIT_IN > GAP_LEAD, SP, NA_real_),
          lx1 = ifelse(dist_in / UNIT_IN > GAP_LEAD, ex, NA_real_), ly1 = ifelse(dist_in / UNIT_IN > GAP_LEAD, ey, NA_real_))]
LD <- LB[!is.na(lx0)]
LD[, `:=`(len = sqrt((lx1 - lx0)^2 + (ly1 - ly0)^2))]
LD[, `:=`(sx = lx0 + (lx1 - lx0) / len * (.016 / UNIT_IN), sy = ly0 + (ly1 - ly0) / len * (.016 / UNIT_IN))]   # leader from the point's edge
# white halos: each label drawn 16 times in white around its place first
HALO <- data.table(a = seq(0, 2 * pi, length.out = 17)[-17])[, `:=`(dx = cos(a) * .011 / UNIT_IN, dy = sin(a) * .011 / UNIT_IN)]
LH <- LB[, .(gene, xc, yc)][, k := 1][HALO[, k := 1], on = "k", allow.cartesian = TRUE][, `:=`(x = xc + dx, y = yc + dy)]

plot_A <- ggplot(pl[order(pattern != "none")], aes(HP, SP)) +
  geom_raster(data = GR, aes(x, y, fill = r), inherit.aes = FALSE, alpha = .07, show.legend = FALSE) +
  scale_fill_manual(values = tint, guide = "none") +
  geom_hline(yintercept = 0, colour = "grey55", linewidth = .25) + geom_vline(xintercept = 0, colour = "grey55", linewidth = .25) +
  geom_segment(data = EDGE, aes(x = x, y = y, xend = xend, yend = yend, linetype = "Closest pair changes"), inherit.aes = FALSE, colour = "grey70", linewidth = .25) +
  geom_abline(aes(slope = 1, intercept = 0, linetype = "Equal difference"), colour = "grey55", linewidth = .3, key_glyph = "path") +
  geom_point(aes(colour = pattern), shape = 16, size = .7, alpha = .8) +
  geom_text(data = WL, aes(x, y, label = lab, hjust = hj), inherit.aes = FALSE, colour = tint[WL$r], size = 7 / .pt,
            family = FONT, fontface = "bold", lineheight = .9, vjust = 1) +
  geom_segment(data = LD, aes(x = sx, y = sy, xend = lx1, yend = ly1), inherit.aes = FALSE, colour = "grey30", linewidth = .2) +
  geom_text(data = LH, aes(x, y, label = gene), inherit.aes = FALSE, size = 7 / .pt, family = FONT, fontface = "italic", colour = "white") +
  geom_text(data = LB, aes(xc, yc, label = gene), inherit.aes = FALSE, size = 7 / .pt, family = FONT, fontface = "italic", colour = "grey10") +
  scale_colour_manual(values = pc, breaks = "none", labels = sprintf("No pattern (%s)", fmt(cnt[["none"]])), name = NULL) +
  scale_linetype_manual(values = c("Equal difference" = "22", "Closest pair changes" = "solid"), name = NULL) +
  labs(x = expression("HIO graft EC minus native intestinal EC ("*log[2]*" CPM)"), y = expression("SC-SIC minus native intestinal EC ("*log[2]*" CPM)")) +
  guides(colour = guide_legend(order = 1, override.aes = list(size = 2.2, alpha = 1)),
         linetype = guide_legend(order = 2, override.aes = list(colour = c("grey55", "grey70"), linewidth = c(.4, .4)))) +
  coord_equal(xlim = XLIM, ylim = LIM, expand = FALSE, clip = "on") +
  scale_x_continuous(breaks = seq(-20, 20, 5), labels = minus) + scale_y_continuous(breaks = seq(-20, 20, 5), labels = minus) +
  theme_rk("none", legend = "right") +
  theme(axis.line.y = element_line(colour = "black", linewidth = .35), axis.ticks.y = element_line(colour = "black", linewidth = .3),
        legend.justification = c(0, 1), legend.key.width = unit(.3, "in"), legend.spacing.y = unit(2, "pt"),
        legend.box.just = "left", legend.margin = margin(0, 0, 0, 4),
        panel.widths = unit(PW, "in"), panel.heights = unit(PH, "in"))
draw_A <- function() {                                   # left- and top-aligned at the plot's own size
  gt <- ggplotGrob(plot_A); wt <- convertWidth(sum(gt$widths), "in", valueOnly = TRUE)
  ht <- convertHeight(sum(gt$heights), "in", valueOnly = TRUE)
  pushViewport(viewport(x = 0, y = 1, width = unit(wt, "in"), height = unit(ht, "in"), just = c("left", "top")))
  grid.draw(gt); popViewport() }

draw_panel <- function() in_place(120, 48, W_IN, H_IN - .47, 3228, draw_A)
