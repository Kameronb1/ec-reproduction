# Figure 3E: nine canonical EC markers. Points = each study's level (triangles intestinal EC studies, circles SC-SIC studies);
# short black ticks behind them = each arm's model estimate from Figure 3D's fit; right axis = 3D's FDR for the gene.
# The four classifier genes sit above the grey rule, under a bracket. Column headers above, two key rows below.

MARKERS <- c("TPH1", "SLC18A1", "FEV", "LMX1A", "DDC", "GCH1", "CHGA", "CHGB", "TAC1")    # top to bottom; first four = the classifier
row_y <- setNames(rev(seq_along(MARKERS)), MARKERS)
ARM_OFF <- c("SC-SIC study" = 0.17, "primary EC study" = -0.17)                        # each arm on its own track within a row
E_COLS <- c("SC-SIC study" = "#AA3377", "primary EC study" = "#43A5B6")
TICK_HALF <- 0.22; TICK_LW <- 0.4                                                       # tick half-length (rows); line width (ggplot mm)
sl <- read_table("panel_e_per_study_levels.tsv.gz")
sl <- sl[sl$gene %in% MARKERS, ]
e_pts <- data.frame(kind = ifelse(sl$arm == "SC_SIC", "SC-SIC study", "primary EC study"), v = sl$level, y = row_y[sl$gene])
e_pts$y <- e_pts$y + ARM_OFF[e_pts$kind]
dv <- read_table("panel_d_donor_fit.tsv.gz")
fit <- dv[match(MARKERS, dv$gene), ]
e_tick <- data.frame(y = rep(row_y[MARKERS], 2), kind = rep(c("SC-SIC study", "primary EC study"), each = length(MARKERS)),
                     v = c(fit$sc_level, fit$primary_level))
e_tick$y0 <- e_tick$y + ARM_OFF[e_tick$kind] - TICK_HALF; e_tick$y1 <- e_tick$y + ARM_OFF[e_tick$kind] + TICK_HALF
e_fdr <- sprintf("%#.2g", fit$FDR)
RULE_Y <- 5.5                                                                            # the grey rule under the classifier genes
CLS_Y <- c(min(row_y[MARKERS[1:4]]) - 0.17 - TICK_HALF, max(row_y[MARKERS[1:4]]) + 0.17 + TICK_HALF)   # their rows' marks, bottom to top

# The ticks are plain black lines drawn before the points, so the points sit on top. The bracket is placed in inches
# from the panel's right edge, so the plot is built, its panel measured, and built again with the bracket.
TICK_W_IN <- TICK_LW * .pt / 96                                                          # tick line width (in)
E_BR <- list(gap = 0.02, arm = 0.03, lab_gap = 0.035)                                   # bracket: in from the panel's right edge, end arms, label gap (in)
e_plot <- function(br = NULL) {
  p <- ggplot(e_pts, aes(v, y)) +
    geom_segment(data = e_tick, aes(x = v, xend = v, y = y0, yend = y1), colour = "black", linewidth = TICK_LW, lineend = "butt",
                 inherit.aes = FALSE) +
    geom_point(aes(fill = kind, shape = kind), colour = "black", size = 1.3, stroke = 0.15) +
    geom_hline(yintercept = RULE_Y, colour = "grey70", linewidth = 0.3)                 # classifier genes above the rule
  if (!is.null(br)) p <- p +
    annotate("path", x = c(br$x - br$arm, br$x, br$x, br$x - br$arm), y = c(CLS_Y[1], CLS_Y[1], CLS_Y[2], CLS_Y[2]), colour = "black", linewidth = 0.3) +
    annotate("text", x = br$x_lab, y = mean(CLS_Y), label = "Four-marker\nclassifier", hjust = 1, vjust = 0.5, lineheight = 0.9,
             size = MM, family = FONT, colour = "black")
  p + scale_fill_manual(values = E_COLS) +
    scale_shape_manual(values = c("SC-SIC study" = 21, "primary EC study" = 24)) +
    scale_x_continuous(limits = c(floor(min(e_pts$v)) - 0.5, ceiling(max(e_pts$v)) + 0.5), breaks = c(4, 8, 12, 16), expand = expansion(add = c(0.2, 0))) +
    scale_y_continuous(breaks = row_y, labels = MARKERS, expand = expansion(add = 0.6), sec.axis = dup_axis(labels = e_fdr, name = NULL)) +
    labs(x = "log2 CPM", y = NULL) +
    theme_fig3() + theme(panel.grid.major.y = element_blank(), axis.text.y.left = element_text(size = 8, face = "italic"),
                         axis.text.y.right = element_text(size = 7.5), plot.margin = margin(3, 3, 1, 1))
}
PLOT_W <- 2.62; PLOT_H <- 1.94
E_PB <- panel_box(e_plot(), PLOT_W, PLOT_H)
E_UPI <- c(x = diff(E_PB$xr) / (E_PB$r - E_PB$l), y = diff(E_PB$yr) / (E_PB$b - E_PB$t))   # data units per inch
E_BR$x <- E_PB$xr[2] - E_BR$gap * E_UPI[["x"]]; E_BR$arm <- E_BR$arm * E_UPI[["x"]]; E_BR$x_lab <- E_BR$x - E_BR$lab_gap * E_UPI[["x"]]
plot_E <- e_plot(E_BR)
# each line of the bracket's label must stay 0.02 in clear of every point and tick beside it, and inside the classifier rows
local({
  png(tempfile(fileext = ".png"), width = 3, height = 2, units = "in", res = DPI, type = "cairo"); lw <- c(str_w("Four-marker"), str_w("classifier")); dev.off()
  lh <- 0.9 * PT / 72 + 0.01; pad <- 0.02; r <- 0.375 * (1.3 * .pt + 0.15 * .stroke / 2) / 72   # one line's height (in, with slack); a point's radius (in)
  for (k in 1:2) { top <- mean(CLS_Y) + (2 - k) * lh * E_UPI[["y"]]; bot <- top - lh * E_UPI[["y"]]       # line 1 above the middle, line 2 below
    l <- E_BR$x_lab - (lw[k] + pad) * E_UPI[["x"]]; b <- bot - pad * E_UPI[["y"]]; t <- top + pad * E_UPI[["y"]]
    pts <- e_pts$v + r * E_UPI[["x"]] > l & e_pts$y + r * E_UPI[["y"]] > b & e_pts$y - r * E_UPI[["y"]] < t
    tks <- e_tick$v + TICK_W_IN * E_UPI[["x"]] > l & e_tick$y1 > b & e_tick$y0 < t
    if (any(pts) || any(tks)) stop("E: line ", k, " of the bracket label would come within ", pad, " in of a point or tick")
    if (b < RULE_Y || t > CLS_Y[2]) stop("E: the bracket label leaves the classifier rows") }
})

# key items: a study's point (filled glyph, thin black edge as ggplot's stroke 0.15) and the model-estimate tick, as long
# as the ticks are drawn
EDGE_LWD <- 0.15 * .stroke / 2
edge_item <- function(label, col, pch) list(label = label, w_sym = .14, draw = function(x, y)
  grid.points(x = unit(x + .04, "in"), y = unit(H - y - .04, "in"), pch = pch, size = unit(.085, "in"), gp = gpar(col = "black", fill = col, lwd = EDGE_LWD)))
TICK_LEN <- local({ pb <- panel_box(plot_E, PLOT_W, PLOT_H); (2 * TICK_HALF) / diff(pb$yr) * (pb$b - pb$t) })   # in
tick_item <- function(label) list(label = label, w_sym = .1, draw = function(x, y)
  grid.lines(x = unit(x + .04, "in"), y = unit(H - y - .04 + c(-1, 1) * TICK_LEN / 2, "in"), gp = gpar(col = "black", lwd = TICK_W_IN * 96, lineend = "butt")))

draw_panel <- function() {
  draw_at(plot_E, 4.78, 1.975, PLOT_W, PLOT_H)
  text_at("EC marker genes", 4.85, 1.92); text_at("FDR", 7.02, 1.92)                   # headers over the two y axes
  key_row(list(edge_item("Intestinal EC study", E_COLS[["primary EC study"]], 24), edge_item("SC-SIC study", E_COLS[["SC-SIC study"]], 21)),
          4.78, PLOT_W, 4.005)
  key_row(list(tick_item("Model estimate (from D)")), 4.78, PLOT_W, 4.165)
}
