# Figure S15C: how many genes are in all six leave-one-study-out top-20 lists (as in B) when the SC-SIC labels are
# swapped at random, as bar charts of the label sets: top, swapped between whole studies (every possible swap); bottom,
# swapped within samples (random draws). The arrow marks the count with the primary study labels. The p value counts
# the swaps reaching the observed count (within samples with the usual +1 for the observed labels). Key below.

swaps <- read_table("s1_marker_stability/b_label_swaps.tsv.gz")
obs <- swaps$genes_in_all[swaps$scope == "STUDY_BLOCK" & swaps$draw == 1]   # draw 1 = the primary study labels
nn <- c(STUDY_BLOCK = sum(swaps$scope == "STUDY_BLOCK"), WITHIN_SAMPLE = sum(swaps$scope == "WITHIN_SAMPLE"))
lv <- c(STUDY_BLOCK = sprintf("Swapped between whole studies (all %d)", nn[["STUDY_BLOCK"]]),
        WITHIN_SAMPLE = sprintf("Swapped within samples (%d draws)", nn[["WITHIN_SAMPLE"]]))
ex <- c(STUDY_BLOCK = sum(swaps$genes_in_all[swaps$scope == "STUDY_BLOCK"] >= obs),
        WITHIN_SAMPLE = sum(swaps$genes_in_all[swaps$scope == "WITHIN_SAMPLE"] >= obs))
pst <- c(STUDY_BLOCK = ex[["STUDY_BLOCK"]] / nn[["STUDY_BLOCK"]], WITHIN_SAMPLE = (1 + ex[["WITHIN_SAMPLE"]]) / (1 + nn[["WITHIN_SAMPLE"]]))
pfr <- c(STUDY_BLOCK = sprintf("%d/%d", ex[["STUDY_BLOCK"]], nn[["STUDY_BLOCK"]]),
         WITHIN_SAMPLE = sprintf("%d/%d", 1L + ex[["WITHIN_SAMPLE"]], 1L + nn[["WITHIN_SAMPLE"]]))
PFMT <- function(p) formatC(p, format = "g", digits = 2)
sub_c <- c(STUDY_BLOCK = sprintf("observed %d; %d/%d swaps ≥ %d (p = %s = %s)", obs, ex[["STUDY_BLOCK"]], nn[["STUDY_BLOCK"]], obs,
                                 pfr[["STUDY_BLOCK"]], PFMT(pst[["STUDY_BLOCK"]])),
           WITHIN_SAMPLE = sprintf("observed %d; %d/%d draws ≥ %d (p = %s = %s)", obs, ex[["WITHIN_SAMPLE"]], nn[["WITHIN_SAMPLE"]], obs,
                                   pfr[["WITHIN_SAMPLE"]], PFMT(pst[["WITHIN_SAMPLE"]])))
XM <- max(c(swaps$genes_in_all, obs)) + 2
key_down <- function(data, params, size)                          # key glyph: the drawn downward arrow
  grid::segmentsGrob(.5, .95, .5, .08, arrow = params$arrow,
                     gp = grid::gpar(col = data$colour, fill = data$colour, lwd = data$linewidth * .pt, lineend = "butt", linejoin = "mitre"))
swap_plot <- function(sc, bottom) {
  d <- swaps[swaps$scope == sc, ]
  k <- as.data.frame(table(d$genes_in_all), stringsAsFactors = FALSE); names(k) <- c("x", "n"); k$x <- as.integer(k$x)
  tl <- k[k$n < .1 * max(k$n), ]                                  # tail bars (under a tenth of the tallest) get their count
  ylo <- .32 * max(k$n)                                           # the arrow stops here, above the tail counts
  ggplot(d, aes(genes_in_all)) + geom_bar(aes(fill = "Labels swapped"), width = .8) +
    geom_text(data = tl, aes(x = x, y = n, label = n), inherit.aes = FALSE, vjust = -.3, size = PT / .pt, family = FONT, colour = "black") +
    geom_segment(data = data.frame(k = "Primary study labels"), aes(x = obs, xend = obs, y = Inf, yend = ylo, colour = k), inherit.aes = FALSE,
                 linewidth = .7, arrow = arrow(length = unit(.055, "in"), angle = 28, type = "closed"), key_glyph = key_down) +
    scale_fill_manual(values = c("Labels swapped" = "#BFC3C6"), name = NULL) + scale_colour_manual(values = c("Primary study labels" = SCC), name = NULL) +
    scale_x_continuous(limits = c(-.6, XM), breaks = scales::breaks_pretty(5), labels = minus) +
    scale_y_continuous(expand = expansion(mult = c(0, .05)), labels = minus) +
    labs(title = lv[[sc]], subtitle = sub_c[[sc]], x = if (bottom) "Genes in all six top-20 lists (as in B)" else NULL, y = "Label sets") +
    theme_rk("y", legend = "bottom") +
    theme(plot.title = element_text(size = 8, face = "bold", colour = "black", margin = margin(0, 0, 1, 0)),
          plot.subtitle = element_text(size = PT, colour = "black", margin = margin(0, 0, 3, 0)),
          axis.text.x = if (bottom) element_text(colour = "black", size = PT) else element_blank(),
          axis.ticks.x = if (bottom) element_line(colour = "black", linewidth = .3) else element_blank(),
          plot.margin = margin(if (bottom) 9 else 2, 6, 2, 2))
}
plot_C <- (swap_plot("STUDY_BLOCK", FALSE) / swap_plot("WITHIN_SAMPLE", TRUE)) +
  plot_layout(guides = "collect", axis_titles = "collect") &
  guides(fill = guide_legend(order = 1), colour = guide_legend(order = 2)) &
  theme(legend.position = "bottom", legend.box.spacing = unit(.05, "in"), legend.key.height = unit(.15, "in"))

draw_panel <- function() in_place(120, 3502, 2.95, 3.35 - .47, 1728, function() print(plot_C, newpage = FALSE))
