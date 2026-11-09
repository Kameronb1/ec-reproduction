# Figure S1C: the four-marker EC classifier refitted with one fitting group left out (and with both Elmentaite
# 2020 groups left out together), then all fitting cells reclassified: cells that lost (left) or gained (right) EC
# classification against the full fit. Row labels on the left, key below.

r <- read_table("cutoff_stability_summary.tsv")
r <- r[r$drawn %in% c(TRUE, "TRUE"), ]
arms <- r[r$design == "Leave one arm out", ]; fam <- r[r$design != "Leave one arm out", ]
arms <- arms[match(c("E-MTAB-8901_paediatric", "E-MTAB-10187_intestinal", "GSE125970", "GSE185224", "E-MTAB-8901_fetal"),
                   arms$held_out), ]                                       # S1A's order
r <- rbind(arms, fam)
r$y <- c(rev(seq_len(nrow(arms))) + 1.1, 1)                                # the study-level row set apart at the bottom
r$study <- c(study(arms$held_out), study("SF_ELMENTAITE2020_DEVGUT"))
r$line2 <- c(sprintf("%s (%s)", study_type(arms$held_out), fmt(arms$held_cells)), sprintf("both groups (%s)", fmt(fam$held_cells)))
HEADER <- c(labels = "Left out (cells)", bars = "Cells that changed EC classification")
SEP_C <- 1.55; LOST <- "#BE654C"; GAINED <- "#397CA5"                      # dotted line above the study-level row
W_C <- 7.3; H_C <- 2.40; TOP_IN <- .13; KEY_BAND <- .33                    # panel width and height, header, key band (in)
KEY_C <- c("Lost (EC in full fit, not refit)", "Gained (EC in refit, not full fit)")
draw_c <- function() {
  lab_w <- max(str_in(r$study, bold = TRUE), str_in(c(r$line2, HEADER[["labels"]]))) + .10
  CW1 <- .02; XL <- c(-44, 12)
  ytop <- max(r$y) + .5
  sep <- geom_hline(yintercept = SEP_C, color = REF, linewidth = .25, linetype = "dotted")
  yax <- scale_y_continuous(limits = c(.5, ytop), expand = c(0, 0))
  pt <- ggplot() + sep +                                 # row labels and their header
    annotate("text", x = 0, y = r$y + .22, label = paste0(r$study, "  "), hjust = 1, vjust = .5, fontface = "bold",
             size = tx, family = FONT, color = "black") +
    annotate("text", x = 0, y = r$y - .22, label = paste0(r$line2, "  "), hjust = 1, vjust = .5, size = tx, family = FONT, color = "black") +
    annotate("text", x = 0, y = ytop, label = paste0(HEADER[["labels"]], "  "), hjust = 1, vjust = -.5, size = tx, family = FONT, color = "black") +
    yax + scale_x_continuous(limits = c(0, CW1), expand = c(0, 0)) +
    coord_cartesian(clip = "off") + theme_void() +
    theme(plot.margin = margin(TOP_IN * 72, 0, 2, 72 * lab_w), plot.background = element_blank())
  bars <- rbind(data.frame(y = r$y, v = -r$lost, n = r$lost, side = "lost"), data.frame(y = r$y, v = r$gained, n = r$gained, side = "gained"))
  pb <- ggplot(bars, aes(xmin = 0, xmax = v, ymin = y - .28, ymax = y + .28, fill = side)) + sep +
    geom_rect() + geom_vline(xintercept = 0, linewidth = .3) +
    geom_text(data = bars[bars$n > 0, ], aes(x = v, y = y, label = n, hjust = ifelse(v < 0, 1.25, -.25)),
              inherit.aes = FALSE, size = tx, family = FONT) +
    scale_fill_manual(values = c(lost = LOST, gained = GAINED)) + yax +
    scale_x_continuous(limits = XL, breaks = c(-30, -20, -10, 0, 10), labels = function(v) abs(v), expand = c(0, 0)) +
    annotate("text", x = mean(XL), y = ytop, label = HEADER[["bars"]], hjust = .5, vjust = -.5, size = tx, family = FONT, color = "black") +
    coord_cartesian(clip = "off") +
    labs(x = NULL, y = NULL) + theme_kit("x") +
    theme(axis.text.y = element_blank(), plot.margin = margin(TOP_IN * 72, 6, 2, 8), plot.background = element_blank())
  # the band under the bars is made exactly KEY_BAND tall, so the key sits in it at measured positions
  build <- function(extra_pt) {
    p <- pt + (pb + theme(plot.margin = margin(TOP_IN * 72, 6, 2 + extra_pt, 8))) +
      plot_layout(widths = unit(c(CW1, 1), c("in", "null"))) + plot_annotation(theme = theme(plot.margin = margin(0, 0, 0, 0)))
    g <- patchworkGrob(p); pr <- unique(g$layout$b[g$layout$name %in% c("panel-1", "panel-2")])
    hin <- vapply(seq_along(g$heights), function(i) in_h(g$heights[i]), 0)
    list(g = g, below = sum(hin[(pr + 1):length(hin)]))
  }
  b0 <- build(0); b1 <- build(max(0, (KEY_BAND - b0$below) * 72))
  grid.draw(b1$g)
  yk <- b1$below - c(.085, .23)                          # key line centers above the content's foot
  for (i in 1:2) {
    grid.rect(x = unit(.09, "in"), y = unit(yk[i], "in"), width = unit(.095, "in"), height = unit(.095, "in"),
              gp = gpar(fill = c(LOST, GAINED)[i], col = NA))
    grid.text(KEY_C[i], x = unit(.21, "in"), y = unit(yk[i], "in"), just = "left", gp = gpar(fontfamily = FONT, fontsize = PT, col = "black"))
  }
}

draw_panel <- function() place(draw_c, 114, 3898, 4380, 1194, W_C, H_C - (.30 + .17))
