# Figure S1B: how the four-marker EC classifier classified the cells of each fitting group and of all of them:
# stacked shares of high-confidence EC, low-confidence EC and cells below the cutoff, with the cell counts in columns on
# the right. Groups in Figure S1A's order (highest share at or above the cutoff first). Key below.

x <- read_table("fitting_outcome_by_arm.tsv")
arm <- x[x$scope == "arm", ]; tot <- x[x$scope == "all_arms", ]
arm <- arm[order(-arm$local_consistency_entering / arm$cells), ]          # S1A's order: share at or above the cutoff
x <- rbind(arm, tot)
k <- x$scope == "arm"
x$study <- "All fitting groups"; x$study[k] <- study(x$arm_id[k])
x$n_lab <- fmt(x$cells)
x$line2 <- NA; x$line2[k] <- sprintf("%s (%s)", study_type(x$arm_id[k]), x$n_lab[k])
x$single <- sprintf("%s (%s)", x$study, x$n_lab)
x$y <- rev(seq_len(nrow(x))) + ifelse(k, .35, 0)                         # the total row set apart below the groups
CL <- names(CALLS)
st <- do.call(rbind, lapply(seq_len(nrow(x)), function(i) {
  n <- c(x$high_confidence_ec[i], x$low_confidence_ec[i], x$insufficient_ec_evidence[i])
  data.frame(y = x$y[i], cl = CL, n = n, x1 = cumsum(n) / x$cells[i], x0 = (cumsum(n) - n) / x$cells[i])
}))
COLS <- c(HIGH_CONFIDENCE_EC = "High-confidence", LOW_CONFIDENCE_EC = "Low-confidence", INSUFFICIENT_EC_EVIDENCE = "Below cutoff")
cnt <- data.frame(y = rep(x$y, 3), cl = rep(CL, each = nrow(x)), n = c(x$high_confidence_ec, x$low_confidence_ec, x$insufficient_ec_evidence))
cnt$lab <- fmt(cnt$n)
SEP_Y <- 1.67                                                            # dotted line above the total row
KEY_B <- c("High-confidence EC (classified as EC): at or above the cutoff, neighbors enriched for score-positive cells",
           "Low-confidence EC (not classified as EC): at or above the cutoff, neighbors not enriched",
           "Insufficient EC evidence: below the cutoff")
W_B <- 7.3; H_B <- 2.90; KEY_H <- .34; HEAD_H <- .05                     # panel width and height, key band, header (in)
draw_b <- function() {
  vh <- convertHeight(unit(1, "npc"), "in", valueOnly = TRUE)
  lab_w <- max(str_in(x$study[1:5], bold = TRUE), str_in(c(x$line2[1:5], x$single[6]))) + .10
  cw <- vapply(CL, function(c) max(str_in(COLS[[c]]), str_in(cnt$lab[cnt$cl == c])), 0)
  right <- cumsum(cw + .16)                              # right edge of each count column
  CW <- max(right) + .02
  ytop <- max(x$y) + .5
  pa <- ggplot(st) + geom_hline(yintercept = SEP_Y, color = REF, linewidth = .25, linetype = "dotted") +
    geom_rect(aes(xmin = x0, xmax = x1, ymin = y - .3, ymax = y + .3, fill = cl)) +
    annotate("text", x = 0, y = x$y[1:5] + .21, label = paste0(x$study[1:5], "  "), hjust = 1, vjust = .5, fontface = "bold",
             size = tx, family = FONT, color = "black") +
    annotate("text", x = 0, y = x$y[1:5] - .21, label = paste0(x$line2[1:5], "  "), hjust = 1, vjust = .5,
             size = tx, family = FONT, color = "black") +
    annotate("text", x = 0, y = x$y[6], label = paste0(x$single[6], "  "), hjust = 1, vjust = .5,
             size = tx, family = FONT, color = "black") +
    scale_fill_manual(values = CALLS) +
    scale_y_continuous(breaks = x$y, labels = NULL, limits = c(.5, ytop), expand = c(0, 0)) +
    scale_x_continuous(limits = c(0, 1), breaks = c(0, .25, .5, .75, 1), labels = c("0%", "25%", "50%", "75%", "100%"), expand = c(0, 0)) +
    coord_cartesian(clip = "off") +
    labs(x = "Share of the group's cells", y = NULL) + theme_kit("x") +
    theme(axis.text.y = element_blank(), plot.margin = margin(HEAD_H * 72, 9, 2, 72 * lab_w), plot.background = element_blank())
  pc <- ggplot(cnt) + geom_hline(yintercept = SEP_Y, color = REF, linewidth = .25, linetype = "dotted") +
    geom_text(aes(x = right[match(cl, CL)], y = y, label = lab), hjust = 1, vjust = .5, size = tx, family = FONT, color = "black") +
    annotate("text", x = right, y = ytop, label = unname(COLS), hjust = 1, vjust = -.45, size = tx, family = FONT, color = "black") +
    scale_y_continuous(limits = c(.5, ytop), expand = c(0, 0)) + scale_x_continuous(limits = c(0, CW), expand = c(0, 0)) +
    coord_cartesian(clip = "off") + theme_void() +
    theme(plot.margin = margin(HEAD_H * 72, 0, 2, 0), plot.background = element_blank())
  p <- pa + pc + plot_layout(widths = unit(c(1, CW), c("null", "in")))
  pushViewport(viewport(y = unit(KEY_H, "in"), height = unit(vh - KEY_H, "in"), just = "bottom"))
  print(p, newpage = FALSE); popViewport()
  # key: high-confidence alone, then low-confidence and insufficient side by side
  k3 <- .04 + .17 + str_in(KEY_B[2]) + .25
  ky <- c(KEY_H - .09, KEY_H - .245, KEY_H - .245); kxx <- c(.04, .04, k3)
  for (i in 1:3) {
    grid.rect(x = unit(kxx[i] + .05, "in"), y = unit(ky[i], "in"), width = unit(.095, "in"), height = unit(.095, "in"),
              gp = gpar(fill = unname(CALLS[i]), col = NA))
    grid.text(KEY_B[i], x = unit(kxx[i] + .17, "in"), y = unit(ky[i], "in"), just = "left",
              gp = gpar(fontfamily = FONT, fontsize = PT, col = "black"))
  }
}

draw_panel <- function() place(draw_b, 120, 2304, 4380, 1494, W_B, H_B - (.30 + .17))
