# Figure S5A: how the 20-gene classifier classifies a cell, worked through for the Elmentaite 2021 classifier (trained on
# the primary study labels of Egozi 2021 and Hickey 2023): (i) each gene's weight, the median of the two training
# studies' EC minus other terminal EEC differences, as a share of the 20 medians; (ii) the score; (iii) the cutoff,
# where the two Gaussian components fitted to the training cells' scores cross; (iv) the three conditions for classifying a cell as EC.

ELM <- "SF_ELMENTAITE2021_GUT_ATLAS"; EGO <- "SF_EGOZI2021_DEVELOPMENTAL_EEC"; HIC <- "SF_HICKEY2023_ADULT_EEC"
pred <- read_table("predictions.tsv.gz")
THR <- unique(pred$threshold[pred$study == ELM & pred$panel_size == 20])               # the classifier's score cutoff
Dt <- read_table("D_effect_heatmap.tsv", "figure2")
ORD <- unique(Dt$gene_symbol[order(Dt$panel_rank)])                                   # Fig. 2D's 20 genes, in its order
w <- read_table("weights.tsv"); w <- w[w$panel_size == 20 & w$study == ELM, ]
TRAIN_FAM <- c(EGO, HIC); ST <- study(TRAIN_FAM)
GENES <- w$gene_symbol; stopifnot(identical(GENES, ORD))
dD <- sapply(TRAIN_FAM, function(f) Dt$effect[Dt$family_id == f][match(GENES, Dt$gene_symbol[Dt$family_id == f])])
hb <- read_table("classifier_explainer_score_histogram.tsv")                          # training-cell scores, bins 0.1 wide
G <- as.list(read_table("classifier_explainer_mixture_parameters.tsv"))               # the two fitted components
brk <- c(hb$lo, hb$hi[nrow(hb)])
r <- which(GENES == "REG4")                                                           # the worked gene
NUM <- list(ego = sprintf("%.2f", dD[r, EGO]), hic = sprintf("%.2f", dD[r, HIC]), med = sprintf("%.2f", w$training_effect[r]),
            tot = sprintf("%.1f", sum(w$training_effect)), pct = sprintf("%.1f%%", 100 * w$weight[r]), cut = sprintf("%.2f", THR))
GEQ <- intToUtf8(8805); TIMES <- intToUtf8(215)
GREY <- "#6E6E6E"; WINE <- PANEL; LOWC <- "#447AA5"; UPC <- "#CC7636"; MIXC <- "#555555"   # mixture colors as in Fig. 2A
# (i) per-study differences and their median; weights as bars
di <- data.frame(gene = rep(GENES, 2), y = rep(21 - seq_along(GENES), 2), v = c(dD[, EGO], dD[, HIC]), st = rep(ST, each = 20))
dm <- data.frame(gene = GENES, y = 21 - seq_along(GENES), med = w$training_effect, wt = w$weight)
p1 <- ggplot() + geom_segment(data = dm, aes(x = pmin(dD[, EGO], dD[, HIC]), xend = pmax(dD[, EGO], dD[, HIC]), y = y, yend = y), color = "#C8C8C8", linewidth = .35) +
  geom_point(data = di, aes(v, y, shape = st), color = GREY, size = 1.25, stroke = .45) +
  geom_point(data = dm, aes(med, y), shape = 124, color = WINE, size = 2.6) +
  scale_shape_manual(values = setNames(c(1, 2), ST), guide = "none") +
  scale_y_continuous(breaks = dm$y, labels = GENES, expand = c(0, 0)) +
  scale_x_continuous(breaks = 0:5, expand = c(0, 0)) + coord_cartesian(xlim = c(-.1, 5.1), ylim = c(.4, 20.6), clip = "off") +
  labs(x = "EC minus other terminal EEC (log expression)", y = NULL) + theme_kit("x") +
  theme(axis.text.y = element_text(face = "italic"), plot.margin = margin(2, 4, 2, 2))
p1w <- ggplot(dm, aes(y = y)) + geom_rect(aes(xmin = 0, xmax = 100 * wt, ymin = y - .32, ymax = y + .32), fill = WINE) +
  annotate("text", x = 100 * w$weight[r] + .6, y = dm$y[r], label = NUM$pct, hjust = 0, size = tx, family = FONT) +
  scale_y_continuous(breaks = NULL, expand = c(0, 0)) + scale_x_continuous(breaks = c(0, 5, 10), labels = function(v) paste0(v, "%"), expand = c(0, 0)) +
  coord_cartesian(xlim = c(0, 13.5), ylim = c(.4, 20.6), clip = "off") + labs(x = "Weight", y = NULL) + theme_kit("x") +
  theme(plot.margin = margin(2, 4, 2, 2))
# (iii) training-cell scores (studies weighted equally) with the two fitted components, their sum and their crossing
BW <- .1; xs <- seq(min(brk), max(brk), length.out = 600)                           # BW: bin width
cu <- rbind(data.frame(x = xs, y = BW * G$lower_weight * dnorm(xs, G$lower_mean, G$shared_sd), k = "Lower component"),
            data.frame(x = xs, y = BW * G$upper_weight * dnorm(xs, G$upper_mean, G$shared_sd), k = "Upper component"),
            data.frame(x = xs, y = BW * (G$lower_weight * dnorm(xs, G$lower_mean, G$shared_sd) + G$upper_weight * dnorm(xs, G$upper_mean, G$shared_sd)), k = "Mixture fit"))
ymax <- max(hb$h, cu$y) * 1.08
p3 <- ggplot() + geom_rect(data = hb, aes(xmin = lo, xmax = hi, ymin = 0, ymax = h), fill = "#D9D9D9", color = "white", linewidth = .15) +
  geom_line(data = cu, aes(x, y, color = k, linetype = k), linewidth = .45) +
  geom_vline(xintercept = THR, linetype = "dashed", linewidth = .4, color = "black") +
  annotate("text", x = THR - .05, y = ymax * .97, label = paste("Cutoff", NUM$cut), hjust = 1, vjust = 1, size = tx, family = FONT) +
  scale_color_manual(values = c("Lower component" = LOWC, "Upper component" = UPC, "Mixture fit" = MIXC), guide = "none") +
  scale_linetype_manual(values = c("Lower component" = "dashed", "Upper component" = "dashed", "Mixture fit" = "solid"), guide = "none") +
  scale_y_continuous(labels = function(v) paste0(round(100 * v), "%"), expand = c(0, 0)) + scale_x_continuous(expand = c(0, 0)) +
  coord_cartesian(xlim = range(brk), ylim = c(0, ymax), clip = "off") +
  labs(x = "Score", y = "Cells (share per bin)") + theme_kit("y") +
  theme(axis.line.y = element_line(color = "black", linewidth = .3), axis.ticks.y = element_line(color = "black", linewidth = .3), plot.margin = margin(2, 4, 2, 2))
G1 <- ggplotGrob(p1); G1w <- ggplotGrob(p1w); G3 <- ggplotGrob(p3)
G1w$heights <- G1$heights                                                             # weight bars line up with the gene rows

W_A <- 7.3; H_A <- 4.98                                                               # panel width and height (in)
N2 <- "Each point: Egozi 2021 or Hickey 2023, EC minus other terminal EEC in donor-averaged mean expression."
N3 <- "Of these steps, only this one uses labels (primary study labels of Egozi 2021 and Hickey 2023);"
N4 <- "the 20 genes were ranked beforehand (Fig. 2D)."
Z1 <- "z-score = (expression - mean) / SD, with mean and SD over all"
Z2 <- "Egozi 2021 and Hickey 2023 cells, each study weighted equally,"
Z3 <- "whatever their label"
C1 <- paste("Egozi 2021 and Hickey 2023 cells with", GEQ, "5 of 20 genes detected,")
C2 <- "each study weighted equally; fitted without labels"
draw_a <- function() {
  vh <- convertHeight(unit(1, "npc"), "in", valueOnly = TRUE)
  at <- function(gr, x, y, w, h) { pushViewport(viewport(x = unit(x, "in"), y = unit(vh - y, "in"), width = unit(w, "in"), height = unit(h, "in"), just = c("left", "top"))); grid.draw(gr); popViewport() }
  tt <- function(s, x, y, bold = FALSE, italic = FALSE, just = c("left", "top"), col = "black")
    grid.text(s, x = unit(x, "in"), y = unit(vh - y, "in"), just = just, gp = gpar(fontfamily = FONT, fontsize = PT, fontface = if (bold) "bold" else if (italic) "italic" else "plain", col = col))
  ex <- function(e, x, y, just = c("left", "top")) grid.text(e, x = unit(x, "in"), y = unit(vh - y, "in"), just = just, gp = gpar(fontfamily = FONT, fontsize = PT, col = "black"))
  # ---- (i), left column ----
  tt("(i) Weights: median EC minus other terminal EEC difference", .04, .02, bold = TRUE)
  at(G1, .02, .20, 2.55, 2.95); at(G1w, 2.70, .20, .90, 2.95)
  ky <- 3.25; kx <- .30
  for (i in 1:2) { grid.points(unit(kx + .05, "in"), unit(vh - ky - .055, "in"), pch = c(1, 2)[i], size = unit(1, "char"), gp = gpar(col = GREY, fontsize = 1.25 * .pt, lwd = .45 * .pt))
    tt(ST[i], kx + .13, ky); kx <- kx + .13 + str_in(ST[i]) + .22 }
  grid.lines(unit(rep(kx + .05, 2), "in"), unit(vh - ky - c(.0, .11), "in"), gp = gpar(col = WINE, lwd = 1.6)); tt("Median", kx + .13, ky)
  ex(bquote(italic("REG4")*": median of "*.(NUM$ego)*" and "*.(NUM$hic)*" = "*.(NUM$med)*"; "*.(NUM$med)*" / "*.(NUM$tot)*" (sum of the 20 medians) = "*.(NUM$pct)), .04, 3.42)
  tt(N2, .04, 3.56); tt(N3, .04, 3.70); tt(N4, .04, 3.84)
  # ---- (ii) and (iii), right column ----
  X2 <- 3.75
  tt("(ii) Score: weighted mean of the 20 z-scores", X2, .02, bold = TRUE)
  tt(Z1, X2 + .1, .20); tt(Z2, X2 + .1, .335); tt(Z3, X2 + .1, .47)
  tt(paste0("score = sum over the 20 genes of (weight ", TIMES, " z-score)"), X2 + .1, .61)
  tt("(iii) Cutoff: where the two fitted score components cross", X2, .78, bold = TRUE)
  tt(C1, X2 + .1, .96); tt(C2, X2 + .1, 1.095)
  at(G3, X2, 1.25, W_A - X2 - .04, 2.05)
  kl <- c("Mixture fit", "Lower component", "Upper component"); kc <- c(MIXC, LOWC, UPC); kt <- c(1, 2, 2); kx <- X2 + 2.25   # key inside the plot, top right
  for (i in 1:3) { ky <- 1.33 + (i - 1) * .15; grid.lines(unit(kx + c(0, .2), "in"), unit(rep(vh - ky - .055, 2), "in"), gp = gpar(col = kc[i], lty = kt[i], lwd = .45 * .pt)); tt(kl[i], kx + .26, ky) }
  # ---- (iv), full width: three condition boxes joined by "and", an arrow, the result ----
  Y4 <- 4.01
  tt("(iv) Call: a cell is called EC when all three hold", .04, Y4, bold = TRUE)
  bx <- c(paste0("Score ", GEQ, " cutoff (", NUM$cut, ")"), paste(GEQ, "5 of the 20 genes detected"), "Score-positive cells enriched\namong its 25 nearest neighbors")
  bwid <- vapply(strsplit(bx, "\n", fixed = TRUE), str_in, 0) + .16; gap <- .34; ec <- "EC call"; ew <- str_in(ec, bold = TRUE) + .2
  tot <- sum(bwid) + 3 * gap + ew
  x0 <- (W_A - tot) / 2; yb <- Y4 + .17; hb2 <- .30
  for (i in 1:3) {
    grid.roundrect(unit(x0, "in"), unit(vh - yb, "in"), unit(bwid[i], "in"), unit(hb2, "in"), r = unit(.04, "in"), just = c("left", "top"), gp = gpar(col = "#8C8C8C", fill = "white", lwd = .6))
    tt(bx[i], x0 + bwid[i] / 2, yb + hb2 / 2, just = c("centre", "centre"))
    x0 <- x0 + bwid[i]
    tt(if (i < 3) "and" else "", x0 + gap / 2, yb + hb2 / 2, just = c("centre", "centre"))
    if (i == 3) grid.lines(unit(x0 + c(.05, gap - .05), "in"), unit(rep(vh - yb - hb2 / 2, 2), "in"), arrow = arrow(length = unit(.05, "in"), type = "closed"), gp = gpar(col = "black", fill = "black", lwd = .7))
    x0 <- x0 + gap }
  grid.roundrect(unit(x0, "in"), unit(vh - yb, "in"), unit(ew, "in"), unit(hb2, "in"), r = unit(.04, "in"), just = c("left", "top"), gp = gpar(col = WINE, fill = "white", lwd = 1.2))
  tt(ec, x0 + ew / 2, yb + hb2 / 2, bold = TRUE, just = c("centre", "centre"))
}

draw_panel <- function() place(draw_a, 120, 48, 4380, 2742, W_A, H_A - (.30 + .17))
