# Figure S7A: Figure 3C's agreement with intestinal EC (Spearman rho between each population's gene fold changes and
# intestinal EC's against other EEC) when every donor or sample counts once (colored dot; open triangle where the study
# has one donor or sample) beside the pooled-cell value of Figure 3C (small grey dot), joined by a grey segment. One row
# per population, one slot per study inside its row; the grey band spans the agreement among intestinal studies per
# donor. Key and the list of slots to the right.

dots <- read_table("s3_per_donor/dots_pooled_vs_balanced.tsv")
POPS <- c("GUT_EC", "SC_SIC", "SC_BETA", "SC_ALPHA", "INDUCED_EEC", "HIO_GRAFT_EC")
ROW <- c(GUT_EC = "Intestinal EC vs\nother EEC", SC_SIC = "SC-SIC vs other\nSC endocrine", SC_BETA = "SC-beta vs other SC\nendocrine (not SC-SIC)",
         SC_ALPHA = "SC-alpha vs other SC\nendocrine (not SC-SIC)", INDUCED_EEC = "Induced EC\n(EC vs other EEC)", HIO_GRAFT_EC = "HIO graft EC\n(EC vs other EEC)")
COL <- c(GUT_EC = "#43A5B6", SC_SIC = SIC, SC_BETA = "#4477AA", SC_ALPHA = "#117733", INDUCED_EEC = "#CCBB44", HIO_GRAFT_EC = "#6943AC")
# slot order inside a row: SC-islet studies in Figure 3A's order from the top; intestinal studies in table order
dots$slot <- NA_integer_
for (p in POPS) {
  k <- which(dots$population == p)
  dots$slot[k] <- if (p %in% c("SC_SIC", "SC_BETA", "SC_ALPHA")) match(dots$study[k], names(STUDY)) else match(dots$study[k], dots$study[k])
}
# vertical layout in slot units: rows top to bottom, each at least 2 slots tall, GAP between rows (a grey separator in it)
GAP <- 2.0; PAD <- 0.7
n_in <- table(factor(dots$population, POPS)); half <- pmax((as.numeric(n_in) - 1) / 2, 1); names(half) <- POPS
centre <- numeric(0); y <- 0
for (p in POPS) { y <- y - half[[p]]; centre[p] <- y; y <- y - half[[p]] - GAP }
dots$y <- centre[dots$population] + ((as.numeric(n_in[dots$population]) + 1) / 2 - dots$slot)
YL <- c(min(dots$y) - PAD, max(dots$y) + PAD)
top_of_row <- tapply(dots$y, factor(dots$population, POPS), max); bot_of_row <- tapply(dots$y, factor(dots$population, POPS), min)
SEP <- (bot_of_row[-length(POPS)] + top_of_row[-1]) / 2          # separators midway between neighboring rows
pts <- data.frame(population = dots$population, y = dots$y, rho = dots$balanced_min5, pooled = dots$pooled, one = dots$flag_min5 == "ONE_UNIT")
r <- range(pts$rho[pts$population == "GUT_EC"]); BAND <- data.frame(lo = r[1], hi = r[2])   # range among intestinal studies, per donor
XL <- c(-.12, .45)
SZ_DOT <- 1.7; SZ_TRI <- 1.1; ST_TRI <- .6; SZ_POOL <- 1.0; COL_POOL <- "grey50"; BAND_FILL <- "#E6E6E6"
KEY_W <- 2.30                                                     # right-hand column for the key and the slot list (inches)
plot_A <- ggplot(pts) +
  geom_rect(data = BAND, aes(xmin = lo, xmax = hi, ymin = -Inf, ymax = Inf), fill = BAND_FILL) +
  geom_hline(yintercept = SEP, colour = "grey85", linewidth = .3) +
  geom_vline(xintercept = 0, colour = REF, linewidth = .3) +
  geom_segment(aes(x = pooled, xend = rho, y = y, yend = y), colour = "grey70", linewidth = .3) +
  geom_point(aes(pooled, y), colour = COL_POOL, size = SZ_POOL, stroke = 0, shape = 16) +
  geom_point(data = pts[!pts$one, ], aes(rho, y, colour = population), size = SZ_DOT, stroke = 0, shape = 16) +
  geom_point(data = pts[pts$one, ], aes(rho, y, colour = population), size = SZ_TRI, stroke = ST_TRI, shape = 24, fill = "white") +
  scale_colour_manual(values = COL, guide = "none") +
  scale_y_continuous(breaks = unname(centre), labels = unname(ROW[names(centre)]), limits = YL, expand = c(0, 0)) +
  scale_x_continuous(limits = XL, breaks = seq(-.1, .4, .1), labels = scales::label_number(accuracy = .1, style_negative = "minus"), expand = c(0, 0)) +
  labs(x = "Agreement with intestinal EC (Spearman ρ of gene fold changes;\nmedian of 3 intestinal studies)", y = NULL) +
  theme_rk("x") +
  theme(axis.text.y = element_text(lineheight = .9), axis.title.x = element_text(lineheight = .95, margin = margin(3, 0, 0, 0)),
        plot.margin = margin(2, KEY_W * 72, 2, 2))

# key and slot list in the right-hand column, each glyph as drawn in the plot
KEY <- list(list(g = "dot", txt = "Donors or samples weighted equally"),
            list(g = "tri", txt = "Study with one donor or sample\n(only the intestinal side\naveraged per donor)"),
            list(g = "pool", txt = "Pooled cells per study (Fig. 3C)"),
            list(g = "band", txt = "Range of agreement among intestinal\nstudies, per donor or sample"))
SLOTS <- c("Slots in each row, top to bottom:",
           "SC-islet rows: Augsornworawat 2020,", "Balboa 2022, Zhu 2023, Siehler 2024,", "Du 2022, Augsornworawat 2023 (nuclei)",
           "Intestinal row: Elmentaite 2020 + 2021,", "Hickey 2023, Egozi 2021")
GLX <- .10; TXT_X <- .24; LINE_H <- .128; ROW_GAP <- .06; SLOT_GAP <- .16
draw_key <- function(x0, ytop) {
  tp <- gpar(fontfamily = FONT, fontsize = PT, col = "black", lineheight = .92)
  y <- ytop
  for (k in KEY) {
    nl <- length(strsplit(k$txt, "\n", fixed = TRUE)[[1]]); hb <- nl * LINE_H; yc <- y - hb / 2
    grid.text(k$txt, unit(x0 + TXT_X, "in"), unit(yc, "in"), just = c("left", "centre"), gp = tp)
    gx <- unit(x0 + GLX, "in"); gy <- unit(yc, "in")
    if (k$g == "dot") grid.points(gx, gy, pch = 16, gp = gpar(col = "grey25", fontsize = SZ_DOT * .pt))
    if (k$g == "tri") grid.points(gx, gy, pch = 24, gp = gpar(col = "grey25", fill = "white", fontsize = SZ_TRI * .pt + ST_TRI * .stroke / 2,
                                                               lwd = ST_TRI * .stroke / 2))
    if (k$g == "pool") grid.points(gx, gy, pch = 16, gp = gpar(col = COL_POOL, fontsize = SZ_POOL * .pt))
    if (k$g == "band") grid.rect(gx, gy, unit(.16, "in"), unit(.11, "in"), gp = gpar(fill = BAND_FILL, col = NA))
    y <- y - hb - ROW_GAP
  }
  y <- y - SLOT_GAP + ROW_GAP
  grid.text(paste(SLOTS, collapse = "\n"), unit(x0 + .04, "in"), unit(y, "in"), just = c("left", "top"), gp = tp)
}
draw_A <- function() {
  gt <- ggplotGrob(plot_A)                                        # built on the open device, which measures the text
  grid.draw(gt)
  vw <- convertWidth(unit(1, "npc"), "in", valueOnly = TRUE); vh <- convertHeight(unit(1, "npc"), "in", valueOnly = TRUE)
  pl <- gt$layout[grepl("^panel", gt$layout$name), ]              # the key starts level with the top of the plot panel
  ptop <- vh - sum(vapply(seq_len(pl$t - 1), function(i) convertHeight(gt$heights[i], "in", valueOnly = TRUE), 0))
  draw_key(vw - KEY_W + .10, ptop - .02)
}

draw_panel <- function() in_place(120, 48, 7.3, 3.20 - .47, 1638, draw_A)
