# Shared by S5A, S5B and S5_key (sourced by each): the two panels of Figure S5 and their place on the page.
#
# For each gene group (G1, G2, G3, as Figures 2D and 2E) and each drawn term: -log10 adjusted p with the group alone (open
# ring; a crossed ring at 0 when no gene of the group is in the term) and with the group widened by its mixed genes
# (filled dot; mixed genes: 2-fold up in EC in two of the three intestinal studies), joined by a line; the term name left
# of 0 (bold: also drawn in Figure 2E); right of the plot, the mixed genes in the term. A: pathway terms; B: cell-type
# signatures. The drawn terms, their names and the bold flags are listed in terms_s5.tsv: for each group, the terms of
# Figure 2E and the terms called with the mixed genes.

BLOCKS <- c("G1", "G2", "G3")                                          # top to bottom
GROUP_COL <- c(G1 = "#43A5B6", G2 = "#A89A2E", G3 = "#B4553F")
GROUP_NAME <- c(G1 = "G1: EC genes low in SC-SICs", G2 = "G2: shared EC genes", G3 = "G3: separate SC-SIC program")
WIDE <- c(G1 = "G1_PLUS_MIXED_G1_LIKE", G2 = "G2_PLUS_MIXED_G2_LIKE", G3 = "G3_PLUS_MIXED_G3_IF_NOT_SHARED")   # each group + its mixed genes
GRID_COL <- "#E4E4E4"

LT <- read_table("lists.tsv")
N_GROUP <- vapply(BLOCKS, function(g) sum(LT$list == g), 0L)
N_ADDED <- vapply(BLOCKS, function(g) length(setdiff(LT$gene[LT$list == WIDE[[g]]], LT$gene[LT$list == g])), 0L)
BLOCK_LAB <- setNames(sprintf("%s (%d + %d mixed genes)", GROUP_NAME, N_GROUP, N_ADDED), BLOCKS)

# The drawn rows of one panel: each term's adjusted p for the group alone (NA: no group gene in the term) and with its
# mixed genes, and the mixed genes in the term. all: one row per gene list and term; key: each row's term as
# terms_s5.tsv names it.
panel_terms <- function(panel, all, key) {
  D <- read_placement("terms_s5.tsv"); D <- D[D$panel == panel, ]
  ia <- match(paste(D$block, D$term), paste(all$list, key))
  ib <- match(paste(WIDE[D$block], D$term), paste(all$list, key))
  stopifnot(!anyNA(ib))
  out <- data.frame(block = D$block, term = D$name, p_without = all$p_adj[ia], p_with = all$p_adj[ib], in2e = D$in_2e,
                    mixed = vapply(seq_along(ib), function(i) paste(sort(setdiff(strsplit(all$genes[ib[i]], "/", fixed = TRUE)[[1]],
                                                                         LT$gene[LT$list == D$block[i]])), collapse = ", "), ""),
                    stringsAsFactors = FALSE)
  out[order(match(out$block, BLOCKS), seq_len(nrow(out))), ]
}
pa <- read_table("pathways_all_terms.tsv.gz")
A <- panel_terms("A", pa, paste0(pa$collection, ": ", pa$term))
ca <- read_table("celltypes_all_terms.tsv")
B <- panel_terms("B", ca, ca$term)

# ---- layout (inches) ----
X0 <- 0.10; PW <- 7.30                                      # plot box left and width
meas <- function(s, face = "plain") { png(tempfile(fileext = ".png"), width = 4, height = 2, units = "in", res = 600, type = "cairo"); on.exit(dev.off())
  vapply(s, function(z) convertWidth(grobWidth(textGrob(z, gp = gpar(fontfamily = FONT, fontsize = PT, fontface = face))), "in", valueOnly = TRUE), 0) }
LAB_GAP <- 0.10; GENE_GAP <- 0.14; GENE_W <- 2.55
TERM_W <- max(meas(c(A$term[A$in2e], B$term[B$in2e]), "bold"), meas(c(A$term, B$term))) + LAB_GAP + 0.04   # widest name fits left of 0
P <- 0.155; LINE2 <- 0.118; HDR <- 0.19; BGAP <- 0.04; COLHDR <- 0.16; BOT <- 0.05
XMAX <- 10; BREAKS <- seq(0, XMAX, 2); THR <- -log10(0.05)
KINDS <- c("Group alone", "Group with its mixed genes", "No group gene in the term")
# points (hairline black outline): open = colored ring over a black ring just outside it; filled = colored fill with a
# black edge; crossed = colored, no black (it sits at 0 and carries no value)
RING <- list(size = 1.8, stroke = .7, under = 2.75, under_stroke = .15); DOT <- list(size = 2.35, stroke = .15)

wrap_genes <- function(s, width) {                          # greedy wrap after a comma, by measured italic width
  if (!nzchar(s)) return("no mixed gene in the term")
  z <- strsplit(s, ", ", fixed = TRUE)[[1]]; lines <- character(); cur <- z[1]
  for (k in seq_along(z)[-1]) { try <- paste0(cur, ", ", z[k])
    if (meas(if (k == length(z)) try else paste0(try, ","), "italic") <= width) cur <- try else { lines <- c(lines, paste0(cur, ",")); cur <- z[k] } }
  paste(c(lines, cur), collapse = "\n") }

panel_data <- function(D) {
  D$mixed_lab <- vapply(D$mixed, wrap_genes, "", width = GENE_W - GENE_GAP - 0.04)
  D$gface <- ifelse(nzchar(D$mixed), "italic", "plain"); D$nl <- lengths(strsplit(D$mixed_lab, "\n", fixed = TRUE))
  D$tface <- ifelse(D$in2e, "bold", "plain"); D$block <- factor(D$block, BLOCKS)
  y <- -COLHDR; HB <- NULL; D$yc <- NA_real_                # y in inches below the plot top
  for (b in BLOCKS) {
    if (b != BLOCKS[1]) y <- y - BGAP
    top <- y - HDR; HB <- rbind(HB, data.frame(block = b, hy = y - HDR / 2, top = top, bottom = NA_real_)); y <- top
    for (i in which(D$block == b)) { h <- P + (D$nl[i] - 1) * LINE2; D$yc[i] <- y - h / 2; y <- y - h }
    HB$bottom[HB$block == b] <- y
  }
  D$v0 <- ifelse(is.na(D$p_without), 0, -log10(D$p_without)); D$v1 <- -log10(D$p_with)
  D$kind0 <- ifelse(is.na(D$p_without), KINDS[3], KINDS[1])
  list(D = D, HB = HB, YMIN = y - BOT)
}
build <- function(pd, upi, left_in) {                       # upi: x units per inch; text at fixed inch offsets from the axis
  D <- pd$D; HB <- pd$HB; YMIN <- pd$YMIN; xl <- function(inch) inch * upi; x_canvas <- -left_in * upi
  GR <- merge(HB, data.frame(x = BREAKS)); O <- D[D$kind0 == KINDS[1], ]; X <- D[D$kind0 == KINDS[3], ]
  ggplot() +
    geom_segment(data = GR, aes(x = x, xend = x, y = top, yend = bottom), colour = GRID_COL, linewidth = .2) +
    geom_segment(data = HB, aes(x = THR, xend = THR, y = top, yend = bottom), colour = "grey45", linewidth = .3, linetype = "22") +
    geom_segment(data = D, aes(x = v0, xend = v1, y = yc, yend = yc), colour = "grey72", linewidth = .35) +
    geom_point(data = O, aes(x = v0, y = yc), shape = 21, fill = "white", colour = "black", size = RING$under - RING$under_stroke, stroke = RING$under_stroke) +   # white inside hides the connecting line
    geom_point(data = O, aes(x = v0, y = yc, colour = block), shape = 1, size = RING$size, stroke = RING$stroke) +
    geom_point(data = X, aes(x = v0, y = yc, colour = block), shape = 13, size = RING$size, stroke = RING$stroke) +
    geom_point(data = D, aes(x = v1, y = yc, fill = block), shape = 21, colour = "black", size = DOT$size, stroke = DOT$stroke) +
    geom_text(data = D, aes(x = -xl(LAB_GAP), y = yc, label = term, fontface = tface), hjust = 1, size = MM, family = FONT, colour = "black") +
    geom_text(data = D, aes(x = XMAX + xl(GENE_GAP), y = yc, label = mixed_lab, fontface = gface), hjust = 0, vjust = .5, lineheight = .95,
              size = MM, family = FONT, colour = "black") +
    annotate("rect", xmin = x_canvas + xl(0.02), xmax = x_canvas + xl(0.105), ymin = HB$hy - .0425, ymax = HB$hy + .0425,
             fill = GROUP_COL[HB$block], colour = NA) +
    annotate("text", x = x_canvas + xl(0.15), y = HB$hy, label = BLOCK_LAB[HB$block], hjust = 0, size = MM, family = FONT, fontface = "bold", colour = "black") +
    annotate("text", x = XMAX + xl(GENE_GAP), y = -COLHDR / 2, label = "Mixed genes in the term", hjust = 0, size = MM, family = FONT, colour = "black") +
    annotate("segment", x = 0, xend = XMAX, y = YMIN, yend = YMIN, colour = "black", linewidth = .35) +
    scale_colour_manual(values = GROUP_COL, guide = "none") + scale_fill_manual(values = GROUP_COL, guide = "none") +
    scale_x_continuous(breaks = BREAKS) +
    coord_cartesian(xlim = c(0, XMAX), ylim = c(YMIN, 0), expand = FALSE, clip = "off") +
    labs(x = expression("−log"[10] ~ "adjusted" ~ italic(p)), y = NULL) +
    theme_s5() + theme(panel.grid.major = element_blank(), panel.grid.minor = element_blank(), axis.line.x = element_blank(), axis.text.y = element_blank(), axis.ticks.y = element_blank(),
                       plot.margin = unit(c(0, GENE_W, .02, TERM_W), "in"))
}
geom_in <- function(g, h) {                                 # fixed parts of the plot around its panel (in)
  png(tempfile(fileext = ".png"), width = PW, height = h, units = "in", res = 600, type = "cairo"); on.exit(dev.off())
  gt <- ggplotGrob(g); pushViewport(viewport(width = unit(PW, "in"), height = unit(h, "in")))
  wd <- convertWidth(gt$widths, "in", valueOnly = TRUE); ht <- convertHeight(gt$heights, "in", valueOnly = TRUE); popViewport()
  nw <- unitType(gt$widths) == "null"; nh <- unitType(gt$heights) == "null"; pl <- gt$layout[gt$layout$name == "panel", ]
  list(left = sum(wd[seq_len(pl$l - 1)]), fixed_w = sum(wd[!nw]), fixed_h = sum(ht[!nh])) }
make_panel <- function(D) {
  pd <- panel_data(D); G0 <- geom_in(build(pd, XMAX / 2, TERM_W), 4)
  g <- build(pd, XMAX / (PW - G0$fixed_w), G0$left)
  list(g = g, h = geom_in(g, 4)$fixed_h - pd$YMIN, pd = pd)
}
PA <- make_panel(A); PB <- make_panel(B)

# ---- page: A from the top, B under it, the shared key under B ----
TOP <- 0.12; GAP_AB <- 0.22
yA <- TOP + 0.02; yB <- yA + PA$h + GAP_AB; yK <- yB + PB$h + 0.10
