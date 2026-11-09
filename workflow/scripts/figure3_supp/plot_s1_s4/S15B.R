# Figure S15B: with each SC-islet study left out of the ranking in turn (columns), is each gene still in the top 20
# (dark), at rank 21-30 (light) or further down? Rows: A's 20 genes in A's rank order, then the other genes that are in
# at least two of the six top-20 lists; right, the number of top-20 lists holding each gene (bold: all six). Key below.

top20 <- read_table("s1_marker_limma/top20.tsv")
D20 <- top20$gene[order(top20$rank)]                             # A's 20 genes, ranked
STATES <- c("Top 20", "Rank 21–30", "Outside top 30")
fold <- read_table("s1_marker_stability/a_fold_ranks.tsv.gz")
fold$state <- ifelse(!is.na(fold$rank) & fold$rank <= 20, STATES[1], ifelse(!is.na(fold$rank) & fold$rank <= 30, STATES[2], STATES[3]))
cnt <- tapply(fold$state == STATES[1], fold$gene, sum)           # top-20 lists holding each gene
others <- setdiff(unique(fold$gene), D20)
others <- others[cnt[others] >= 2]                               # genes in one list only are not drawn
others <- others[order(-cnt[others], others)]
stopifnot(all(D20 %in% fold$gene))
fold <- fold[fold$gene %in% c(D20, others), ]

ROW <- c(setNames(seq_along(D20), D20), setNames(length(D20) + 1 + seq_along(others), others))   # row 21 holds "Other genes"
NR <- max(ROW)
XC <- XT + NS * CW + .20                                          # center of the count column
WB <- XC + .10 + .12                                              # plot width and height (inches)
HB <- HEAD + NR * P + .06 + LBLH + .14 + .14 + .04
y0 <- HB - HEAD; yb <- y0 - NR * P
ry <- function(r) y0 - (r - .5) * P
fold$x <- XT + (match(fold$left_out, names(STUDY)) - .5) * CW; fold$y <- ry(ROW[fold$gene])
fold$state <- factor(fold$state, levels = STATES)
lab <- data.frame(gene = names(ROW), y = ry(ROW), n = cnt[names(ROW)])
lab$txt <- paste0(lab$n, "/", NS); lab$face <- ifelse(lab$n == NS, "bold", "plain")
cols <- setNames(c(SCC, SCC_LIGHT, "#F1F1F1"), STATES)
KT <- yb - .06 - LBLH - .14; KS <- .12; KG <- .05; KSEP <- .18    # key: three chips in one row, centered
kw <- str_w_at(STATES, 100); kx <- WB / 2 - (sum(KS + KG + kw) + KSEP * 2) / 2 + c(0, cumsum(KS + KG + kw + KSEP))[1:3]
key <- data.frame(state = factor(STATES, levels = STATES), x = kx + KS / 2, tx = kx + KS + KG, y = KT - KS / 2)

plot_B <- ggplot() +
  geom_tile(data = fold, aes(x, y, fill = state), width = CW * .95, height = P * .9, colour = "white", linewidth = .4) +
  geom_tile(data = key[key$state != STATES[3], ], aes(x, y, fill = state), width = KS, height = KS, colour = NA) +
  geom_tile(data = key[key$state == STATES[3], ], aes(x, y, fill = state), width = KS, height = KS, colour = "grey80", linewidth = .3) +
  scale_fill_manual(values = cols, breaks = STATES, guide = "none") +
  tx(data = data.frame(x = XG, y = ry(ROW), lab = names(ROW)), aes(x, y, label = lab), hjust = 1, vjust = .5, fontface = "italic") +
  tx(data = lab[lab$face == "plain", ], aes(x = XC, y = y, label = txt), hjust = .5, vjust = .5) +
  tx(data = lab[lab$face == "bold", ], aes(x = XC, y = y, label = txt), hjust = .5, vjust = .5, fontface = "bold") +
  tx(data = data.frame(x = XT + (seq_len(NS) - .5) * CW, y = yb - .06, lab = unname(STUDY)), aes(x, y, label = lab), angle = ANG, hjust = 1, vjust = .5) +
  tx(data = data.frame(x = c(XL, XT + NS * CW / 2, XC), y = y0 + .05, lab = c("A's 20 genes", "Study left out", "Lists"), h = c(0, .5, .5)),
     aes(x, y, label = lab, hjust = h), vjust = 0) +
  annotate("segment", x = XL, xend = XC + .14, y = y0 - length(D20) * P, yend = y0 - length(D20) * P, linewidth = .35, colour = "black") +
  tx(data = data.frame(x = XL, y = ry(length(D20) + 1), lab = "Other genes in ≥ 2 lists"), aes(x, y, label = lab), hjust = 0, vjust = .5) +
  tx(data = key, aes(x = tx, y = y, label = state), hjust = 0, vjust = .5) +
  scale_x_continuous(limits = c(0, WB), expand = c(0, 0)) + scale_y_continuous(limits = c(0, HB), expand = c(0, 0)) +
  coord_cartesian(clip = "off") + theme_void(base_family = FONT) + theme(plot.margin = margin(0, 0, 0, 0), legend.position = "none")

draw_panel <- function() in_place(2220, 48, WB, HB, 3394, function() print(plot_B, newpage = FALSE))
