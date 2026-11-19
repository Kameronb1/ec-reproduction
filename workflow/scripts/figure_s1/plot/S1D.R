# Supplementary Figure S1D: the genes of Figure 2C (FDR < 0.05, at least 2-fold) before (x) and after (y) each check,
# signed toward the group the gene is higher in: after ambient-RNA correction, and after removing cells with many reads
# of other cell types' marker genes. Grey band: below 2-fold after the check. Named: the genes below 2-fold after a
# check (bold) and the two largest drops in each plot.

GD <- read_table("genes_before_after.tsv.gz")
GD <- GD[GD$core %in% c(TRUE, "TRUE"), ]
SIDE <- c(NATIVE_EC_HIGHER = "Higher in intestinal EC", SC_SIC_HIGHER = "Higher in SC-SIC")
CHECK <- c(soupx = "after ambient-RNA correction", doublets = "after removing cells with many\nreads of other cell types")
SIDECOL <- c(NATIVE_EC_HIGHER = COL[["gut_ec"]], SC_SIC_HIGHER = COL[["sc_sic"]])
LD <- do.call(rbind, lapply(names(SIDE), function(sd) do.call(rbind, lapply(names(CHECK), function(ck) {
  g <- GD[GD$call == sd, ]; s <- if (sd == "SC_SIC_HIGHER") 1 else -1
  data.frame(side = sd, check = ck, gene = g$gene, before = s * g$logFC, after = s * g[[paste0("logFC_", ck)]]) }))))
# plot headers: side, check, and how many of the plotted genes stay at 2-fold or more
STRIP_D <- do.call(rbind, lapply(split(LD, list(LD$side, LD$check)), function(d) data.frame(side = d$side[1], check = d$check[1],
  kept = sum(d$after >= 1), n = nrow(d), lab = sprintf("%s,\n%s\n%s of %s still 2-fold or more", SIDE[[d$side[1]]], CHECK[[d$check[1]]], fmt(sum(d$after >= 1)), fmt(nrow(d))))))
STRIP_D <- STRIP_D[order(match(STRIP_D$side, names(SIDE)), match(STRIP_D$check, names(CHECK))), ]
LD$panel <- factor(STRIP_D$lab[match(paste(LD$side, LD$check), paste(STRIP_D$side, STRIP_D$check))], STRIP_D$lab)
LD$col <- SIDECOL[LD$side]; LD$below <- LD$after < 1
# named genes: those below 2-fold after a check, and the two largest drops in each plot
LABD <- do.call(rbind, lapply(split(LD, LD$panel), function(d) { d$drop <- d$before - d$after
  d[d$gene %in% unique(c(d$gene[d$below], head(d$gene[order(-d$drop)], 2))), ] }))
LIM_D <- c(0, ceiling(max(c(LD$before, LD$after)))); YLIM_D <- c(min(0, floor(min(LD$after))), LIM_D[2])
pD <- ggplot(LD, aes(before, after)) +
  annotate("rect", xmin = -Inf, xmax = Inf, ymin = -Inf, ymax = 1, fill = "grey93") +
  geom_abline(slope = 1, intercept = 0, colour = "grey50", linewidth = .3, linetype = "22") +
  geom_hline(yintercept = 1, colour = "grey45", linewidth = .3) +
  geom_point(data = LD[!LD$below, ], aes(fill = I(col)), shape = 21, colour = "black", size = .7, stroke = .08) +
  geom_point(data = LD[LD$below, ], aes(fill = I(col)), shape = 21, colour = "black", size = 1.6, stroke = .4) +
  geom_text_repel(data = LABD, aes(label = gene, fontface = ifelse(below, "bold.italic", "italic")), size = MM, family = FONT, colour = "black",
                  bg.colour = "white", bg.r = .12, min.segment.length = 0, segment.size = .2, box.padding = .3, point.padding = .1,
                  nudge_x = 1.2, nudge_y = -.6, max.overlaps = Inf, seed = 1, max.time = 3) +
  facet_wrap(~panel, nrow = 1) +
  scale_x_continuous(limits = LIM_D, breaks = c(0, 5, 10), labels = minus, expand = expansion(add = .3)) +
  scale_y_continuous(limits = YLIM_D, breaks = c(0, 5, 10), labels = minus, expand = expansion(add = .3)) +
  labs(x = "Fold change in Fig. 2C (log2, toward the higher group)", y = "Fold change after\nthe check (log2)") +
  theme_pg("both") + theme(panel.spacing.x = unit(.12, "in"))

draw_panel <- function() draw_after_72dpi(function() {
  letter("D", LX, ROW3); place(pD, X0, ROW3, XR - X0, D_H)
  key_row(list(item("Higher in intestinal EC", function(x, y) sym_dot(x, y, COL[["gut_ec"]], .06), .13), item("Higher in SC-SIC", function(x, y) sym_dot(x, y, COL[["sc_sic"]], .06), .13),
               item("No change", function(x, y) sym_line(x, y, "grey50", lty = "22", lwd = .3, w = .26), .32),
               item("Below 2-fold after the check", function(x, y) { sym_box(x, y, "grey93", .2, .12); grid.points(unit(x, "in"), unit(H - y, "in"), pch = 21, size = unit(.08, "in"), gp = gpar(col = "black", fill = "white", lwd = .4 * .stroke / 2)) }, .26)),
          X0, XR - X0, KEY_D_Y)
})
