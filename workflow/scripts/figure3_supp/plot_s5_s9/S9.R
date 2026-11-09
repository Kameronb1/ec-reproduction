# Figure S9: Figure 3D's core genes, fold change before (x) and after (y) SoupX ambient-RNA correction or removal of
# cells with other-lineage reads, signed so the side each gene favors is positive. One facet per side and check; the grey band holds
# genes that fall below 2-fold after correction or removal (outlined points).

suppressPackageStartupMessages(library(ggrepel))

G <- read_table("s3_artifact_checks/genes_before_after.tsv.gz")
G <- G[G$core, ]                                        # Figure 3D's core genes: FDR < 0.05 and at least 2-fold
SIDE <- c(NATIVE_EC_HIGHER = "Higher in intestinal EC", SC_SIC_HIGHER = "Higher in SC-SIC")
CHECK <- c(soupx = "After ambient-RNA correction", doublets = "After removing cells with many reads from\nother cell types’ marker genes")
COL <- c(NATIVE_EC_HIGHER = "#43A5B6", SC_SIC_HIGHER = "#AA3377")
fmt <- function(n) format(n, big.mark = ",", trim = TRUE)

# one row per gene and check: fold change before and after, signed toward the side the gene favors
L <- do.call(rbind, lapply(names(SIDE), function(sd) do.call(rbind, lapply(names(CHECK), function(ck) {
  g <- G[G$call == sd, ]; s <- if (sd == "SC_SIC_HIGHER") 1 else -1
  data.frame(side = sd, check = ck, gene = g$gene, before = s * g$logFC, after = s * g[[paste0("logFC_", ck)]]) }))))
# facet strips: side, check, and how many genes stay 2-fold or more (log2 >= 1)
strip <- do.call(rbind, lapply(split(L, list(L$side, L$check)), function(d)
  data.frame(side = d$side[1], check = d$check[1],
             lab = sprintf("%s,\n%s\n%s/%s still 2-fold or more", SIDE[[d$side[1]]], sub("^After ", "after ", CHECK[[d$check[1]]]),
                           fmt(sum(d$after >= 1)), fmt(nrow(d))))))
L$panel <- factor(strip$lab[match(paste(L$side, L$check), paste(strip$side, strip$check))],
                  strip$lab[order(match(strip$side, names(SIDE)), match(strip$check, names(CHECK)))])

# labeled genes per facet: the (up to 6) genes below 2-fold with the largest drop, plus the 3 largest drops overall
lab <- do.call(rbind, lapply(split(L, L$panel), function(d) { d$drop <- d$before - d$after
  d[d$gene %in% unique(c(head(d$gene[d$after < 1][order(-d$drop[d$after < 1])], 6), head(d$gene[order(-d$drop)], 3))), ] }))
# two label layers over every point (unlabeled points carry "" and act as obstacles): other genes, and genes below
# 2-fold, which are kept inside the band to the right
L$key <- paste(L$gene, L$side, L$check); lab$key <- paste(lab$gene, lab$side, lab$check)
RP <- L; RP$lbl <- ifelse(RP$key %in% lab$key, RP$gene, ""); RP$below <- RP$after < 1
RB <- RP; RB$lbl[!RB$below] <- ""; RP$lbl[RP$below] <- ""
RP$nx <- ifelse(nzchar(RP$lbl), .5, 0); RP$ny <- ifelse(nzchar(RP$lbl), -.6, 0)
LIM <- c(0, ceiling(max(c(L$before, L$after))))
YLIM <- c(min(0, floor(min(L$after))), LIM[2])
BAND <- "Below 2-fold after correction or removal"
draw_key_band <- function(data, params, size) {         # key: the grey band with an outlined point on it
  grobTree(rectGrob(gp = gpar(fill = "grey93", col = "grey80", lwd = .5)),
           pointsGrob(.5, .5, pch = 21, size = unit(5, "pt"), gp = gpar(col = "black", fill = "white", lwd = .9))) }
minus <- scales::label_number(style_negative = "minus")

plot_S9 <- ggplot(L, aes(before, after)) +
  geom_rect(data = data.frame(k = BAND), aes(xmin = -Inf, xmax = Inf, ymin = -Inf, ymax = 1, fill = k), inherit.aes = FALSE, key_glyph = draw_key_band) +
  geom_abline(aes(slope = 1, intercept = 0, linetype = "No change"), colour = "grey55", linewidth = .3) +
  geom_hline(yintercept = 1, colour = "grey45", linewidth = .3) +
  geom_point(data = L[L$after >= 1, ], aes(colour = side), size = .8, alpha = .55) +
  geom_point(data = L[L$after < 1, ], aes(fill = side), shape = 21, colour = "black", stroke = .5, size = 1.7, show.legend = FALSE) +
  geom_text_repel(data = RP, aes(label = lbl), size = 7 / .pt, family = FONT, fontface = "italic", colour = "black",
                  bg.color = "white", bg.r = .12, nudge_x = RP$nx, nudge_y = RP$ny, min.segment.length = 0, segment.size = .2,
                  segment.colour = "grey30", box.padding = .3, point.padding = .1, max.overlaps = Inf, seed = 1, force = 1.5, max.time = 5) +
  geom_text_repel(data = RB, aes(label = lbl), size = 7 / .pt, family = FONT, fontface = "bold.italic", colour = "black",
                  bg.color = "white", bg.r = .12, nudge_x = ifelse(nzchar(RB$lbl), 1.2, 0), min.segment.length = 0, segment.size = .2, direction = "y",
                  segment.colour = "grey30", box.padding = .15, point.padding = .15, max.overlaps = Inf, seed = 1, max.time = 5,
                  xlim = c(1.9, NA), ylim = c(-.15, .9)) +
  facet_wrap(~ panel, ncol = 2) +
  scale_colour_manual(values = COL, labels = SIDE, name = NULL) +
  scale_fill_manual(values = c(COL, setNames("grey93", BAND)), breaks = BAND, name = NULL) +
  scale_linetype_manual(values = c("No change" = "22"), name = NULL) +
  scale_x_continuous(limits = LIM, expand = expansion(add = .2), labels = minus) +
  scale_y_continuous(limits = YLIM, expand = expansion(add = .2), labels = minus) +
  labs(x = expression("Fig. 3D fold change ("*log[2]*", toward the favored group)"), y = expression("Fold change after ambient-RNA correction or cell removal ("*log[2]*")")) +
  panel_theme("both", legend = "bottom") +
  guides(colour = guide_legend(order = 1, override.aes = list(size = 2, alpha = 1), theme = theme(legend.key.width = unit(.16, "in"))),
         linetype = guide_legend(order = 2, theme = theme(legend.key.width = unit(.32, "in"))),
         fill = guide_legend(order = 3, theme = theme(legend.key.width = unit(.26, "in")))) +
  theme(strip.text = element_text(face = "bold", hjust = 0, size = 7.5, lineheight = 1), legend.justification = "left",
        legend.spacing.x = unit(.18, "in"))

draw_panel <- function() plot_panel(plot_S9, 60, 48, PW, 9.0)
