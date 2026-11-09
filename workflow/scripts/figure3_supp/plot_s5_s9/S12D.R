# Figure S12D: STAC per sample (samples with at least 10 cells of the group): share of the sample's cells detecting it
# and its level, across EC cells, SC-islet cells, adult islet cells (adult donors only; adult islet alpha includes
# dividing alpha cells) and brain neurons; dark tick: median of the samples. Open circles: single nuclei; filled: whole
# cells. Row labels carry the number of samples.

suppressPackageStartupMessages(library(data.table))

drawn_tab <- as.data.table(read_table("s9_stac/stac_by_sample_drawn.tsv"))
ord <- c("Primary gut EC", "HIO EC", "Induced EC", "SC-SIC", "SC-alpha", "SC-beta", "SC-delta", "Adult islet alpha", "Adult islet beta", "Adult islet delta",
         "Adult raphe serotonin neurons", "Adult dopamine neurons", "Fetal serotonin neurons")
SHOWN_D <- c("Native intestinal EC", "HIO graft EC", "Induced EC", "SC-SIC", "SC-alpha", "SC-beta", "SC-delta", "Adult islet alpha", "Adult islet beta",
             "Adult islet delta", "Brain serotonin neurons, adult", "Brain dopamine neurons, adult", "Brain serotonin neurons, fetal")
gcol <- setNames(c("#43A5B6", "#6943AC", "#CCBB44", "#AA3377", "#117733", "#4477AA", "#AE8C44", "#0B4F26", "#1C3A63", "#6B5220",
                   "#222222", "#888888", "#555555"), ord)
NS_D <- drawn_tab[, .N, by = grp][match(ord, grp), N]
lab5 <- setNames(sprintf("%s (%d)", SHOWN_D, NS_D), ord)                  # row label with the number of samples
FAC <- c(pct = 'bold("Cells with ")*bolditalic("STAC")*bold(" detected (% of the sample\'s cells)")',
         lcp = 'bolditalic("STAC")*bold(" expression, ln(1 + CPM)")')
L <- rbind(drawn_tab[, .(grp, nuclei, facet = FAC[["pct"]], v = pct)], drawn_tab[, .(grp, nuclei, facet = FAC[["lcp"]], v = lcp)])
L[, `:=`(grp = factor(grp, rev(ord)), facet = factor(facet, FAC), shp = factor(ifelse(nuclei, "Single nuclei", "Whole cells"), c("Whole cells", "Single nuclei")))]
med <- L[, .(v = median(v)), by = .(grp, facet)]
Y_SEP_D <- length(ord) + 1 - c(3, 7, 10) - .5          # below the EC rows, the SC-islet rows and the adult islet rows
BLANK <- data.table(grp = factor(ord[1], rev(ord)), facet = factor(FAC[["pct"]], FAC), v = c(0, 100))   # % axis runs to 100
MIX <- drawn_tab[, .(mix = uniqueN(nuclei) == 2L), by = grp][mix == TRUE, grp]   # rows holding both nuclei and whole cells
# vertical position carries no data: in mixed rows nuclei sit 0.2 row above the line and whole cells 0.15 below, plus jitter
set.seed(1)
Lp <- copy(L)[, mix := as.character(grp) %in% MIX]
Lp[, yy := as.integer(grp) + ifelse(mix, ifelse(nuclei, .20, -.15), 0) + ifelse(mix, runif(.N, -.06, .06), runif(.N, -.15, .15))]
Lp <- Lp[order(nuclei)]                                 # nucleus points drawn last, on top
xlab_fun <- function(b) vapply(b, function(z) if (is.na(z)) NA_character_ else if (max(b, na.rm = TRUE) > 10) paste0(format(z, trim = TRUE), "%")
                                else gsub("-", "−", format(z, trim = TRUE)), "")
gd <- ggplot(L, aes(v, grp)) + geom_blank() +          # trains the y scale on all 13 groups first, in order
  geom_hline(yintercept = Y_SEP_D, colour = "grey80", linewidth = .3) +
  geom_blank(data = BLANK) +
  geom_errorbar(data = med, aes(xmin = v, xmax = v, linetype = "Median of samples"), width = .8, linewidth = .6, colour = "grey15",
                orientation = "y", key_glyph = "vpath") +
  geom_point(data = Lp, aes(v, yy, colour = grp, shape = shp), size = 1.5, stroke = .6, fill = "white") +
  facet_wrap(~facet, ncol = 2, scales = "free_x", labeller = label_parsed) +
  scale_colour_manual(values = gcol, guide = "none") +
  scale_shape_manual(values = c("Whole cells" = 16, "Single nuclei" = 21), name = NULL, drop = FALSE) +
  scale_linetype_manual(values = c("Median of samples" = 1), name = NULL) +
  scale_y_discrete(labels = lab5) + scale_x_continuous(expand = expansion(mult = c(.02, .04)), labels = xlab_fun) +
  labs(x = NULL, y = NULL) + panel_theme("x", legend = "bottom") +
  guides(shape = guide_legend(order = 1, override.aes = list(size = 1.8)), linetype = guide_legend(order = 2, override.aes = list(colour = "grey15", linewidth = .6))) +
  theme(strip.text = element_text(size = 7.5, face = "bold", hjust = 0), legend.justification = "left", panel.spacing.x = unit(.32, "in"),
        legend.margin = margin(0, 0, 0, 0), legend.key.height = unit(.13, "in"))
# a header over the row labels, on the facet headers' line and in their style, says what the bracketed numbers are
ROW_HEADER <- expression(bold("Group (samples)"))      # plotmath like the facet headers, so the baselines match
with_row_header <- function() {                         # gd drawn as a gtable with the header in the strip row over the row labels
  gt <- ggplotGrob(gd)
  s <- gt$layout[gt$layout$name == "strip-t-1-1", ]; a <- gt$layout[gt$layout$name == "axis-l-1-1", ]
  hdr <- element_grob(calc_element("strip.text.x.top", gd$theme), label = ROW_HEADER, hjust = 1, margin_y = TRUE,
                      x = unit(1, "npc") - unit(1.3 + 1.5, "pt"))   # right end where the row labels end: tick length + label margin
  grid.draw(gtable::gtable_add_grob(gt, hdr, t = s$t, l = a$l, clip = "off", name = "row-header"))
}

draw_panel <- function() plot_panel(with_row_header, 120, 3366, PW, 3.40)
