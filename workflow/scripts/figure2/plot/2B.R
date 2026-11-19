# Figure 2B: correlation (Pearson r) of each sample with three references, primary intestinal EC, adult islet beta
# cells and adult brain serotonin neurons, one row per group. A sample of a reference group is compared with the mean of
# the other samples of its group. Black bar: group median. Filled points: whole cells; rings: single nuclei. Key below.
# positions_b_rows.tsv gives each sample's drawing order and its fixed vertical step within its row (layout only).

suppressPackageStartupMessages(library(data.table))
mlab <- function(b) vapply(b, function(z) if (is.na(z)) NA_character_ else gsub("-", "−", format(z, trim = TRUE, drop0trailing = TRUE)), "")
SHP <- c("Whole cells" = 16, "Single nuclei" = 1)
NUC_STUDY <- c("GSE199636", "SF_HICKEY2023_ADULT_EEC", "Siletti 2023")   # single-nucleus studies: rings
# the group values of the table (names) and their drawn row labels and colors, top row first
NAMES <- c("Brain serotonin neurons" = "Brain serotonin neurons, adult", "Primary gut EC" = "Primary intestinal EC", "HIO EC" = "HIO graft EC",
           "SC-SIC" = "SC-SIC", "SC-beta" = "SC-beta", "SC-alpha" = "SC-alpha", "Adult islet beta" = "Adult islet beta")
GCOL <- setNames(c("#222222", "#43A5B6", "#6943AC", "#AA3377", "#4477AA", "#117733", "#1C3A63"), NAMES)

panel_theme_b <- function() {
  gl <- element_line(color = "#E4E4E4", linewidth = .2)
  theme_minimal(base_family = FONT, base_size = PT) +
    theme(text = element_text(color = "black", size = PT),
          axis.text = element_text(color = "black", size = PT),
          axis.title = element_text(color = "black", size = PT),
          strip.text = element_text(color = "black", size = 8, face = "bold", hjust = 0),
          plot.title = element_text(color = "black", size = 8, face = "bold"),
          plot.subtitle = element_text(color = "black", size = PT),
          legend.text = element_text(color = "black", size = PT),
          legend.title = element_text(color = "black", size = PT),
          legend.position = "bottom", legend.key.size = unit(.12, "in"),
          panel.grid.minor = element_blank(), panel.grid.major.x = gl, panel.grid.major.y = element_blank(),
          axis.line.x = element_line(color = "black", linewidth = .35),
          axis.ticks.x = element_line(color = "black", linewidth = .3),
          axis.ticks.length = unit(1.3, "pt"),
          panel.spacing = unit(.16, "in"),
          plot.margin = margin(2, 6, 2, 2))
}

S <- as.data.table(read_table("brain_similarity.tsv"))
pos <- as.data.table(read_placement("positions_b_rows.tsv"))
stopifnot(setequal(S$group, pos$group), !anyDuplicated(pos$group))
S <- S[match(pos$group, group)][, step := pos$step]
S[, base := sub(" \\(nuclei study\\)$", "", grp)]
S[, nuclei := study %in% NUC_STUDY]
stopifnot(all(S$base %in% names(NAMES)))
REF3 <- c(r_gut_ec = "Primary intestinal EC", r_beta = "Adult islet beta", r_serotonin = "Brain serotonin neurons")
L <- melt(S, id.vars = c("group", "base", "nuclei", "step"), measure.vars = names(REF3), variable.name = "ref", value.name = "r")
L[, `:=`(ref = factor(REF3[as.character(ref)], REF3), col = factor(NAMES[base], NAMES),
         shp = factor(ifelse(nuclei, "Single nuclei", "Whole cells"), names(SHP)))]
L[, row := as.integer(factor(col, rev(NAMES)))]
L[, yj := row + step * .07]                             # fixed vertical spread within a row
L <- L[order(nuclei, -as.integer(col))]                 # filled first, open (nuclei) points on top
med <- L[, .(r = median(r)), by = .(row, ref)]
plot_B <- ggplot(L) +
  geom_vline(xintercept = 0, colour = "grey60", linewidth = .25) +
  # whole cells: colored fill, hairline black edge; single nuclei: a colored ring with a hairline black ring just outside it;
  # the invisible layer only feeds the shape key below the panel
  geom_point(data = L[nuclei == FALSE], aes(r, yj, fill = col), shape = 21, colour = "black", size = 1.55, stroke = .15) +
  geom_point(data = L[nuclei == TRUE], aes(r, yj), shape = 1, colour = "black", size = 1.9, stroke = .15) +
  geom_point(data = L[nuclei == TRUE], aes(r, yj, colour = col), shape = 1, size = 1.2, stroke = .5) +
  geom_point(aes(r, yj, shape = shp), alpha = 0) +
  geom_segment(data = med, aes(x = r, xend = r, y = row - .36, yend = row + .36), linewidth = .5) +
  facet_wrap(~ref, nrow = 1, labeller = as_labeller(function(x) sub("Brain serotonin neurons", "Brain serotonin\nneurons", x, fixed = TRUE))) +
  scale_colour_manual(values = GCOL, drop = FALSE, guide = "none") + scale_fill_manual(values = GCOL, drop = FALSE, guide = "none") +
  scale_shape_manual(values = c("Whole cells" = 21, "Single nuclei" = 1), drop = FALSE, name = NULL) +
  guides(shape = guide_legend(override.aes = list(alpha = 1, size = c(2.3, 2), colour = "black", fill = c("grey35", NA), stroke = c(.15, .5)))) +
  scale_y_continuous(breaks = seq_along(NAMES), labels = rev(NAMES), expand = expansion(add = .5)) +
  scale_x_continuous(limits = c(-.35, 1), breaks = c(0, .5, 1), labels = mlab) +
  labs(x = "Correlation (r) with the reference; black bar, group median", y = NULL) +
  panel_theme_b() +
  theme(strip.text = element_text(size = PT, face = "bold", hjust = .5, vjust = 0, lineheight = .9), legend.margin = margin(0, 0, 0, 0),
        legend.box.spacing = unit(3, "pt"), legend.spacing.x = unit(2, "pt"),
        panel.spacing = unit(.08, "in"))

draw_panel <- function() plot_at(plot_B, "B_plot")
