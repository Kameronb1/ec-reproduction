# Figure S12A: correlation of each sample with three references, native intestinal EC, adult islet beta and adult brain
# serotonin neurons, one row per group: SC-islet cells, adult islet beta, native intestinal EC and brain neurons. A sample
# of a reference group is compared with the mean of the other samples of its group. Black bar: group median. Open
# circles: single nuclei; filled: whole cells.

suppressPackageStartupMessages(library(data.table))

A <- as.data.table(read_table("s9_brain/similarity_units.tsv"))
A[, base := sub(" \\(nuclei study\\)$", "", grp)]
A[, nuclei := study %in% NUC_STUDY]
REF <- c(r_gut_ec = "Native intestinal EC", r_beta = "Adult islet beta", r_serotonin = "Brain serotonin neurons")
L <- melt(A, id.vars = c("group", "base", "nuclei"), measure.vars = names(REF), variable.name = "ref", value.name = "r")
L[, `:=`(ref = factor(REF[as.character(ref)], REF), col = factor(NAMES[base], NAMES),
         shp = factor(ifelse(nuclei, "Single nuclei", "Whole cells"), names(SHP)))]
L[, row := as.integer(factor(col, rev(NAMES)))]
L[, yj := row + (frank(group, ties.method = "dense") %% 7 - 3) * .07]      # fixed vertical spread within a row
L <- L[order(nuclei, -as.integer(col))]                 # filled first, open (nuclei) points on top
med <- L[, .(r = median(r)), by = .(row, ref)]
plot_A <- ggplot(L) +
  geom_rect(data = data.frame(ref = factor("Adult islet beta", REF)), aes(xmin = -Inf, xmax = Inf, ymin = -Inf, ymax = Inf),
            fill = "grey97", inherit.aes = FALSE) +                 # middle column shaded
  geom_vline(xintercept = c(.5, 1), colour = "#E4E4E4", linewidth = .2) +   # gridlines again, over the shading
  geom_vline(xintercept = 0, colour = "grey60", linewidth = .25) +
  geom_point(aes(r, yj, colour = col, shape = shp), size = 1.2, stroke = .5) +
  geom_segment(data = med, aes(x = r, xend = r, y = row - .36, yend = row + .36), linewidth = .5) +
  facet_wrap(~ref, nrow = 1) +
  scale_colour_manual(values = GCOL, drop = FALSE, guide = "none") +
  scale_shape_manual(values = SHP, drop = FALSE, name = NULL) +
  guides(shape = guide_legend(override.aes = list(size = 2, colour = "black"))) +
  scale_y_continuous(breaks = seq_along(NAMES), labels = rev(NAMES), expand = expansion(add = .5)) +
  scale_x_continuous(limits = c(-.35, 1), breaks = c(0, .5, 1), labels = mlab) +
  labs(x = "Correlation (r) with the reference; black bar, group median", y = NULL) +
  panel_theme("x", legend = "right") +
  theme(strip.text = element_text(size = PT, face = "bold", hjust = .5), legend.margin = margin(0, 0, 0, 4))

draw_panel <- function() plot_panel(plot_A, 120, 48, PW, 2.85)
