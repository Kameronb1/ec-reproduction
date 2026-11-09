# Figure 3B: the four-marker EC classifier's marker score of SC-islet cells, one row per study. In each row the SC-SIC,
# SC-beta and other endocrine cells (mostly SC-alpha and SC-delta) each get a density scaled to its own peak; the dashed
# line is Figure 2A's cutoff; at the right, the share of SC-SIC and of SC-beta cells scoring above it. Rows are ordered
# by SC-SIC cells, most at the bottom. Two key rows below.
suppressPackageStartupMessages(library(ggridges))

B_COLS <- c("SC-SIC" = "#AA3377", "SC-beta" = "#4477AA", "Other endocrine" = "#8C8C8C")
bc <- read_table("panel_b_scores_by_cell.tsv.gz")
cutoff <- read_table("panel_b_cutoff.tsv")$cutoff
bc$who <- factor(ifelse(bc$class == "SC_SIC", "SC-SIC", ifelse(bc$class == "SC_BETA", "SC-beta", "Other endocrine")), rev(names(B_COLS)))
# share of each study's SC-SIC and SC-beta cells above the cutoff (whole percent), studies listed most SC-SIC cells first
above <- read_table("panel_b_above_cutoff.tsv")
bc$study <- factor(bc$study, unique(above$study))                                       # first level = bottom row
share <- function(cls) { a <- above[above$class == cls, ]; a$pct[match(levels(bc$study), a$study)] }
pct <- data.frame(study = factor(levels(bc$study), levels(bc$study)), sic = share("SC_SIC"), beta = share("SC_BETA"))
stopifnot(!anyNA(pct$sic), !anyNA(pct$beta))
B_TOP <- length(levels(bc$study)) + 1.0                                                  # room for the top row's curves
X_SIC <- 2.12; X_BETA <- 2.62                                                            # the two share columns (score units)

plot_B <- ggplot(bc, aes(score, study, fill = who, colour = who, group = interaction(study, who))) +
  geom_vline(xintercept = cutoff, linetype = "22", colour = "grey20", linewidth = 0.4) +
  geom_density_ridges(aes(height = after_stat(ndensity)), stat = "density", trim = TRUE, scale = 0.92, alpha = 0.25, linewidth = 0.3, bw = 0.06) +
  geom_text(data = pct, aes(x = X_SIC, y = study, label = paste0(sic, "%")), nudge_y = 0.3, hjust = 0, size = MM, family = FONT,
            colour = B_COLS[["SC-SIC"]], inherit.aes = FALSE) +
  geom_text(data = pct, aes(x = X_BETA, y = study, label = paste0(beta, "%")), nudge_y = 0.3, hjust = 0, size = MM, family = FONT,
            colour = B_COLS[["SC-beta"]], inherit.aes = FALSE) +
  scale_fill_manual(values = B_COLS, guide = "none") + scale_colour_manual(values = B_COLS, guide = "none") +
  scale_x_continuous(limits = c(-1.15, 3.1), breaks = -1:2, expand = c(0, 0)) +
  scale_y_discrete(expand = expansion(add = c(0.1, 0))) + coord_cartesian(ylim = c(1, B_TOP), clip = "off") +
  labs(x = "Four-marker score", y = NULL) +
  theme_fig3() + theme(panel.grid.major.y = element_line(colour = "#E4E4E4", linewidth = 0.2), panel.grid.major.x = element_blank())

# key item: the dashed cutoff line
dash_item <- function(label) list(label = label, w_sym = .14, draw = function(x, y)
  grid.lines(x = unit(x + .045, "in"), y = unit(H - y - c(-.01, .09), "in"), gp = gpar(col = "grey20", lwd = 0.45 * .pt, lty = "22")))

draw_panel <- function() {
  draw_at(plot_B, 3.50, .199, 3.89, 1.091)
  key_row(list(square_item("SC-SIC", B_COLS[[1]]), square_item("SC-beta", B_COLS[[2]]),
               square_item("Other endocrine cells (mostly SC-alpha, -delta)", B_COLS[[3]])), 3.50, 3.89, 1.40)
  key_row(list(dash_item(sprintf("Cutoff (Fig. 2A), %.4f", cutoff)), text_item("Right: % of SC-SIC | SC-beta above the cutoff")),
          3.50, 3.89, 1.56)
}
