# Figure 3H: the fetal–adult score of every row (age_scores.py: each profile's own expression of the adult-high and
# fetal-high genes, centered on the four intestinal arms; adult-high mean minus fetal-high mean). Grey bars join each
# SC-islet study's SC-beta and SC-SIC scores. The McCauley HIO graft row sits in its own unlabeled block. Dots have thin
# black outlines. Key column at the right, centered on the plot.

H_COLS <- c("fetal intestinal EC" = "#E8601C", "adult intestinal EC" = "#43A5B6", "SC-SIC" = "#AA3377", "SC-beta" = "#4477AA",
            "induced EC" = "#CCBB44", "HIO graft EC" = "#6943AC")
h <- read_table("panel_h_drawn.tsv")
stopifnot(sum(h$what == "induced EEC") == 1)          # Singh 2024 row: its EC cells are scored, so "induced EC"
h$what[h$what == "induced EEC"] <- "induced EC"; h$row <- sub("induced EEC", "induced EC", h$row, fixed = TRUE)
h$what <- sub(" gut EC$", " intestinal EC", h$what); h$block[h$block == "Gut tissue"] <- "Intestinal tissue"   # drawn names
# '*' marks the four intestinal tissue rows whose donors were used to pick the genes; the three test sources (starred in the
# table) carry a dagger instead
H_PICK <- c("Elmentaite 2020 + 2021 fetal", "Egozi 2021 fetal", "Elmentaite 2021 adult", "Hickey 2023 adult (nuclei)")
H_TEST <- c("Yu 2021 fetal", "Wang 2020 adult", "Burclaff 2022 adult")
gut <- h$block == "Intestinal tissue"
stopifnot(setequal(h$row[gut], c(H_PICK, paste(H_TEST, "*"))), !any(grepl("\\*", h$row[!gut])))
h$row[gut] <- ifelse(h$row[gut] %in% H_PICK, paste(h$row[gut], "*"), sub(" *", " †", h$row[gut], fixed = TRUE))
h$block <- factor(h$block, c("Intestinal tissue", "Grown in vitro", " "))
h$row <- factor(h$row, unique(h$row[order(h$row_from_top, decreasing = TRUE)]))
sic <- h[h$what == "SC-SIC", ]
seg <- data.frame(block = sic$block, row = sic$row, x = h$score[h$what == "SC-beta"][match(sic$row, h$row[h$what == "SC-beta"])], xend = sic$score)

plot_H <- ggplot(h, aes(score, row)) +
  geom_vline(xintercept = 0, colour = "grey70", linewidth = 0.3) +
  geom_segment(data = seg, aes(x = x, xend = xend, y = row, yend = row), colour = "grey70", linewidth = 0.5) +   # drawn first: it orders the in-vitro rows
  geom_point(aes(fill = what), size = 2.4, shape = 21, colour = "black", stroke = .1) +
  facet_grid(block ~ ., scales = "free_y", space = "free_y") +
  scale_fill_manual(values = H_COLS, guide = "none") +
  scale_y_discrete(expand = expansion(add = c(0.65, 0.45))) +
  labs(x = "Fetal–adult score (own expression; adult-high minus fetal-high genes, log2)", y = NULL) +
  theme_fig3() + theme(panel.grid.major.y = element_blank(), panel.spacing.y = unit(3, "pt"), strip.clip = "off",
                       strip.text.y = element_text(size = PT, angle = -90, hjust = 0.5))

LINE <- 1.2 * PT / 72                                 # one line of body text (in)
# A column of colored dots with labels, then one note per mark ('*', the dagger), centered on y_mid.
key_column <- function(labels, cols, left, y_mid, notes, step = .16, gap = .06) {
  n_lines <- sum(lengths(lapply(notes, `[[`, "lines")))
  y0 <- y_mid - (step * (length(labels) - 1) + .10 + gap * length(notes) + LINE * n_lines) / 2
  for (i in seq_along(labels)) { y <- y0 + (i - 1) * step; dot_item(labels[i], cols[i], outline = TRUE)$draw(left, y); text_at(labels[i], left + .14, y - .005) }
  ys <- y0 + (length(labels) - 1) * step + .10 + gap
  for (nt in notes) {
    grid.text(nt$mark, x = unit(left + .04, "in"), y = unit(H - ys + .005, "in"), just = c("centre", "top"), gp = gpar(fontfamily = FONT, fontsize = PT, col = "black"))
    for (k in seq_along(nt$lines)) text_at(nt$lines[k], left + .14, ys - .005 + (k - 1) * LINE)
    ys <- ys + LINE * length(nt$lines) + gap
  }
}

draw_panel <- function() {
  draw_at(plot_H, .10, 7.05, 6.14, 2.185)
  key_column(c("Fetal intestinal EC", "Adult intestinal EC", "SC-SIC", "SC-beta", "Induced EC", "HIO graft EC"), unname(H_COLS), 6.34, 7.05 + 2.185 / 2 - 0.10,
             list(list(mark = "*", lines = c("Used to pick", "the genes")), list(mark = "†", lines = c("Not used to", "pick the genes"))))
}
