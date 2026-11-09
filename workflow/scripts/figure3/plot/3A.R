# Figure 3A: the primary-study-labeled cells of each 10x SC-islet study as shares of five populations (SC-SIC,
# SC-beta, SC-alpha, SC-delta and every other label), with the study's total at the right. Two key rows below.

# Display names of the six studies, by GEO series; Augsornworawat 2023 is single-nucleus RNA.
A_STUDY <- c(GSE151117 = "Augsornworawat 2020", GSE167880 = "Balboa 2022", GSE202497 = "Zhu 2023",
             GSE226346 = "Siehler 2024", GSE185036 = "Du 2022", GSE199636 = "Augsornworawat 2023 (nuclei)")
census <- read_table("panel_a_census.tsv")
census <- census[census$block == "main" & census$part == "cells", ]          # Veres 2019 (inDrops) is in the table, not drawn
A_COLS <- c("SC-SIC" = "#AA3377", "SC-beta" = "#4477AA", "SC-alpha" = "#117733", "SC-delta" = "#AE8C44", "Other cells" = "#9FB6C4")
a <- data.frame(study = unname(A_STUDY[census$dataset]), class = census$key, cells = census$value)
a$class[!a$class %in% names(A_COLS)] <- "Other cells"                         # every other primary study label
a$study <- factor(a$study, rev(unique(a$study)))                              # the table's first study on top
a$class <- factor(a$class, rev(unique(a$class)))
a$share <- a$cells / ave(a$cells, a$study, FUN = sum)
a_tot <- aggregate(cells ~ study, a, sum); a_tot$lab <- format(a_tot$cells, big.mark = ",", trim = TRUE)

plot_A <- ggplot(a, aes(share, study, fill = class)) +
  geom_col(orientation = "y", width = 0.62, colour = "grey55", linewidth = 0.22) +
  geom_text(data = a_tot, aes(x = 1.02, y = study, label = lab), inherit.aes = FALSE, hjust = 0, size = MM, family = FONT, colour = "black") +
  scale_fill_manual(values = A_COLS) +
  scale_x_continuous(limits = c(0, 1.32), breaks = c(0, .5, 1), labels = c("0%", "50%", "100%"), expand = c(0, 0)) +
  guides(x = guide_axis(cap = "both")) +                                       # axis line from 0% to 100%; totals sit past it
  labs(x = "Share of cells (total at right)", y = NULL) +
  theme_fig3() + theme(panel.grid.major.y = element_blank())

draw_panel <- function() {
  draw_at(plot_A, .10, .107, 3.23, 1.213)
  key_row(lapply(c("SC-SIC", "SC-beta", "SC-alpha", "SC-delta"), function(s) square_item(s, A_COLS[[s]])), .10, 3.23, 1.40)
  key_row(list(square_item("Other cells (mostly non-endocrine, unresolved, progenitor)", A_COLS[["Other cells"]])), .10, 3.23, 1.56)
}
