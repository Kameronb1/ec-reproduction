# Figure S17G: Zhu 2023's Figure 5E genes after CDX2 knockout: knockout minus wild type (equal-depth pseudobulk log2
# CPM, all transferred nuclei) of four serotonin-pathway and four beta-cell identity genes in SC-SIC and SC-beta, from
# the same knockout test as Figure 5F (cdx2ko_group_test.py). It shares its row with S17F, to its right.

G <- read_table("cdx2_knockout/group_test/figure5e_genes.tsv")
G <- G[G$variant == "all transferred" & G$population %in% c("SC_EC", "SC_beta"), ]   # SC_alpha rows (5F) not drawn here
SER <- c("TPH1", "SLC18A1", "LMX1A", "DDC"); BET <- c("INS", "IAPP", "PDX1", "NKX6-1")
G$pop <- factor(c(SC_EC = "SC-SIC", SC_beta = "SC-beta")[G$population], c("SC-SIC", "SC-beta"))
G$gene <- factor(G$gene, c(SER, BET))
G$block <- factor(ifelse(G$gene %in% SER, "Serotonin pathway", "Beta-cell identity"), c("Serotonin pathway", "Beta-cell identity"))
pG <- ggplot(G, aes(gene, ko_minus_wt_log2, fill = pop)) + geom_hline(yintercept = 0, colour = REF_COL, linewidth = .25) +
  geom_point(position = position_dodge(.5), size = 2, shape = 21, colour = "black", stroke = STROKE_THIN) + facet_grid(. ~ block, scales = "free_x", space = "free_x") +
  scale_fill_manual(values = c("SC-SIC" = SC_SIC, "SC-beta" = SC_BETA), name = NULL) +
  labs(x = NULL, y = "CDX2 knockout minus\nwild type (log2)") + theme_panel("y", legend = "bottom") +
  theme(axis.text.x = element_text(face = "italic"), axis.line.y = element_line(colour = "black", linewidth = .35),
        axis.ticks.y = element_line(colour = "black", linewidth = .3), strip.text = element_text(hjust = .5, size = PT, face = "bold"),
        legend.margin = margin(0, 0, 0, 0), legend.box.spacing = unit(.04, "in"),
        legend.key.spacing.x = unit(.14, "in"), plot.margin = margin(2, 8, 4, 2))

draw_panel <- function() in_row("FG", 2, function(w, vh) print(pG, newpage = FALSE))
