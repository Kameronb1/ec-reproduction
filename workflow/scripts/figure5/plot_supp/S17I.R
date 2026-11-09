# Figure S17I: CDX2 knockout (Zhu 2023) and the 78 reproducible scorable Figure 5E regulons in SC-SIC (pseudobulk): each
# regulon's knockout-minus-wild-type change minus the median of its own 500 random gene sets, by the direction of its
# factor's targets (higher or lower in SC-SIC); the two CDX2(+) regulons named by network; black ticks: row medians;
# counts: regulons below 0. Keys below. It shares its row with S17H, to its right.

e <- regulon_panels()

draw_panel <- function() in_row("HI", 2, function(w, vh) e$drawI(w))
