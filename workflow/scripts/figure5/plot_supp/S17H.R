# Figure S17H: CDX2 knockout (Zhu 2023) and the targets of CDX2's own regulons. For each CDX2(+) regulon (network, measurable
# targets), knockout minus wild type in SC-SIC, SC-beta and SC-alpha, pseudobulk (top) and per-cell detection (bottom), against
# 500 expression-matched random gene sets (grey bar: middle 95%; tick: their median); the number right of each row is the
# fraction of random sets as low or lower. Keys below. It shares its row with S17I.

e <- regulon_panels()

draw_panel <- function() in_row("HI", 1, function(w, vh) e$drawH(w))
