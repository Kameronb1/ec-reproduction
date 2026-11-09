# Figure S17B: promoter accessibility of each gene group relative to the typical gene of the same cell type in Hickey
# 2023's native EC and other enteroendocrine cells of the small intestine and colon (donors pooled): each group's median
# gene (line and dots) and the middle half of its genes (light bar), on a fold axis. Keys in the S17A-C key below. A, B
# and C are built together (helpers.R).

t <- top_row()

draw_panel <- function() draw_top(t, 2)
