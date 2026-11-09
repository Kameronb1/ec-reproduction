# Figure S17C: promoter accessibility of each gene group relative to the typical gene of the same cell type in Singh
# 2024's induced EC, D, X and T4 cells on days 5 and 6 of induction: each group's median gene (line and dots) and the
# middle half of its genes (light bar), on a fold axis. Keys in the S17A-C key below. A, B and C are built together
# (helpers.R).

t <- top_row()

draw_panel <- function() draw_top(t, 3)
