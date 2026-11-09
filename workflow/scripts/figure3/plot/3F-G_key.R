# Figures 3F and 3G, shared key: the color of each gene group and of the mixed genes with what the group holds, in two
# rows of two, centered under both panels.

# Key items in a grid, row by row (ncol columns, each as wide as its widest item), centered in [left, left + width];
# first row top at y.
key_grid <- function(items, left, width, y, ncol, gap = .25, step = .16) {
  widths <- vapply(items, function(it) it$w_sym + str_w(it$label), 0)
  col <- (seq_along(items) - 1) %% ncol + 1; row <- (seq_along(items) - 1) %/% ncol + 1
  cw <- vapply(seq_len(ncol), function(k) max(widths[col == k]), 0)
  total <- sum(cw) + gap * (ncol - 1); stopifnot(total <= width)
  x0 <- left + (width - total) / 2 + c(0, cumsum(cw + gap))[col]
  for (i in seq_along(items)) { yy <- y + (row[i] - 1) * step; items[[i]]$draw(x0[i], yy); text_at(items[[i]]$label, x0[i] + items[[i]]$w_sym, yy - .005) }
}

draw_panel <- function()
  key_grid(list(dot_item("Group 1 (G1): EC genes low in SC-SIC", GROUP_COLS[["G1"]]), dot_item("Group 2 (G2): EC genes SC-SIC have", GROUP_COLS[["G2"]]),
                dot_item("Group 3 (G3): SC-SIC's own program", GROUP_COLS[["G3"]]), dot_item("Mixed: 2 of 3 intestinal studies", GROUP_COLS[["Mixed"]])),
           .10, 7.29, 6.595, ncol = 2)
