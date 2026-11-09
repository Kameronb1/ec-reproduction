# Figure 4D: cholesterol synthesis in each row: the change (log2) of each of the 26 cholesterol synthesis genes in SC-SIC,
# SC-beta and SC-alpha (one dot per gene) and their median (black tick). Rows as in A, with grey bands on alternate rows.

D_OFF <- c(SC_SIC = -.3, SC_BETA = 0, SC_ALPHA = .3)            # each population's offset within its row
PCH <- read_table("D_cholesterol_medians.tsv")
DG <- read_table("D_cholesterol_per_gene.tsv"); DG <- DG[order(DG$comparison, DG$population, DG$gene, method = "radix"), ]
# The gene dots are spread over +-0.1 row in a fixed order (golden-ratio steps, ranked), so they do not sit on one line.
spread <- function(n, half) if (n < 2) 0 else (rank(seq_len(n) * 0.6180339887 %% 1) - (n + 1) / 2) / ((n - 1) / 2) * half
DG$y <- rowy(DG$comparison) + D_OFF[DG$population] + ave(DG$change, DG$comparison, DG$population, FUN = function(v) spread(length(v), .1))
PCH$y <- rowy(PCH$comparison) + D_OFF[PCH$population]
# the tick is the median of the drawn dots
dmed <- aggregate(change ~ comparison + population, DG, median)
stopifnot(nrow(dmed) == nrow(PCH), max(abs(dmed$change - PCH$change[match(paste(dmed$comparison, dmed$population), paste(PCH$comparison, PCH$population))])) < 1e-12)

pD <- d_frame() + scale_color_manual(values = POP_COL) +
  geom_point(data = DG, aes(x = change, y = y, color = population), size = .5, alpha = .85, stroke = 0) +
  geom_segment(data = PCH, aes(x = change, xend = change, y = y - .14, yend = y + .14), color = "black", linewidth = .5)

draw_panel <- function() shifted_up(ROW_UP[["cd"]], {
  g <- ggplotGrob(nobg(pD)); P <- row_place(g, CD_TOP, CD_GROW)
  b <- panels_in(g, D_X, P[["y"]], RIGHT - D_X, P[["h"]])
  band_rects(b)
  draw_at(g, D_X, P[["y"]], RIGHT - D_X, P[["h"]])
  cd_xtitle("Expression change (log2)", b[[1]], b[[1]], P[["y"]] + P[["h"]] + CD_XTITLE_GAP)
  col_header("26 cholesterol-synthesis genes", b[[1]])
})
