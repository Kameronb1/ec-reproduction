# Supplementary Figure S1B: agreement with intestinal EC (Spearman rho, median over the intestinal studies) with each
# donor or sample weighted equally (colored; triangle: a study with one donor or sample), against the pooled value of
# Figure 2A (grey). Rows: intestinal studies, the six SC-islet studies (SC-SIC, SC-beta, SC-alpha in three lanes), HIO
# graft. Grey band: the range of the intestinal studies against each other.

BB <- read_table("dots_pooled_vs_balanced.tsv")
BB$study_name <- ifelse(BB$study %in% names(GSE), GSE[BB$study], BB$study)
BB$study_name[BB$study_name == "McCauley 2023 both organoids"] <- "McCauley 2023"
INT_ORDER <- c("Elmentaite 2020 + 2021", "Egozi 2021", "Hickey 2023")
LANE <- c(SC_SIC = .22, SC_BETA = 0, SC_ALPHA = -.22)                      # lanes inside an SC-islet study row
SEC_B <- c("Intestinal EC vs other EEC", "SC-islet studies", "HIO graft EC vs other EEC")
y <- 0; ysB <- numeric(0); rowsB <- data.frame()
ysB[SEC_B[1]] <- y; y <- y - 1
for (s in INT_ORDER) { rowsB <- rbind(rowsB, data.frame(name = s, y = y, sec = 1)); y <- y - 1 }
y <- y - .25; ysB[SEC_B[2]] <- y; y <- y - 1.15
for (s in SC_ORDER) { rowsB <- rbind(rowsB, data.frame(name = s, y = y, sec = 2)); y <- y - 1.3 }
y <- y + .15 - .25; ysB[SEC_B[3]] <- y; y <- y - 1
rowsB <- rbind(rowsB, data.frame(name = "McCauley 2023", y = y, sec = 3))
BB$y <- rowsB$y[match(BB$study_name, rowsB$name)] + ifelse(BB$population %in% names(LANE), LANE[BB$population], 0)
if (anyNA(BB$y)) stop("panel B: a study without a row")
BB$one <- BB$flag_min5 == "ONE_UNIT"; BB$col <- POPCOL[BB$population]
BAND_B <- range(BB$balanced_min5[BB$population == "GUT_EC"])
YL_B <- c(min(rowsB$y) - .6, max(ysB) + .6)
SEP_B <- ysB[-1] + .62
SECR_B <- do.call(rbind, lapply(1:3, function(q) { r <- rowsB$y[rowsB$sec == q]; data.frame(lo = min(r) - .5, hi = max(r) + .5) }))
GRID_B <- do.call(rbind, lapply(c(0, .2, .4), function(v) data.frame(v = v, SECR_B, col = if (v == 0) "grey55" else "#E4E4E4", lw = if (v == 0) .25 else .2)))
pB_dots <- ggplot(BB) +
  geom_rect(data = SECR_B, aes(xmin = BAND_B[1], xmax = BAND_B[2], ymin = lo, ymax = hi), fill = "#E6E6E6", inherit.aes = FALSE) +   # within the sections, not over their headers
  geom_hline(yintercept = SEP_B, colour = "grey85", linewidth = .25) +
  geom_segment(data = GRID_B, aes(x = v, xend = v, y = lo, yend = hi, colour = I(col), linewidth = I(lw)), inherit.aes = FALSE) +   # zero line and grid, within the sections
  geom_segment(aes(x = pooled, xend = balanced_min5, y = y, yend = y), colour = "grey70", linewidth = .3) +
  geom_point(aes(pooled, y), colour = "grey45", size = .9, stroke = 0, shape = 16) +
  geom_point(data = BB[!BB$one, ], aes(balanced_min5, y, fill = I(col)), shape = 21, colour = "black", size = 1.55, stroke = EDGE) +
  geom_point(data = BB[BB$one, ], aes(balanced_min5, y, fill = I(col)), shape = 24, colour = "black", size = 1.45, stroke = EDGE) +
  scale_x_continuous(limits = c(-.1, .45), breaks = c(0, .2, .4), labels = c("0", "0.2", "0.4"), expand = c(0, 0)) +
  scale_y_continuous(limits = YL_B, expand = c(0, 0)) +
  labs(x = "Agreement with intestinal
EC (Spearman ρ)", y = NULL) +
  theme_pg("none") + theme(axis.text.y = element_blank())
LAB_B <- rbind(data.frame(y = ysB, txt = SEC_B, face = "bold", x = 0, hj = 0), data.frame(y = rowsB$y, txt = rowsB$name, face = "plain", x = 1, hj = 1))
pB_lab <- ggplot(LAB_B) + geom_text(aes(x, y, label = txt, hjust = hj, fontface = face), size = MM, family = FONT, vjust = .5) + coord_cartesian(clip = "off") +
  scale_x_continuous(limits = c(0, 1), expand = c(0, 0)) + scale_y_continuous(limits = YL_B, expand = c(0, 0)) + theme_void() +
  labs(x = " ") + theme(plot.margin = margin(1, 3, 1, 0), axis.title.x = element_text(size = PT, family = FONT))
B_LAB_W <- 1.17
pB <- pB_lab + pB_dots + plot_layout(widths = unit(c(B_LAB_W, 1), c("in", "null")))

draw_panel <- function() draw_after_72dpi(function() {
  letter("B", B_X - .15 - .05, ROW1); place(pB, B_X, ROW1, B_W, B_H)
  key_row(list(item("SC-SIC", function(x, y) sym_box(x, y, COL[["sc_sic"]], .1), .14), item("SC-beta", function(x, y) sym_box(x, y, COL[["beta"]], .1), .14),
               item("SC-alpha", function(x, y) sym_box(x, y, COL[["alpha"]], .1), .14)), B_X, B_W, KEY_B_Y)
  key_row(list(item("Donors or samples weighted equally", function(x, y) sym_dot(x, y, "grey45"), .14)), B_X, B_W, KEY_B_Y + .15)
  # the triangle entry says why one-sample studies sit off the pooled dot; two lines, the second under the first's text
  TRI1 <- "One donor or sample: only the intestinal"; TRI2 <- "side is reweighted"
  key_row(list(item(TRI1, function(x, y) sym_tri(x, y, "grey45"), .14)), B_X, B_W, KEY_B_Y + .30)
  tx(TRI2, B_X + (B_W - (.14 + str_in(TRI1))) / 2 + .14, KEY_B_Y + .30 + KEY_B_LINE)
  key_row(list(item("Pooled cells (Fig. 2A)", function(x, y) grid.points(unit(x, "in"), unit(H - y, "in"), pch = 16, size = unit(.045, "in"), gp = gpar(col = "grey45")), .12)),
          B_X, B_W, KEY_B_Y + .45 + KEY_B_LINE)
  key_row(list(item("Range among intestinal studies", function(x, y) sym_box(x, y, "#E6E6E6"), .22)), B_X, B_W, KEY_B_Y + .60 + KEY_B_LINE)
})
