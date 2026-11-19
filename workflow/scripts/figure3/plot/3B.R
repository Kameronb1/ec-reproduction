# Figure 3B: the 28 G1 genes (columns) in each study (rows): log2(1 + CPM) of counts pooled over the study's samples.
B <- as.data.table(read_table("g1_pooled_levels.tsv"))
COLS <- data.table(study = c("Elmentaite 2020 + 2021", "Egozi 2021", "Hickey 2023", "McCauley 2023", "Augsornworawat 2020",
                             "Balboa 2022", "Du 2022", "Augsornworawat 2023", "Zhu 2023", "Siehler 2024"),
                   pop = c(rep("Primary intestinal EC", 3), "HIO graft EC", rep("SC-SIC", 6)))
COLS[, label := study][study %in% c("Hickey 2023", "Augsornworawat 2023"), label := paste(study, "(nuclei)")]
HDR <- c("Primary intestinal EC" = "Primary\nintestinal EC", "HIO graft EC" = "HIO graft EC", "SC-SIC" = "SC-SIC")
COLS[, grp := as.integer(factor(pop, unique(pop)))]
LV <- as.matrix(B[, COLS$study, with = FALSE])                        # 28 genes x 10 studies (drawn transposed)
LV_RANGE <- c(0, ceiling(max(LV)))
MAGMA <- viridisLite::magma(256, begin = .08, end = .96); MAGMA_AT <- seq(0, 1, length.out = 256)^1.5
lv_col <- function(v) scales::gradient_n_pal(MAGMA, MAGMA_AT, "Lab")(scales::rescale(v, from = LV_RANGE))
B_W <- 7.3
LB <- with_measure(list(gene_w = max(str_in(B$gene, PT_GENE, "italic")), lab_w = max(str_in(COLS$label, PT_GENE)),
                        hdr_w = max(unlist(lapply(strsplit(HDR, "\n"), function(s) str_in(s, PT, "bold")))),
                        t1_w = str_in("Expression, log2(1 + CPM),"), t2_w = str_in("counts pooled per study")))
B_GX <- .04 + LB$hdr_w                                               # right edge of the group names
B_BARX <- B_GX + .06                                                 # group color bar (0.07 in wide)
HM_X <- B_BARX + .07 + .06 + LB$lab_w + .05                          # left edge of the heatmap
CW <- min(.19, (B_W - HM_X - .02) / nrow(B))                         # column width; rows the same (near-square cells)
R <- CW; GAPR <- .07
if (R < .16) stop("panel B cells would be narrower than 0.16 in")
B_HEAD <- .04; B_T <- .26                                            # G1 header line, then the grid
ROW_Y <- B_T + (seq_len(nrow(COLS)) - 1) * R + (COLS$grp - 1) * GAPR  # top of each study row
HM_BOT <- max(ROW_Y) + R
GENE_X <- HM_X + (seq_len(nrow(B)) - .5) * CW
KEY_H <- .14 + .14 + .10 + .085 + .16                                # two title lines, gap, ramp, tick labels
if (.04 + max(LB$t1_w, LB$t2_w, 1.05) > HM_X - .05) stop("panel B color key reaches the gene names")
B_H <- HM_BOT + .04 + max(LB$gene_w, KEY_H) + .04
draw_B <- function() {
  .VH <<- B_H
  grid.points(XU(.09), YU(B_HEAD + .06), pch = 16, size = unit(.09, "in"), gp = gpar(col = GROUP_COL[["G1"]]))
  grid.text(GLAB[["G1"]], XU(.17), YU(B_HEAD + .06), just = c("left", "centre"), gp = gp(PT, "bold"))
  for (j in seq_len(nrow(COLS)))
    grid.rect(XU(GENE_X), YU(ROW_Y[j] + R / 2), width = unit(CW, "in"), height = unit(R, "in"),
              gp = gpar(fill = lv_col(LV[, j]), col = "white", lwd = .2 * 72.27 / 25.4))
  for (k in unique(COLS$grp)) { i <- which(COLS$grp == k); yt <- min(ROW_Y[i]); hh <- length(i) * R
    grid.rect(XU(HM_X), YU(yt), width = unit(nrow(B) * CW, "in"), height = unit(hh, "in"), just = c("left", "top"), gp = gpar(fill = NA, col = "grey75", lwd = .45))
    grid.rect(XU(B_BARX), YU(yt), width = unit(.07, "in"), height = unit(hh, "in"), just = c("left", "top"), gp = gpar(fill = POP_COL[[COLS$pop[i[1]]]], col = NA))
    grid.text(HDR[[COLS$pop[i[1]]]], XU(B_GX), YU(yt + hh / 2), just = c("right", "centre"), gp = gp(PT, "bold", lineheight = .9)) }
  grid.text(COLS$label, XU(HM_X - .05), YU(ROW_Y + R / 2), just = c("right", "centre"), gp = gp(PT_GENE))
  grid.text(B$gene, XU(GENE_X), YU(HM_BOT + .04), rot = 90, just = c("right", "centre"), gp = gp(PT_GENE, "italic"))
  # color key under the row labels, beside the gene names
  y0 <- HM_BOT + .06
  grid.text(expression("Expression, log"[2]*"(1 + CPM),"), XU(.04), YU(y0), just = c("left", "top"), gp = gp())
  grid.text("counts pooled per study", XU(.04), YU(y0 + .14), just = c("left", "top"), gp = gp())
  bw <- 1.05; yb <- y0 + .38
  ramp_bar(lv_col(seq(LV_RANGE[1], LV_RANGE[2], length.out = 256)), .04, yb, bw)
  br <- seq(0, LV_RANGE[2], 4); bx <- .04 + (br - LV_RANGE[1]) / diff(LV_RANGE) * bw
  grid.segments(XU(bx), YU(yb + .085), XU(bx), YU(yb + .085 - .025), gp = gpar(col = "white", lwd = .8))
  grid.text(br, XU(bx), YU(yb + .085 + .035), just = c("centre", "top"), gp = gp())
}
draw_panel <- function() panel_canvas(120, 2136, B_W, B_H, draw_B)
