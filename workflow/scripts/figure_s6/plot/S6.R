# Figure S6: G3 and serotonin-gene scores per published cell type of human fetal pancreas, one point per specimen or
# embryo, with the mean (tick) and the 95% range of the matched random gene sets (grey bar).
P <- as.data.table(read_table("points.tsv")); S <- as.data.table(read_table("rows.tsv"))
EDGE_LWD <- 0.15 * .stroke / 2                                    # hairline black edge of the points
STUDY <- c(ol = "Olaniru 2023, 12 to 20 weeks", yu = "Yu 2021b, 8 to 19 weeks")
ROWS_B <- list(ol = c("Endocrine progenitors", "Beta", "Alpha", "Delta", "Schwann", "Mesenchymal", "Acinar", "Ductal", "Endothelial", "Immune"),
               yu = c("Trunk", "Tip", "Endocrine progenitors", "Beta", "Alpha", "Delta", "Neuron", "Mesenchymal", "Acinar", "Ductal", "Endothelial", "Immune"))
N_ENDO <- c(ol = 4L, yu = 6L)                                     # endocrine and progenitor rows come first
for (s in names(ROWS_B)) { k <- unique(S[study == s, .(row, kind)])[match(ROWS_B[[s]], row)]
  if (anyNA(k$kind) || any(k$kind[seq_len(N_ENDO[[s]])] != "endo") || any(k$kind[-seq_len(N_ENDO[[s]])] != "other")) stop("row kinds") }
COLB <- c(G3 = "G3", serotonin = "serotonin")
OTHER_COL <- "#A6A6A6"; BAND_COL <- "#E3E3E3"
pt_col <- function(group, kind) ifelse(kind == "endo", GROUP_COL[ifelse(group == "G3", "G3", "G2")], OTHER_COL)
LIM <- rbind(P[, .(group, v = delta)], S[, .(group, v = null_lo)], S[, .(group, v = null_hi)])[, .(lo = min(v), hi = max(v)), by = group]
LIM[, `:=`(lo = lo - .04 * (hi - lo), hi = hi + .04 * (hi - lo))]
BRK <- list(G3 = seq(-.5, 1, .5), serotonin = -1:3)
B_W <- 7.3; RB <- .16
LBB <- with_measure(list(lab = max(str_in(paste0(unlist(ROWS_B), " (12)"), PT_GENE)), sec = max(str_in(STUDY, PT, "bold")),
                         hdr = str_in("Cell type (specimens or embryos)", PT, "bold")))
LAB_R <- .04 + max(LBB$lab, LBB$hdr) + .02
if (LBB$sec > LAB_R + .1) stop("a section header runs into the plots")
X1L <- LAB_R + .14; GAPC <- .32; PW <- (B_W - X1L - GAPC - .06) / 2
XL <- c(G3 = X1L, serotonin = X1L + PW + GAPC)
xmap <- function(v, g) { l <- LIM[group == g]; XL[[g]] + (v - l$lo) / (l$hi - l$lo) * PW }
HDR_Y <- .04; CHIP_Y <- .17
SEC_Y <- c(ol = .31, yu = NA); ROW0 <- c(ol = SEC_Y[["ol"]] + .15, yu = NA)
SEC_Y[["yu"]] <- ROW0[["ol"]] + length(ROWS_B$ol) * RB + .07; ROW0[["yu"]] <- SEC_Y[["yu"]] + .15
AX_Y <- ROW0[["yu"]] + length(ROWS_B$yu) * RB + .03
B_KEY_Y <- AX_Y + .38
B_H <- B_KEY_Y + .14
# axis title with log2 written as log with a subscript 2
AXT <- expression("Expression relative to matched random gene sets (log"[2] * "); one point per specimen or embryo")
draw_B <- function() {
  .VH <<- B_H
  set.seed(1)
  # column headers with the group chip
  grid.text("G3: separate SC-SIC program", XU(XL[["G3"]]), YU(HDR_Y), just = c("left", "top"), gp = gp(PT, "bold"))
  w1 <- str_in("Serotonin genes (G2): ", PT, "bold")
  grid.text("Serotonin genes (G2): ", XU(XL[["serotonin"]]), YU(HDR_Y), just = c("left", "top"), gp = gp(PT, "bold"))
  grid.text("TPH1, DDC, SLC18A1", XU(XL[["serotonin"]] + w1), YU(HDR_Y), just = c("left", "top"), gp = gp(PT, "bold.italic"))
  if (w1 + str_in("TPH1, DDC, SLC18A1", PT, "bold.italic") > PW) stop("serotonin header wider than its column")
  for (g in COLB) grid.rect(XU(XL[[g]]), YU(CHIP_Y), width = unit(PW, "in"), height = unit(.07, "in"), just = c("left", "top"),
                            gp = gpar(fill = GROUP_COL[[if (g == "G3") "G3" else "G2"]], col = NA))
  grid.text("Cell type (specimens or embryos)", XU(LAB_R), YU(HDR_Y), just = c("right", "top"), gp = gp(PT, "bold"))
  for (s in names(ROWS_B)) {
    rows <- ROWS_B[[s]]; y0 <- ROW0[[s]]; yc <- y0 + (seq_along(rows) - .5) * RB; ybot <- y0 + length(rows) * RB
    grid.text(STUDY[[s]], XU(.04), YU(SEC_Y[[s]]), just = c("left", "top"), gp = gp(PT, "bold"))
    nn <- S[study == s & group == "G3"][match(rows, row)]$n
    grid.text(sprintf("%s (%d)", rows, nn), XU(LAB_R), YU(yc), just = c("right", "centre"), gp = gp(PT_GENE))
    ysep <- y0 + N_ENDO[[s]] * RB
    grid.segments(XU(.04), YU(ysep), XU(XL[["serotonin"]] + PW), YU(ysep), gp = gpar(col = "grey55", lwd = .5, lty = 3))
    for (g in COLB) {
      for (b in BRK[[g]]) grid.segments(XU(xmap(b, g)), YU(y0), XU(xmap(b, g)), YU(ybot), gp = gpar(col = GRID_COL, lwd = .5))
      sb <- S[study == s & group == g][match(rows, row)]
      grid.rect(XU(xmap(sb$null_lo, g)), YU(yc), width = unit(xmap(sb$null_hi, g) - xmap(sb$null_lo, g), "in"), height = unit(.72 * RB, "in"),
                just = c("left", "centre"), gp = gpar(fill = BAND_COL, col = NA))
      grid.segments(XU(xmap(0, g)), YU(y0), XU(xmap(0, g)), YU(ybot), gp = gpar(col = "grey45", lwd = .7))
      pp <- P[study == s & group == g]; pp[, y := yc[match(row, rows)] + runif(.N, -.2, .2) * RB]
      grid.points(XU(xmap(pp$delta, g)), YU(pp$y), pch = 21, size = unit(.052, "in"),
                  gp = gpar(fill = pt_col(g, pp$kind), col = "black", lwd = EDGE_LWD))
      grid.segments(XU(xmap(sb$mean_delta, g)), YU(yc - .32 * RB), XU(xmap(sb$mean_delta, g)), YU(yc + .32 * RB),
                    gp = gpar(col = "grey15", lwd = 1.6, lineend = "butt"))
    }
  }
  for (g in COLB) {                                              # one x axis per column, under the last block
    grid.segments(XU(XL[[g]]), YU(AX_Y), XU(XL[[g]] + PW), YU(AX_Y), gp = gpar(col = "black", lwd = .7))
    bx <- xmap(BRK[[g]], g)
    grid.segments(XU(bx), YU(AX_Y), XU(bx), YU(AX_Y + .018), gp = gpar(col = "black", lwd = .6))
    grid.text(format(BRK[[g]]), XU(bx), YU(AX_Y + .035), just = c("centre", "top"), gp = gp())
  }
  grid.text(AXT, XU((XL[["G3"]] + XL[["serotonin"]] + PW) / 2), YU(AX_Y + .17), just = c("centre", "top"), gp = gp())
  two_dots <- key_item("Endocrine and progenitor cells", .24, function(x, y) {
    grid.points(XU(x + .0425 + c(0, .1)), YU(rep(y + .0425, 2)), pch = 21, size = unit(.065, "in"),
                gp = gpar(fill = GROUP_COL[c("G3", "G2")], col = "black", lwd = EDGE_LWD)) })
  grey_dot <- key_item("Other cell types", .14, function(x, y)
    grid.points(XU(x + .0425), YU(y + .0425), pch = 21, size = unit(.065, "in"), gp = gpar(fill = OTHER_COL, col = "black", lwd = EDGE_LWD)))
  key_row(list(two_dots, grey_dot, tick_item("Mean"), sq_item("95% range of the random gene sets", BAND_COL)),
          XL[["G3"]] - .3, XL[["serotonin"]] + PW - XL[["G3"]] + .3, B_KEY_Y)
}
draw_panel <- function() panel_canvas(60, 48, B_W, B_H, draw_B)
