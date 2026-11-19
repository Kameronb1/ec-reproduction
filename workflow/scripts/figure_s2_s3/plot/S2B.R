# Figure S2B: percent of cells classified as EC within each published label of the four test studies, with the
# detection (dot size) and mean log expression among detecting cells (dot color) of the four markers.
lab_B <- read_table("C_primary_labels.tsv")
hio <- read_table("hio_BC_classifier_calls_by_label.tsv")
hio_calls <- read_table("intermediate/hio_four_marker_calls.tsv")
hio_tab <- aggregate(cbind(n, rule_called) ~ author_label, hio, sum)
# the McCauley 2023 per-label sums equal the HIO graft cells and the HIO graft cells classified as EC
stopifnot(sum(hio_tab$rule_called) == sum(hio_calls$arm_id == "GSE214852_WT" & hio_calls$presentation_class == "HIGH_CONFIDENCE_EC"),
          sum(hio_tab$n) == sum(hio_calls$arm_id == "GSE214852_WT"))
lab_B <- rbind(lab_B[c("family", "author_label_original", "author_label_class", "n", "called")],
               data.frame(family = "hio", author_label_original = hio_tab$author_label,
                          author_label_class = ifelse(hio_tab$author_label == "Endocrine_1", "EXPLICIT_EC", "TERMINAL_NON_EC_EEC"),
                          n = hio_tab$n, called = hio_tab$rule_called))
lab_B$rate <- 100 * lab_B$called / lab_B$n
FAM_B <- c("gut_cell_atlas", "egozi", "hickey", "hio")
EC_FILL_B <- c(gut_cell_atlas = EC_GUT, egozi = EC_GUT, hickey = EC_GUT, hio = EC_HIO)
CLASS_FILL <- c(TERMINAL_NON_EC_EEC = OTHER_EEC, BROAD_EEC_PARENT_UNRESOLVED = UNRES, ENDOCRINE_PROGENITOR = PROG)
lab_B$fill <- ifelse(lab_B$author_label_class == "EXPLICIT_EC", EC_FILL_B[lab_B$family], CLASS_FILL[lab_B$author_label_class])
stopifnot(all(lab_B$family %in% FAM_B), !anyNA(lab_B$fill), all(lab_B$rate[lab_B$author_label_class == "ENDOCRINE_PROGENITOR"] == 0))
lab_B <- lab_B[order(match(lab_B$family, FAM_B)), ]       # rows keep the table's order within each study
mm <- read_table("C_marker_detection.tsv")
GENES_B <- c("TPH1", "SLC18A1", "FEV", "LMX1A")
FILL_LIM <- c(3, 8); RAMP <- grDevices::colorRampPalette(c("#FEE5D9", "#FB6A4A", "#99000D"))(240)
fill_of <- function(v) RAMP[pmax(1, pmin(240, 1 + round(239 * (pmin(pmax(v, FILL_LIM[1]), FILL_LIM[2]) - FILL_LIM[1]) / diff(FILL_LIM))))]
DOT_MAX <- .085                                            # diameter (in) of a dot for 100% of cells with a read

PITCH_B <- .11; TITLE_H <- .15
draw_B <- function(x0, y0, w) {
  lw <- max(vapply(lab_B$author_label_original, str_w, 0)) + .06      # label column
  col_p <- .21
  xb <- x0 + lw; xd1 <- x0 + w - DOT_MAX / 2 - .02 - 3 * col_p          # first dot column centre
  bw <- (xd1 - DOT_MAX / 2 - .14 - xb) - max(str_w("123/260")) - .05    # inches for 100%
  # gene names, italic, at 45 degrees above the dot columns, ending just above the first study title
  ytop <- y0 + .40
  for (j in seq_along(GENES_B)) txt(GENES_B[j], xd1 + (j - 1) * col_p - .03, ytop - .02, just = c("left", "bottom"), rot = 45, face = "italic")
  y <- ytop; rows <- list(); k <- 0; blocks <- list()
  for (f in FAM_B) {
    txt(TITLE_B[[f]][1], x0, y + .015, face = "bold")
    txt(TITLE_B[[f]][2], x0 + str_w(TITLE_B[[f]][1], face = "bold") + .06, y + .015)
    y <- y + TITLE_H
    d <- lab_B[lab_B$family == f, ]
    for (i in seq_len(nrow(d))) {
      k <- k + 1; yc <- y + (i - .5) * PITCH_B
      if (i %% 2 == 1) rect_at(x0, yc - PITCH_B / 2, w, PITCH_B, fill = "#F2F2F2")
      rows[[k]] <- list(d = d[i, ], yc = yc)
    }
    blocks[[f]] <- c(y, y + nrow(d) * PITCH_B); y <- y + nrow(d) * PITCH_B + .03
    seg_at(x0, y - .015, x0 + w, y - .015, lwd = .5)
  }
  ybot <- y - .015
  for (bl in blocks) for (v in c(0, 50, 100)) seg_at(xb + bw * v / 100, bl[1], xb + bw * v / 100, bl[2], col = "#DADADA", lwd = .4)   # grid only behind rows
  for (r in rows) {
    d <- r$d; yc <- r$yc
    txt(d$author_label_original, xb - .05, yc, just = c("right", "centre"))
    if (d$rate > 0) rect_at(xb, yc - .036, bw * d$rate / 100, .072, fill = d$fill)
    txt(sprintf("%s/%s", format(d$called, big.mark = ","), format(d$n, big.mark = ",")), xb + bw * d$rate / 100 + .04, yc,
        just = c("left", "centre"))
    for (j in seq_along(GENES_B)) {
      m <- mm[mm$family == d$family & mm$author_label_original == d$author_label_original & mm$gene == GENES_B[j], ]
      stopifnot(nrow(m) == 1, m$n == d$n)
      if (m$cells_detecting > 0)
        grid.points(unit(xd1 + (j - 1) * col_p, "in"), unit(H - yc, "in"), pch = 21,
                    size = unit(DOT_MAX * sqrt(m$share_detecting), "in"),
                    gp = gpar(col = "black", fill = fill_of(m$mean_log_detecting), lwd = .3))
    }
  }
  seg_at(xb, ybot, xb + bw, ybot, lwd = .6)
  for (v in c(0, 50, 100)) { seg_at(xb + bw * v / 100, ybot, xb + bw * v / 100, ybot + .03, lwd = .6)
    txt(paste0(v, "%"), xb + bw * v / 100, ybot + .045, just = c("centre", "top")) }
  txt("Classified as EC (% of cells with that label)", xb + bw / 2, ybot + .045 + 1.15 * LINE, just = c("centre", "top"))
  # keys
  yk <- ybot + .045 + 2.45 * LINE + .06
  key_row(list(k_multi("EC, published label", c(EC_GUT, EC_HIO)), k_square("Other terminal EEC", OTHER_EEC),
               k_square("Unresolved EEC", UNRES)), x0, yk, width = w)
  yk2 <- yk + 1.6 * LINE
  t1 <- "Cells with a read:"; t2 <- "Mean log expression among them:"; labs <- c("25%", "50%", "100%"); sz <- c(.25, .5, 1)
  tot <- str_w(t1) + .06 + sum(vapply(labs, str_w, 0) + .17) + .2 + str_w(t2) + .06 + str_w("3") + .04 + .45 + .04 + str_w("8+")
  x <- x0 + (w - tot) / 2
  txt(t1, x, yk2, just = c("left", "centre")); x <- x + str_w(t1) + .06
  for (i in 1:3) { grid.points(unit(x + .05, "in"), unit(H - yk2, "in"), pch = 21, size = unit(DOT_MAX * sqrt(sz[i]), "in"),
                               gp = gpar(col = "black", fill = "#333333", lwd = .3))
    txt(labs[i], x + .12, yk2, just = c("left", "centre")); x <- x + str_w(labs[i]) + .17 }
  x <- x + .2; txt(t2, x, yk2, just = c("left", "centre")); x <- x + str_w(t2) + .06
  txt("3", x, yk2, just = c("left", "centre")); x <- x + str_w("3") + .04
  grid.raster(matrix(RAMP, nrow = 1), x = unit(x, "in"), y = unit(H - yk2, "in"), width = unit(.45, "in"), height = unit(.08, "in"),
              just = c("left", "centre"), interpolate = FALSE)
  txt("8+", x + .49, yk2, just = c("left", "centre"))
  yk2 + .07
}

draw_panel <- function() {
  letter("B", 2.62, .12)
  stopifnot(draw_B(2.80, .14, R - 2.80) == S2_B_BOTTOM)
}
