# Figure S2B: the 13 clusters of the fitting cells, one row each: cluster number and its color on the S2A map,
# assigned identity and cell count, the share of its cells from each fitting group and in each four-marker class
# (stacked bars), and the reference genes that mark it. Keys for both bar colorings below.

cc <- read_table("cell_clusters.tsv.gz")
summ <- read_table("cluster_summary.tsv")
GROUP_ORDER <- c("E-MTAB-8901_paediatric", "E-MTAB-10187_intestinal", "GSE125970", "GSE185224", "E-MTAB-8901_fetal")   # S1A order
GROUP_COL <- c("E-MTAB-8901_paediatric" = "#4F5D1C", "E-MTAB-10187_intestinal" = "#2B3F6E", "GSE125970" = "#E48AAE",
               "GSE185224" = "#7A7A7A", "E-MTAB-8901_fetal" = "#B9CC8F")
ID_SHORT <- c("EC (enterochromaffin)" = "EC", "unassigned" = "Unassigned")
summ$ident <- ifelse(summ$assigned_identity %in% names(ID_SHORT), ID_SHORT[summ$assigned_identity], summ$assigned_identity)
summ$name <- sprintf("%s (%d)", summ$ident, summ$n_cells)
summ$genes <- gsub(" \\([0-9.]+\\)", "", summ$support_genes)          # "" = no supporting gene, drawn as "none"
CLASSES <- names(CALLS)
stack <- function(col, levels) {                                      # each cluster's shares, stacked left to right
  sh <- unclass(table(factor(cc$cluster, levels = summ$cluster), factor(cc[[col]], levels = levels))) / summ$n_cells
  x1 <- t(apply(sh, 1, cumsum)); list(x0 = x1 - sh, x1 = x1)
}
SG <- stack("arm", GROUP_ORDER); SC <- stack("presentation_class", CLASSES)
HDR <- c(name = "Cluster, identity (cells)", group = "Fitting group", class = "Four-marker classifier calls", genes = "Reference genes marking the cluster")
W_B <- 7.3; RP <- .22; BAR_H <- .155; N <- nrow(summ)                  # panel width, row pitch, bar height (in)
X_NUM <- .24; X_CHIP <- .30; X_NAME <- .45                           # number right-aligned at X_NUM, name left at X_NAME
B1 <- 1.62; BW <- 1.36; GAP <- .40; B2 <- B1 + BW + GAP             # the two bar columns
key_rows <- list(list(title = "Fitting group", labels = study_short(GROUP_ORDER), cols = unname(GROUP_COL[GROUP_ORDER])),
                 list(title = "Four-marker classifier calls", labels = unname(CLASS_KEY[CLASSES]), cols = unname(CALLS[CLASSES])))
KEY_GAP <- .16; SW <- .17                                             # gap between key entries, swatch and its space (in)
G0 <- B2 + max(BW + .18, str_in(HDR[["class"]], bold = TRUE) + .16)   # genes column clears the bars and the class header
key_w <- vapply(key_rows, function(k) sum(vapply(k$labels, str_in, 0) + SW) + KEY_GAP * (length(k$labels) - 1), 0)
TITLE_W <- max(vapply(key_rows, function(k) str_in(k$title, bold = TRUE), 0)) + .14
KEY_X0 <- (W_B - TITLE_W - max(key_w)) / 2                            # keys centered as a block
gp <- function(face = "plain") gpar(fontfamily = FONT, fontsize = PT, fontface = face, col = "black")
draw_b <- function() {
  vh <- convertHeight(unit(1, "npc"), "in", valueOnly = TRUE)
  hdr_y <- vh - .05                                                   # header text top
  y_row <- function(i) vh - .33 - (i - 1) * RP                        # row centers
  for (i in seq_len(N)) if (i %% 2 == 1)                              # alternating row bands
    grid.rect(unit(.02, "in"), unit(y_row(i), "in"), unit(W_B - .06, "in"), unit(RP, "in"), just = c("left", "centre"), gp = gpar(fill = "#F0F0F0", col = NA))
  grid.text(HDR[["name"]], unit(.04, "in"), unit(hdr_y, "in"), just = c("left", "top"), gp = gp("bold"))
  grid.text(HDR[["group"]], unit(B1, "in"), unit(hdr_y, "in"), just = c("left", "top"), gp = gp("bold"))
  grid.text(HDR[["class"]], unit(B2, "in"), unit(hdr_y, "in"), just = c("left", "top"), gp = gp("bold"))
  grid.text(HDR[["genes"]], unit(G0, "in"), unit(hdr_y, "in"), just = c("left", "top"), gp = gp("bold"))
  ybot <- y_row(N) - RP / 2; ytop <- y_row(1) + RP / 2
  for (b in c(B1, B2)) for (t in c(0, .5, 1))                          # light gridlines behind the bars
    grid.lines(unit(rep(b + BW * t, 2), "in"), unit(c(ybot, ytop), "in"), gp = gpar(col = GRID_COL, lwd = .6))
  for (i in seq_len(N)) {
    y <- y_row(i)
    grid.text(as.character(summ$cluster[i]), unit(X_NUM, "in"), unit(y, "in"), just = c("right", "centre"), gp = gp())
    grid.rect(unit(X_CHIP, "in"), unit(y, "in"), unit(.095, "in"), unit(.095, "in"), just = c("left", "centre"),
              gp = gpar(fill = CLUSTER_COL[[as.character(summ$cluster[i])]], col = NA))   # the cluster's color on the S2A map
    grid.text(summ$name[i], unit(X_NAME, "in"), unit(y, "in"), just = c("left", "centre"), gp = gp())
    for (k in seq_along(GROUP_ORDER)) { w <- SG$x1[i, k] - SG$x0[i, k]; if (w > 0)
      grid.rect(unit(B1 + BW * SG$x0[i, k], "in"), unit(y, "in"), unit(BW * w, "in"), unit(BAR_H, "in"), just = c("left", "centre"),
                gp = gpar(fill = GROUP_COL[[GROUP_ORDER[k]]], col = NA)) }
    for (k in seq_along(CLASSES)) { w <- SC$x1[i, k] - SC$x0[i, k]; if (w > 0)
      grid.rect(unit(B2 + BW * SC$x0[i, k], "in"), unit(y, "in"), unit(BW * w, "in"), unit(BAR_H, "in"), just = c("left", "centre"),
                gp = gpar(fill = CALLS[[CLASSES[k]]], col = NA)) }
    if (nzchar(summ$genes[i])) grid.text(summ$genes[i], unit(G0, "in"), unit(y, "in"), just = c("left", "centre"), gp = gp("italic"))
    else grid.text("none", unit(G0, "in"), unit(y, "in"), just = c("left", "centre"), gp = gp())
  }
  for (b in c(B1, B2)) {                                               # x axes: line, ticks, 0 / 50 / 100%
    grid.lines(unit(c(b, b + BW), "in"), unit(rep(ybot, 2), "in"), gp = gpar(col = "black", lwd = 1))
    for (t in c(0, .5, 1)) {
      grid.lines(unit(rep(b + BW * t, 2), "in"), unit(c(ybot, ybot - .018), "in"), gp = gpar(col = "black", lwd = .85))
      grid.text(c("0", "50", "100%")[match(t, c(0, .5, 1))], unit(b + BW * t, "in"), unit(ybot - .04, "in"), just = c("centre", "top"), gp = gp())
    }
  }
  ky <- ybot - .36                                                     # keys: one row per coloring, titles in one column
  for (r in seq_along(key_rows)) {
    k <- key_rows[[r]]; y <- ky - (r - 1) * .19; x <- KEY_X0
    grid.text(k$title, unit(x, "in"), unit(y, "in"), just = c("left", "centre"), gp = gp("bold"))
    x <- x + TITLE_W
    for (i in seq_along(k$labels)) {
      grid.rect(unit(x + .05, "in"), unit(y, "in"), unit(.095, "in"), unit(.095, "in"), gp = gpar(fill = k$cols[i], col = NA))
      grid.text(k$labels[i], unit(x + SW, "in"), unit(y, "in"), just = c("left", "centre"), gp = gp())
      x <- x + str_in(k$labels[i]) + SW + KEY_GAP
    }
  }
}
H_B <- .47 + .33 + (N - 1) * RP + RP / 2 + .36 + .19 + .12           # panel height: band, rows, axis and keys (in)

draw_panel <- function() place(draw_b, 120, 2940, 4380, 2286, W_B, H_B - (.30 + .17))
