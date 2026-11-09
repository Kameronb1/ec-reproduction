# Figure 1C: each study's endocrine cells on its own UMAP, colored by primary study label (Singh 2024's EC are induced
# EC, McCauley 2023's Endocrine_1 are HIO graft EC); in Wang 2020, Burclaff 2022 and Yu 2021, whose primary study
# labels have no EC, the intestinal EC are the cells the four-marker classifier classifies as EC (Figure 2). One tile per study, in Figure 1B's
# order, four per row, tinted by sample type; color key below.

studies_C <- read_table("study_roles_drawn.tsv")[c("study", "context")]
studies_C$author <- sub(" [0-9].*$", "", studies_C$study); studies_C$year <- sub("^[^0-9]* ", "", studies_C$study)
umap_C <- read_table("display_embedding_cells.tsv.gz")[c("study", "umap_1", "umap_2", "group")]
CTX_COL <- c("Native intestine" = "#43A5B6", "Induced EEC" = "#CCBB44", "SC-islet" = "#AA3377", "Adult islet" = "#4477AA",
             "HIO graft" = "#6943AC")
KEY_C <- c(gut_ec = "Intestinal EC", induced_ec = "Induced EC", hio_ec = "HIO graft EC", sc_sic = "SC-SIC", beta = "(SC-)beta",
           alpha = "(SC-)alpha", delta = "(SC-)delta", other = "Other")

# One study's UMAP: point size shrinks with the number of cells; small groups are drawn last, on top.
umap_thumb <- function(d) {
  sz <- min(.55, max(.05, 6 / sqrt(nrow(d))))
  d <- d[order(match(d$group, c("other", "alpha", "delta", "beta", "induced_ec", "hio_ec", "gut_ec", "sc_sic"))), ]
  ggplot(d, aes(umap_1, umap_2, colour = group)) + geom_point(size = sz, shape = 16, stroke = 0) +
    scale_colour_manual(values = COL, guide = "none") +
    theme_void() + theme(panel.border = element_rect(fill = NA, colour = "grey85", linewidth = .25),
                         panel.background = element_rect(fill = "white", colour = NA), plot.margin = margin(0, 0, 0, 0))
}

# 7 pt text anchored at (x, y) in the row.
tx <- function(s, x, y, PH, just = c("left", "top"))
  grid.text(s, unit(x, "in"), unit(PH - y, "in"), just = just, gp = gpar(fontfamily = "Arial", fontsize = 7, col = "black"))

# Panel C sits in the right part of the B/C row: left edge X0, content from Y0 (above the row's top edge) to Y0 + HC.
X0 <- 4.25; Y0 <- .02 - .31; WC <- 3.10; HC <- 4.10

draw_panel <- function() in_row(ROW_BC, function(PH) {
  top <- Y0 + .19 + 2 * .13 + .05                          # first tile row
  NC <- 4; NR <- 4; KS <- .12; KH <- 2 * KS + .04; cw <- WC / NC; rh <- (Y0 + HC - top - KH) / NR; TH <- .27
  for (i in seq_len(nrow(studies_C))) {
    s <- studies_C[i, ]; x <- X0 + ((i - 1) %% NC) * cw; y <- top + ((i - 1) %/% NC) * rh
    # study name and year centered between the tile top (y - .01) and the UMAP top (y + TH), on a tinted tile
    ym <- (y - .01 + y + TH) / 2
    tx(s$author, x + cw / 2, ym - .055, PH, just = c("centre", "centre"))
    tx(s$year, x + cw / 2, ym + .055, PH, just = c("centre", "centre"))
    grid.rect(unit(x + .02, "in"), unit(PH - y + .01, "in"), unit(cw - .04, "in"), unit(rh - .02, "in"), just = c("left", "top"),
              gp = gpar(fill = CTX_COL[[s$context]], col = NA, alpha = .07))
    pushViewport(viewport(x = unit(x + .02, "in"), y = unit(PH - y - TH, "in"), width = unit(cw - .04, "in"),
                          height = unit(rh - .03 - TH, "in"), just = c("left", "top")))
    grid.draw(ggplotGrob(umap_thumb(umap_C[umap_C$study == s$study, ])))
    popViewport()
  }
  # color key, two rows of four
  ky <- top + NR * rh + .04
  for (k in seq_along(KEY_C)) {
    x <- X0 + ((k - 1) %% 4) * cw + .04; y <- ky + ((k - 1) %/% 4) * KS
    grid.points(unit(x + .03, "in"), unit(PH - y - .045, "in"), pch = 16, size = unit(.075, "in"), gp = gpar(col = COL[[names(KEY_C)[k]]]))
    tx(KEY_C[[k]], x + .10, y, PH)
  }
})
