# Figure 1C: the six SC-islet studies, one UMAP each (its own embedding), cells colored by their published label, with
# the SC-SIC share of the study's labeled endocrine cells printed over each UMAP. Three tiles per row, tinted, color
# key and the meaning of the percentage below.

STUDIES_C <- c("Augsornworawat 2020", "Balboa 2022", "Du 2022", "Augsornworawat 2023", "Zhu 2023", "Siehler 2024")
TINT_C <- "#AA3377"                                        # the SC-islet source color, at 7% opacity
KEY_C <- c(sc_sic = "SC-SIC", beta = "SC-beta", alpha = "SC-alpha", delta = "SC-delta", other = "Other endocrine")
# panel C's left edge and width on the page (in), right of panel B
X_C <- 3.3533333333333339; W_C <- 3.9966666666666657
SHIFT_C <- .21                                             # the content sits this much above the row's top margin

umap_C <- read_table("display_embedding_cells.tsv.gz")
share_C <- read_table("sc_sic_share.tsv")
stopifnot(identical(share_C$study, STUDIES_C), setequal(unique(umap_C$group), names(KEY_C)))

# one study's UMAP: point size shrinks with the number of cells; small groups are drawn last, on top
umap_thumb <- function(d) {
  sz <- min(.55, max(.05, 6 / sqrt(nrow(d)))) * 2.2
  d <- d[order(match(d$group, c("other", "alpha", "delta", "beta", "sc_sic"))), ]
  ggplot(d, aes(umap_1, umap_2, colour = group)) + geom_point(size = sz, shape = 16, stroke = 0) +
    scale_colour_manual(values = COL, guide = "none") +
    theme_void() + theme(panel.border = element_rect(fill = NA, colour = "grey85", linewidth = .25),
                         panel.background = element_rect(fill = "white", colour = NA), plot.margin = margin(0, 0, 0, 0))
}

# x0, y0 = panel top-left (in); w, h = panel width and height; shift = how far the content moves up
draw_C <- function(x0, y0, w, h, PH, shift) {
  panel_heading("C", x0, y0, PH)
  y0 <- y0 - shift
  top <- y0 + .40                                          # first tile row
  NC <- 3; NR <- 2; KS <- .15; KH <- 2 * KS + .06; cw <- w / NC; rh <- (y0 + h - top - KH) / NR; TH <- .30
  for (i in seq_along(STUDIES_C)) { s <- STUDIES_C[i]; x <- x0 + ((i - 1) %% NC) * cw; y <- top + ((i - 1) %/% NC) * rh
    grid.rect(unit(x + .02, "in"), unit(PH - y + .01, "in"), unit(cw - .04, "in"), unit(rh - .02, "in"), just = c("left", "top"),
              gp = gpar(fill = TINT_C, col = NA, alpha = .07))
    # study name, then the SC-SIC share, centred in the band above the UMAP
    tx(s, x + cw / 2, y + .02, PH, 7, just = c("centre", "top"))
    tx(paste("SC-SIC", share_C$printed[i]), x + cw / 2, y + .15, PH, 7, bold = TRUE, col = COL[["sc_sic"]], just = c("centre", "top"))
    pushViewport(viewport(x = unit(x + .02, "in"), y = unit(PH - y - TH, "in"), width = unit(cw - .04, "in"), height = unit(rh - .03 - TH, "in"), just = c("left", "top")))
    grid.draw(ggplotGrob(umap_thumb(umap_C[umap_C$study == s, ]))); popViewport() }
  # color key (one row), then what the percentage means
  ky <- top + NR * rh + .05
  key_row(lapply(names(KEY_C), function(k) dot_key(KEY_C[[k]], COL[[k]])), x0, w, ky, PH, gap = .12, size = 7)
  key_row(list(text_item("Percent: SC-SIC share of the labeled endocrine cells in each study")), x0, w, ky + KS, PH, size = 7)
}

draw_panel <- function() draw_after_72dpi(function() in_row("BC", function(PH)
  draw_C(X_C, .02, W_C, PH - .02 + SHIFT_C - .03, PH, SHIFT_C)))
