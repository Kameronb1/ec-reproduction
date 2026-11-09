# Figure S10B: gene group G2, EC genes SC-SIC has, in two parts by the Figure 3D call (rows in each part by the Figure 3D
# difference, SC-SIC minus intestinal EC, largest on top), in S10A's columns; then the column labels and the color keys that A
# and B share.

draw_panel <- function() {
  L <- s6_layout(); Y <- L$Y
  heatmap_canvas(120, 48, L$H, function() {
    block_heading(Y$h2, GROUP_COL[["G2"]], sprintf("G2: EC genes SC-SIC have (%d genes)", nrow(L$g2h) + nrow(L$g2n)))
    grid.text(sprintf("Higher in SC-SIC, Fig. 3D FDR < 0.05 (%d)", nrow(L$g2h)), XU(.04), YU(Y$sb_h), just = c("left", "top"), gp = gp(PT_STRIP, "bold"))
    draw_block(cells(L$g2h, FALSE, L$rng), L$g2h$gene, Y$b2h, L$xs, L$p, L$R, L$gl)
    grid.text(sprintf("No significant Fig. 3D difference (%d)", nrow(L$g2n)), XU(.04), YU(Y$sb_n), just = c("left", "top"), gp = gp(PT_STRIP, "bold"))
    draw_block(cells(L$g2n, FALSE, L$rng), L$g2n$gene, Y$b2n, L$xs, L$p, L$R, L$gl)
    draw_col_labels(Y$bot, L$xs, 60)
    draw_keys(Y$key, L$left, 2.55, 5.5, FALSE, L$rng)
  })
}
