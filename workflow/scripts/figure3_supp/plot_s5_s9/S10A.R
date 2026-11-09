# Figure S10A: gene group G1, EC genes low in SC-SIC. Per gene (rows, strongest Figure 3D evidence on top): the per-study
# level in each intestinal and SC-islet study (magma ramp) and five fold changes (ramps colored by where the gene is higher),
# under the column-group headers that S10B shares. The column labels and keys are drawn with S10B.

draw_panel <- function() {
  L <- s6_layout()
  heatmap_canvas(120, 48, L$H, function() {
    block_heading(L$Y$h1, GROUP_COL[["G1"]], sprintf("G1: EC genes low in SC-SIC (%d genes)", nrow(L$g1)))
    draw_strips(L$Y$strip, L$xs, L$p, "EC vs other EEC")
    draw_block(cells(L$g1, FALSE, L$rng), L$g1$gene, L$Y$b1, L$xs, L$p, L$R, L$gl)
  })
}
