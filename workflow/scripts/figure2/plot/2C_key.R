# Figure 2C key, below the panel: the colors of the genes higher in intestinal EC and higher in SC-SICs, the grey of the
# genes with FDR >= 0.05 (with the counts of each, from the fit table), and the ring marking the four EC marker genes.

source(file.path(dirname(PLOT_SCRIPT), "volcano_and_groups.R"))
fmt <- function(n) format(n, big.mark = ",")
N_EC <- sum(dv$call == "NATIVE_EC_HIGHER"); N_SIC <- sum(dv$call == "SC_SIC_HIGHER"); N_NS <- sum(dv$call == "NO_CALL")
ring_stroke <- unique(layer_data(plot_C, 5)$stroke)                        # the rings' stroke, so the key ring matches
KEY_C <- list(dot_item(sprintf("Higher in intestinal EC (%s)", fmt(N_EC)), "#43A5B6"),
              dot_item(sprintf("Higher in SC-SIC (%s)", fmt(N_SIC)), "#AA3377"),
              dot_item(sprintf("FDR ≥ 0.05 (%s of %s tested)", fmt(N_NS), fmt(nrow(dv))), "#C7C7C7"),
              ring_item("Four EC markers", "grey10", 2 * .pt + ring_stroke * .stroke / 2, ring_stroke * .stroke / 2))

draw_panel <- function() { k <- at("C_key"); key_grid(KEY_C, k[["x"]], k[["w"]], k[["y"]], ncol = 1) }
