# Figure S8B: the top pathway terms (GO biological process, Reactome, Hallmark; left) and MSigDB C8 cell-type
# signatures (right) called for each of the six HIO gene groups of A (three patterns, each split by direction), as bars
# of -log10 adjusted p, the term's name on its bar. Left of each bar: the term's genes in the group over its genes in
# the background. A group without a called term says so.

suppressPackageStartupMessages(library(data.table))
genes <- as.data.table(read_table("s4_hio/genes_all.tsv.gz"))[expressed == TRUE]
PW <- as.data.table(read_table("s4_hio/enrichment/enrichment_called.tsv"))
CT <- as.data.table(read_table("s4_hio/enrichment_celltype/celltype_called.tsv"))
for (d in list(PW, CT)) {
  d[, K := as.integer(sub("/.*", "", bg_ratio))]                  # the term's genes in the background
  stopifnot(!is.na(d$K), d$K >= d$genes_in_list_and_term, d$K >= 10L, d$K <= 500L, d$called, d$p_adj < .05)
}
GRP <- c(stem_up = "HIO graft EC and\nSC-SIC higher", stem_down = "HIO graft EC and\nSC-SIC lower",
         gut_up = "SC-SIC higher;\nHIO graft EC closer to\nnative intestinal EC", gut_down = "SC-SIC lower;\nHIO graft EC closer to\nnative intestinal EC",
         hio_up = "HIO graft EC\nhigher; SC-SIC closer\nto native intestinal EC", hio_down = "HIO graft EC\nlower; SC-SIC closer\nto native intestinal EC")
# group sizes: the expressed genes of each pattern, split by direction
n_grp <- c(stem_up = "higher in HIO and SC-SIC", stem_down = "lower in HIO and SC-SIC", gut_up = "higher in SC-SIC", gut_down = "lower in SC-SIC",
           hio_up = "higher in HIO", hio_down = "lower in HIO")
GLAB <- sprintf("%s (%d)", GRP, vapply(n_grp, function(d) sum(genes$direction == d), 0L)); names(GLAB) <- names(GRP)
GCOL <- c(stem_up = SIC, stem_down = SIC, gut_up = "#43A5B6", gut_down = "#43A5B6", hio_up = "#6943AC", hio_down = "#6943AC")
NONE <- "none at adjusted p < 0.05"

# ---- the drawn terms: called terms by adjusted p (then name), a term whose gene set equals a better one's dropped; the
# top 3 pathways, and the top 3 cell types with at most one per source study and tissue
dedup <- function(d) { d <- d[order(p_adj, term)]; d[!duplicated(vapply(strsplit(genes, "/", fixed = TRUE), function(z) paste(sort(z), collapse = "/"), ""))] }
SRC <- list("^MANNO_MIDBRAIN_" = "La Manno | midbrain", "^HE_LIM_SUN_FETAL_LUNG_" = "He, Lim, Sun | fetal lung", "^FAN_OVARY_" = "Fan | ovary",
            "^FAN_EMBRYONIC_CTX_" = "Fan | embryonic cortex", "^HU_FETAL_RETINA_" = "Hu | fetal retina", "^TRAVAGLINI_LUNG_" = "Travaglini | lung",
            "^ZHONG_PFC_" = "Zhong | fetal prefrontal cortex", "^GAO_LARGE_INTESTINE_" = "Gao | large intestine", "^LAKE_ADULT_KIDNEY_" = "Lake | adult kidney")
source_of <- function(t) { for (k in names(SRC)) if (grepl(k, t)) return(SRC[[k]])
  m <- regmatches(t, regexec("^(DESCARTES_FETAL|BUSSLINGER)_([A-Z]+)_", t))[[1]]; if (length(m) == 3) return(paste(m[2], m[3], sep = " | "))
  stop("cell-type term without a source study and tissue: ", t) }
# display names: sentence case, acronyms kept, Hallmark names in their usual form
PW_NAME <- c("Coagulation pathway" = "Coagulation pathway", "homophilic cell-cell adhesion" = "Homophilic cell–cell adhesion",
             HALLMARK_TNFA_SIGNALING_VIA_NFKB = "TNF-α signaling via NF-κB", "T cell differentiation" = "T cell differentiation",
             "lymphocyte differentiation" = "Lymphocyte differentiation", "anterograde trans-synaptic signaling" = "Anterograde trans-synaptic signaling",
             "trans-synaptic signaling" = "Trans-synaptic signaling", "modulation of chemical synaptic transmission" = "Modulation of chemical synaptic transmission",
             "regulation of cell adhesion" = "Regulation of cell adhesion", "response to bacterium" = "Response to bacterium",
             "cytokine production" = "Cytokine production", HALLMARK_INTERFERON_ALPHA_RESPONSE = "Interferon-α response",
             HALLMARK_INTERFERON_GAMMA_RESPONSE = "Interferon-γ response", HALLMARK_EPITHELIAL_MESENCHYMAL_TRANSITION = "Epithelial–mesenchymal transition",
             "Attachment of bacteria to epithelial cells" = "Attachment of bacteria to epithelial cells", "ECM proteoglycans" = "ECM proteoglycans")
CT_NAME <- c(FAN_OVARY_CL7_ANGEIOGENIC_ENDOTHELIAL_CELL = "Ovary endothelial", HU_FETAL_RETINA_FIBROBLAST = "Fetal retina fibroblast",
             TRAVAGLINI_LUNG_EREG_DENDRITIC_CELL = "Lung dendritic", MANNO_MIDBRAIN_NEUROTYPES_HDA1 = "Midbrain dopamine neurons",
             DESCARTES_FETAL_INTESTINE_ENS_NEURONS = "Fetal enteric neurons", ZHONG_PFC_C3_MICROGLIA = "Fetal cortex microglia",
             DESCARTES_FETAL_INTESTINE_INTESTINAL_EPITHELIAL_CELLS = "Fetal intestinal epithelium", BUSSLINGER_DUODENAL_MATURE_ENTEROCYTES = "Duodenal enterocytes",
             DESCARTES_FETAL_KIDNEY_URETERIC_BUD_CELLS = "Fetal kidney ureteric bud", HE_LIM_SUN_FETAL_LUNG_C1_PULMONARY_NEUROENDOCRINE_CELL = "Fetal lung neuroendocrine",
             DESCARTES_FETAL_CEREBELLUM_OLIGODENDROCYTES = "Fetal cerebellum oligodendrocytes", GAO_LARGE_INTESTINE_ADULT_CA_ENTEROENDOCRINE_CELLS = "Large intestine enteroendocrine",
             MANNO_MIDBRAIN_NEUROTYPES_HNBML5 = "Embryonic midbrain neuroblasts", GAO_LARGE_INTESTINE_ADULT_CJ_IMMUNE_CELLS = "Large intestine immune",
             LAKE_ADULT_KIDNEY_C26_MESANGIAL_CELLS = "Kidney mesangial")
none_row <- function(l) data.table(g = l, name = NONE, v = NA_real_, n = NA_integer_, K = NA_integer_, term = NA_character_)
terms_of <- function(d, l, NM) {
  nm <- unname(NM[d$term]); stopifnot(!anyNA(nm))                # every drawn term has a display name
  data.table(g = l, name = nm, v = -log10(d$p_adj), n = d$genes_in_list_and_term, K = d$K, term = d$term)
}
pick_pw <- function(l) { d <- dedup(PW[list == l]); if (!nrow(d)) none_row(l) else terms_of(head(d, 3), l, PW_NAME) }
pick_ct <- function(l) { d <- dedup(CT[list == l]); if (!nrow(d)) return(none_row(l))
  out <- NULL; seen <- character()
  for (i in seq_len(nrow(d))) { if (!is.null(out) && nrow(out) == 3) break
    s <- source_of(d$term[i]); if (s %in% seen) next; seen <- c(seen, s); out <- rbind(out, d[i]) }
  terms_of(out, l, CT_NAME) }
A <- rbindlist(lapply(names(GRP), pick_pw)); B <- rbindlist(lapply(names(GRP), pick_ct))
stopifnot(nrow(CT[list == "stem_up"]) == 0L)                      # higher in both: no cell-type signature

# ---- rows: per group the larger of the two columns' term counts; row k of each column side by side
A[, r := seq_len(.N), by = g]; B[, r := seq_len(.N), by = g]
A[, nk := ifelse(is.na(n), NA_character_, sprintf("%d / %d", n, K))]; B[, nk := ifelse(is.na(n), NA_character_, sprintf("%d / %d", n, K))]
NR <- sapply(names(GRP), function(l) max(sum(A$g == l), sum(B$g == l)))
GRID <- rbindlist(lapply(names(GRP), function(l) data.table(g = l, r = seq_len(NR[[l]])))); GRID[, id := paste(g, r)]
rows_of <- function(d) { d <- merge(GRID, d, by = c("g", "r"), all.x = TRUE, sort = FALSE); d <- d[match(GRID$id, paste(d$g, d$r))]
  d[is.na(name), name := ""]; d[, key := factor(paste(g, r), rev(GRID$id))]; d[, gl := factor(GLAB[g], GLAB)]; d }
AD <- rows_of(A); BD <- rows_of(B)
LP <- .02; G2 <- .04                                     # inches: left pad; count to bar origin
TB <- c(P = LP + max(str_w_at(c(na.omit(AD$nk), "group / term genes"), 600, 7)) + G2,
        C = LP + max(str_w_at(c(na.omit(BD$nk), "group / term genes"), 600, 7)) + G2)   # count block widths (names sit on the bars)
HI <- c(P = max(AD$v, na.rm = TRUE) * 1.04, C = max(BD$v, na.rm = TRUE) * 1.04)
BRK <- list(P = c(0, 5, 10), C = c(0, 20, 40))
col_plot <- function(d, colname, key, Wcol, strips) {
  BW <- Wcol - TB[[key]]; u <- HI[[key]] / BW; lo <- -TB[[key]] * u     # data units per inch; text block left of 0
  lead <- data.table(gl = factor(GLAB[[1]], GLAB), y = NR[[1]] + .95)
  axis_seg <- data.table(gl = factor(GLAB[[length(GLAB)]], GLAB))
  p <- ggplot(d, aes(v, key)) + geom_blank(aes(x = 0)) +
    geom_col(aes(fill = g), width = .82, alpha = .32, na.rm = TRUE) +   # pale, so the names read over them (as Fig. 3G)
    geom_text(data = d[nzchar(name) & name != NONE], aes(label = name), x = .03 * u, hjust = 0, size = 7 / .pt, family = FONT, colour = "black") +
    geom_text(data = d[name == NONE], label = 'none~at~adjusted~italic(p)<0.05', parse = TRUE, x = .03 * u, hjust = 0, size = 7 / .pt, family = FONT, colour = "black") +
    geom_text(data = d[!is.na(n)], aes(label = nk), x = -G2 * u, hjust = 1, size = 7 / .pt, family = FONT, colour = "black") +
    geom_text(data = lead, aes(x = -G2 * u, y = y), label = "group / term genes", inherit.aes = FALSE, hjust = 1, vjust = 0, size = 7 / .pt, family = FONT, colour = "black") +
    geom_segment(data = axis_seg, aes(x = 0, xend = HI[[key]], y = -Inf, yend = -Inf), inherit.aes = FALSE, colour = "black", linewidth = .35) +
    coord_cartesian(clip = "off") +
    facet_grid(gl ~ ., scales = "free_y", space = "free_y", switch = "y") +
    scale_fill_manual(values = GCOL, guide = "none") +
    scale_x_continuous(limits = c(lo, HI[[key]]), breaks = BRK[[key]], labels = minus, expand = c(0, 0)) +
    labs(x = expression("−"*log[10]*" adjusted "*italic(p)), y = NULL, title = colname) + theme_rk("x", legend = "bottom") +
    theme(axis.text.y = element_blank(), axis.ticks.y = element_blank(), panel.grid.major.y = element_blank(), axis.line.x = element_blank(),
          plot.title = element_text(face = "bold", size = PT, hjust = 0, margin = margin(b = 11)), plot.title.position = "panel",
          strip.placement = "outside", strip.text.y.left = element_text(angle = 0, hjust = 1, vjust = .5, face = "bold", size = 7, lineheight = .9,
                                                                        margin = margin(0, 6, 0, 0)),
          panel.spacing.y = unit(3.5, "pt"), legend.key.width = unit(.16, "in"), plot.margin = margin(2, 4, 2, 2),
          axis.title.x = element_text(hjust = 1, size = PT, colour = "black"))
  if (!strips) p <- p + theme(strip.text.y.left = element_blank())
  p }
build <- function(Wp, Wc) col_plot(AD, "Pathways", "P", Wp, TRUE) + col_plot(BD, "Cell-type signatures", "C", Wc, FALSE) +
  plot_layout(widths = c(Wp, Wc), guides = "collect") & theme(legend.position = "bottom", legend.justification = "left",
                                                              legend.margin = margin(0, 0, 0, 0))
# the same bar width in both columns: the layout's fixed widths, measured on a device of the plot's size
W_IN <- 7.3; H_PLOT <- 3.90 - .47
fixed_w <- function(p) { png(tempfile(fileext = ".png"), width = W_IN, height = H_PLOT, units = "in", res = 600, type = "cairo"); on.exit(dev.off())
  g <- patchworkGrob(p); pushViewport(viewport(width = unit(W_IN, "in"), height = unit(H_PLOT, "in")))
  w <- fixed_part(g$widths, convertWidth); popViewport(); w }
AVAIL <- W_IN - fixed_w(build(TB[["P"]] + 1, TB[["C"]] + 1)) - .02
BW <- (AVAIL - TB[["P"]] - TB[["C"]]) / 2
plot_B <- build(TB[["P"]] + BW, TB[["C"]] + BW)

draw_panel <- function() in_place(120, 3336, W_IN, H_PLOT, 2058, function() print(plot_B, newpage = FALSE))
