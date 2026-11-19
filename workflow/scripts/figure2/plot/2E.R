# Figure 2E: the top called terms of each gene group, pathways (top) and cell-type signatures (bottom). Bars = -log10
# adjusted p; the dot left of each bar has its area proportional to the group's genes in the term. Both columns use one x
# scale, measured on a device like the page's. Dot-size key below the panel.

# ---- the terms drawn per group --------------------------------------------------------------------------------------
# Pathways: called terms by adjusted p (then name), a term with the same gene set as a better one dropped, the top 6 for
# G3 and 3 for G2 and G1. Cell types (MSigDB C8): the same, plus at most one term per source study and tissue, top 3.
# Term size = the term's genes in the test background.
GKEY <- c(G3_SC_SIC_OWN_PROGRAM = "G3", G2_EC_GENES_SC_SIC_HAVE = "G2", G1_EC_GENES_SC_SIC_LACK = "G1")   # table values
TOP <- c(G3 = 6, G2 = 3, G1 = 3)
n_group <- table(read_table("gene_groups.tsv.gz")$group)
called <- function(d) { d$called <- d$called %in% c(TRUE, "TRUE", "True"); d }
same_set <- function(genes) duplicated(vapply(strsplit(genes, "/", fixed = TRUE), function(z) paste(sort(z), collapse = "/"), ""))
term_size <- function(bg_ratio) as.integer(sub("/.*$", "", bg_ratio))
pw <- called(read_table("pathways_called.tsv"))
pathways <- do.call(rbind, lapply(names(GKEY), function(k) {
  d <- pw[pw$list == k & pw$called, ]; d <- d[order(d$p_adj, d$term), ]; d <- head(d[!same_set(d$genes), ], TOP[[GKEY[[k]]]])
  if (!nrow(d)) return(data.frame(column = "Pathways", group = GKEY[[k]], row = 1L, name = sprintf("no term called (%d genes)", n_group[[k]]),
                                  p_adj = NA_real_, genes_in_group = NA_integer_, genes_in_term = NA_integer_))
  nm <- sub("^HALLMARK_", "", d$term); nm <- gsub("_", " ", nm); up <- grepl("^[A-Z0-9 ]+$", nm); nm[up] <- tolower(nm[up])
  data.frame(column = "Pathways", group = GKEY[[k]], row = seq_len(nrow(d)), name = nm, p_adj = d$p_adj,
             genes_in_group = d$genes_in_list_and_term, genes_in_term = term_size(d$bg_ratio)) }))
# each considered C8 term's source study and tissue, and the display name of each drawn term
SOURCE <- c(BUSSLINGER_DUODENAL_EC_CELLS = "Busslinger | adult duodenum", GAO_LARGE_INTESTINE_ADULT_CA_ENTEROENDOCRINE_CELLS = "Gao | adult large intestine",
            DESCARTES_FETAL_INTESTINE_CHROMAFFIN_CELLS = "Descartes | fetal intestine", HE_LIM_SUN_FETAL_LUNG_C1_PULMONARY_NEUROENDOCRINE_CELL = "He, Lim, Sun | fetal lung",
            HE_LIM_SUN_FETAL_LUNG_C1_INTERM_NEUROENDOCRINE_CELL = "He, Lim, Sun | fetal lung", HE_LIM_SUN_FETAL_LUNG_C1_PULMONARY_NE_PRECURSOR_CELL = "He, Lim, Sun | fetal lung",
            DESCARTES_FETAL_STOMACH_NEUROENDOCRINE_CELLS = "Descartes | fetal stomach", DESCARTES_FETAL_LUNG_VISCERAL_NEURONS = "Descartes | fetal lung",
            MANNO_MIDBRAIN_NEUROTYPES_HNBML5 = "La Manno | midbrain", MANNO_MIDBRAIN_NEUROTYPES_HDA2 = "La Manno | midbrain",
            MANNO_MIDBRAIN_NEUROTYPES_HDA = "La Manno | midbrain", DESCARTES_FETAL_LIVER_MESOTHELIAL_CELLS = "Descartes | fetal liver")
NAME <- c(BUSSLINGER_DUODENAL_EC_CELLS = "duodenal EC", GAO_LARGE_INTESTINE_ADULT_CA_ENTEROENDOCRINE_CELLS = "adult colon enteroendocrine",
          DESCARTES_FETAL_INTESTINE_CHROMAFFIN_CELLS = "fetal intestine chromaffin", HE_LIM_SUN_FETAL_LUNG_C1_PULMONARY_NEUROENDOCRINE_CELL = "fetal lung neuroendocrine",
          DESCARTES_FETAL_LUNG_VISCERAL_NEURONS = "fetal lung visceral neurons", MANNO_MIDBRAIN_NEUROTYPES_HNBML5 = "embryonic midbrain neuroblasts",
          DESCARTES_FETAL_LIVER_MESOTHELIAL_CELLS = "fetal liver mesothelium")
ct <- called(read_table("celltypes_called.tsv"))
celltypes <- do.call(rbind, lapply(names(GKEY), function(k) {
  d <- ct[ct$list == k & ct$called, ]; d <- d[order(d$p_adj, d$term_id), ]; d <- d[!same_set(d$genes), ]
  s <- SOURCE[d$term_id]; take <- which(!duplicated(s))[1:3]
  if (anyNA(take) || anyNA(s[seq_len(max(take))])) stop("fewer than 3 cell-type terms, or a considered term has no source key: ", k)
  d <- d[take, ]
  data.frame(column = "Cell-type signatures", group = GKEY[[k]], row = seq_len(nrow(d)), name = unname(NAME[d$term_id]), p_adj = d$p_adj,
             genes_in_group = d$genes_in_list_and_term, genes_in_term = term_size(d$bg_ratio)) }))
stopifnot(!anyNA(celltypes$name))
terms <- rbind(pathways, celltypes)

# ---- the plot --------------------------------------------------------------------------------------------------------
G_BOX <- at("E_plot")
G_FILL <- GROUP_COLS[c("G1", "G2", "G3")]
DOT_MAX <- 2.6; DOT_X <- 0.075                           # largest dot (ggplot size); dot center, in left of 0
LAB_OFF <- 0.04; R_PAD <- 0.04                           # in: names start LAB_OFF right of 0; padding right of the longest name
# one row per term in each column (a group with no called term: one row, "no term called (n genes)", in the group's middle
# row); row heights are equal across the two columns
g_rows <- function(column) {
  d <- terms[terms$column == column, ]
  n_rows <- setNames(pmax(1, tapply(d$row, d$group, max)[c("G1", "G2", "G3")]), c("G1", "G2", "G3"))
  grid_id <- paste(rep(names(n_rows), n_rows), unlist(lapply(n_rows, seq_len)))
  none <- is.na(d$p_adj); d$row[none] <- as.integer((n_rows[d$group[none]] + 1) %/% 2)
  m <- d[match(grid_id, paste(d$group, d$row)), ]
  data.frame(g = factor(sub(" .*", "", grid_id), c("G1", "G2", "G3")), name = ifelse(is.na(m$name), "", m$name), v = -log10(m$p_adj),
             n = m$genes_in_group, key = factor(grid_id, rev(grid_id)), col = factor(column))
}
gp_d <- g_rows("Pathways"); gc_d <- g_rows("Cell-type signatures")

on_device <- function(w, h, f) { png(tempfile(fileext = ".png"), width = w, height = h, units = "in", res = DPI, type = "cairo"); on.exit(dev.off()); f() }
g_col <- function(d, Wcol, hi, LL) {
  upi <- hi / (Wcol - LL); lo <- -LL * upi                 # x units per inch; the dot area is LL in wide
  ggplot(d, aes(v, key)) + geom_blank(aes(x = 0)) +
    geom_col(aes(fill = g), width = 0.78, alpha = 0.35, na.rm = TRUE) +
    geom_text(data = d[nzchar(d$name), ], aes(label = name), x = LAB_OFF * upi, hjust = 0, size = MM, family = FONT, colour = "black") +
    geom_point(data = d[!is.na(d$n), ], aes(size = n, fill = g), x = -DOT_X * upi, shape = 21, colour = "black", stroke = 0.15) +
    scale_size_area(max_size = DOT_MAX, limits = c(0, max(terms$genes_in_group, na.rm = TRUE)), guide = "none") +
    facet_grid(g ~ col, scales = "free_y", space = "free_y", switch = "y") +
    scale_fill_manual(values = G_FILL, guide = "none") +
    scale_x_continuous(limits = c(lo, hi), breaks = seq(0, floor(hi / 2) * 2, by = 2), expand = c(0, 0)) +
    labs(x = expression("−log10 adjusted" ~ italic(p)), y = NULL) +
    theme_fig2() + theme(axis.text.y = element_blank(), axis.ticks.y = element_blank(), panel.grid.major.y = element_blank(), panel.grid.major.x = element_blank(),
      strip.placement = "outside", strip.text.y.left = element_text(angle = 0, hjust = 0.5, vjust = 0.5, face = "bold", size = PT, colour = "black"),
      strip.text.x = element_text(face = "bold", size = PT, colour = "black", hjust = 0, margin = margin(1, 0, 2, 0)),
      panel.spacing.y = unit(2, "pt"), plot.margin = margin(1, 0, 1, 1))
}
# stacked, pathways on top, on one x scale; tick labels and the axis title only under the cell-type column
g_pair <- function(Wpanel, hi) (g_col(gp_d, Wpanel, hi, LLP) + labs(x = NULL) + theme(axis.text.x = element_blank())) /
  g_col(gc_d, Wpanel, hi, LLC) +
  plot_layout(heights = c(nrow(gp_d), nrow(gc_d))) + plot_annotation(theme = theme(plot.margin = margin(0, 0, 0, 0)))

# both columns span the box (minus its fixed parts: the group strip and margins), on one x scale that reaches the largest
# bar of either column; the longest term name must fit right of 0
w_txt <- function(s) on_device(3, 2, function() vapply(s, str_w, 0))
LLP <- 0.16; LLC <- 0.16                                 # the dot area left of 0 (in)
fixed_w <- on_device(G_BOX[["w"]], G_BOX[["h"]], function() {
  g <- patchworkGrob(g_pair(1, 4)); pushViewport(viewport(width = unit(G_BOX[["w"]], "in"), height = unit(G_BOX[["h"]], "in")))
  w <- convertWidth(g$widths, "in", valueOnly = TRUE); popViewport(); sum(w[unitType(g$widths) != "null"]) })
G_PANEL_W <- G_BOX[["w"]] - fixed_w
G_NAME_W <- max(w_txt(c(gp_d$name, gc_d$name)[nzchar(c(gp_d$name, gc_d$name))]))
if (LLP + LAB_OFF + G_NAME_W + R_PAD > G_PANEL_W) stop("2E: the longest term name does not fit in the panel: ", round(G_NAME_W, 3), " in")
plot_E <- g_pair(G_PANEL_W, max(c(gp_d$v, gc_d$v), na.rm = TRUE) * 1.04)

# key item: a dot of the size drawn for n group genes in a term
size_item <- function(n) list(label = as.character(n), w_sym = .13, draw = function(x, y)
  grid.points(x = unit(x + .045, "in"), y = unit(H - y - .04, "in"), pch = 21, size = unit(1, "char"),
              gp = gpar(col = "black", fill = "grey80", fontsize = DOT_MAX * sqrt(n / max(terms$genes_in_group, na.rm = TRUE)) * .pt + 0.15 * .stroke / 2,
                        lwd = 0.15 * .stroke / 2)))

draw_panel <- function() {
  plot_at(plot_E, "E_plot")
  k <- at("E_key"); key_row(c(list(text_item("Group genes in term:")), lapply(c(2, 5, 10, 15), size_item)), k[["x"]], k[["w"]], k[["y"]])
}
