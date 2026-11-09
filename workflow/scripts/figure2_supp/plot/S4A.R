# Figure S4A: for each gene (rows: Fig. 2D's 20 genes above the black line, then the other genes reaching any
# top-20 list) and each gene-ranking dataset left out (columns, in Fig. 2D's order), whether that leave-one-out
# ranking puts the gene in its top 20, at ranks 21-30 or below 30, with the number of top-20 lists holding the gene.
# The rows above the line sit higher by the room the line's label needs.

gs <- read_table("gene_selection_by_fold.tsv")
dt <- read_table("D_effect_heatmap.tsv", "figure2")
D <- unique(dt$gene_symbol[order(dt$panel_rank)]); stopifnot(length(D) == 20L)       # Fig. 2D's 20 genes, in its order
genes <- unique(gs$gene)
cnt <- sapply(genes, function(g) sum(gs$state[gs$gene == g] == "Top 20"))             # top-20 lists holding each gene
rest <- setdiff(genes, D)
ord <- c(D[order(-cnt[D], match(D, D))], rest[order(-cnt[rest], rest)])               # by number of lists, then rank or name
# columns in Fig. 2D's order: datasets labeled by the four-marker classifier, then Elmentaite 2021, Hickey 2023,
# Egozi 2021 (primary study labels)
FOLD_ORD <- c("SF_GSE185224_ADULT_REGIONAL_INTESTINE", "SF_ELMENTAITE2020_DEVGUT", "SF_GSE125970_ADULT_INTESTINE", "SF_YU2021_ENDODERM",
              "SF_ELMENTAITE2021_GUT_ATLAS", "SF_HICKEY2023_ADULT_EEC", "SF_EGOZI2021_DEVELOPMENTAL_EEC")
stopifnot(setequal(FOLD_ORD, gs$fold))
fl <- setNames(study(sub("^SF_YU2021_ENDODERM$", "SF_YU2021_ENDODERM:fetal", FOLD_ORD)), FOLD_ORD)
xpos <- setNames(c(1:4, 4.3 + 1:3), FOLD_ORD)                                         # 0.3 gap between the two groups
R2130 <- "Ranks 21–30"
grid_df <- expand.grid(gene = ord, fold = FOLD_ORD, stringsAsFactors = FALSE)
grid_df$state <- gs$state[match(paste(grid_df$gene, grid_df$fold), paste(gs$gene, gs$fold))]
stopifnot(!anyNA(grid_df$state))
grid_df$state[grid_df$state == "Ranks 21-30"] <- R2130                                # the key prints an en dash
grid_df$x <- unname(xpos[grid_df$fold]); grid_df$y <- match(grid_df$gene, rev(ord))
ny <- length(ord); XC <- max(xpos) + .85                                             # XC: count column
lab <- data.frame(y = match(ord, rev(ord)), n = paste0(cnt[ord], "/7"))
cols <- setNames(c(PANEL, PANEL_LIGHT, "#F1F1F1"), c("Top 20", R2130, "Below 30"))
LINE <- sum(!ord %in% D) + .5                                                         # the black line under Fig. 2D's 20 genes
grp <- data.frame(x0 = c(xpos[[1]], xpos[[5]]) - .4, x1 = c(xpos[[4]], xpos[[7]]) + .4,
                  lab = c("Four-marker classifier calls", "Primary study labels"))
NOTE_A <- "Fig. 2D's 20 genes"
# `gap` (data units) opens the rows above the line for its label, which sits at note_y
plot_a <- function(gap = 0, note_y = NULL) {
  up <- function(y) y + ifelse(y > LINE, gap, 0)
  gd <- grid_df; gd$y <- up(gd$y); lb <- lab; lb$y <- up(lb$y)
  YB <- ny + .8 + gap; YT <- ny + 1.0 + gap                                           # bracket line; bracket and header baseline
  g <- ggplot(gd, aes(x, y, fill = state)) + geom_tile(color = "white", linewidth = .4, width = .95, height = .9) +
    geom_text(data = lb, aes(x = XC, y = y, label = n), inherit.aes = FALSE, hjust = 0, size = tx, family = FONT) +
    geom_segment(data = grp, aes(x = x0, xend = x1, y = YB, yend = YB), inherit.aes = FALSE, linewidth = .35, color = "black") +
    geom_text(data = grp, aes(x = (x0 + x1) / 2, y = YT, label = lab), inherit.aes = FALSE, vjust = 0, size = 7 / .pt, family = FONT, color = "black") +
    annotate("text", x = XC, y = YT, label = "Lists", hjust = 0, vjust = 0, size = 7 / .pt, family = FONT, color = "black") +
    scale_fill_manual(values = cols, breaks = names(cols), name = NULL) +
    scale_x_continuous(breaks = unname(xpos), labels = unname(fl), limits = c(.5, XC + .95), expand = c(0, 0)) +
    scale_y_continuous(breaks = up(seq_len(ny)), labels = rev(ord), limits = c(.5, ny + 1.65 + gap), expand = c(0, 0)) +
    geom_hline(yintercept = LINE, linewidth = .35, color = "black") +
    coord_cartesian(clip = "off") +
    labs(x = "Gene-ranking dataset left out", y = NULL) + theme_kit("none", legend = "bottom") +
    theme(axis.text.x = element_text(angle = 40, hjust = 1), axis.text.y = element_text(face = "italic"),
          axis.line.x = element_blank(), axis.ticks.x = element_blank())
  if (!is.null(note_y)) g <- g + annotate("text", x = Inf, y = note_y, label = NOTE_A, hjust = 1, vjust = 0, size = 7 / .pt,
                                          family = FONT, color = "black")
  g
}
# The line's label, 7 pt, sits at the line's right end with its baseline NOTE_BASE in above it. The rows above the line
# move up GAP_PX pixels and the panel grows as much, so every tile keeps its size. Data units per inch are measured on
# the plot drawn without the gap, 6.1 in tall.
GAP_PX <- 78; NOTE_BASE <- .036
unit_in <- panel_h(plot_a(), 6.7 - .60) / (ny + 1.15)                               # inches per gene row
g10 <- plot_a(gap = GAP_PX / 600 / unit_in, note_y = LINE + NOTE_BASE / unit_in)

draw_panel <- function() place(g10, 120, 48, 2502, 3774, 4.17, 6.7 - (.30 + (.17 + .13)) + GAP_PX / 600)
