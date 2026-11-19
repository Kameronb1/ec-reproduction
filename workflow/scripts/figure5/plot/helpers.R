# Drawing helpers of Figure 5: text size, colors, theme, keys, the three panels' plots and their places on the page.
# Sourced by workflow/scripts/common/draw_panel.R, which opens the page device (7.5 x 10 in, 600 dpi, cairo) and defines
# read_table() and result_path(). Every position is in inches from the top-left corner of the current viewport (the page,
# or panel B's slot). Each plot script calls render_figure5() with its panel's letter: the page is laid out as a whole and
# only that panel is drawn.

if (.Platform$OS.type == "windows") stopifnot(nzchar(Sys.setlocale("LC_CTYPE", "English_United States.utf8")))

W <- 7.5; H <- 10                                   # page size (in)
pt <- 7.5; tm <- 7 / ggplot2::.pt                   # page text size (pt); geom_text size for 7 pt
pct <- function(d) 100 * (2^d - 1)                  # log2 difference -> % more

# colors and names shared with Figures 2 and 4
GCOL <- c(G1 = "#43A5B6", G2 = "#A89A2E", G3 = "#B4553F")          # gene groups G1-G3
POP_COL <- c("SC-SIC" = "#AA3377", "SC-beta" = "#4477AA", "SC-alpha" = "#117733")
BAND <- "#D9D9D9"                                   # random-range bars and panel B row bands

# ggplot theme for panel A
th <- function(grid_x = TRUE, grid_y = FALSE) theme_minimal(base_family = "Arial", base_size = pt) + theme(
  text = element_text(color = "black", size = pt), axis.text = element_text(color = "black", size = 7),
  axis.title = element_text(color = "black", size = 7), plot.title = element_text(color = "black", size = 7.5, face = "bold", margin = margin(0, 0, 2, 0)),
  strip.text = element_text(color = "black", size = 7, face = "bold", margin = margin(1, 0, 2, 0)),
  panel.grid.minor = element_blank(),
  panel.grid.major.x = if (grid_x) element_line(color = "#E4E4E4", linewidth = .2) else element_blank(),
  panel.grid.major.y = if (grid_y) element_line(color = "#E4E4E4", linewidth = .2) else element_blank(),
  axis.line.x = element_line(color = "black", linewidth = .35), axis.ticks.x = element_line(color = "black", linewidth = .3),
  axis.ticks.length = unit(1.3, "pt"), plot.margin = margin(2, 4, 2, 2), legend.position = "none", plot.title.position = "panel")
NOXL <- theme(axis.line.x = element_blank(), axis.ticks.x = element_blank())
zero_h <- geom_hline(yintercept = 0, color = "#777777", linewidth = .3)

# ---- text -------------------------------------------------------------------------------------------------------------
text_at <- function(s, x, y, size = pt, bold = FALSE, col = "black", just = c("left", "top"), vp_h = H)
  grid.text(s, x = unit(x, "in"), y = unit(vp_h - y, "in"), just = just, gp = gpar(fontfamily = "Arial", fontsize = size, fontface = if (bold) "bold" else "plain", col = col))
str_in <- function(s, size = pt) { pushViewport(viewport(gp = gpar(fontfamily = "Arial", fontsize = size))); v <- convertWidth(stringWidth(s), "in", valueOnly = TRUE); popViewport(); v }
# panel letter, 11 pt bold, top-left at (x, y)
heading <- function(letter, x, y, vp_h = H) text_at(letter, x, y, 11, TRUE, vp_h = vp_h)

# ---- placing plots and drawing keys -----------------------------------------------------------------------------------
draw_at <- function(p, x, y, w, h) { pushViewport(viewport(x = x / W, y = 1 - y / H, width = w / W, height = h / H, just = c("left", "top")))
  grid.draw(if (inherits(p, "gtable")) p else ggplotGrob(p)); popViewport() }
# one row of key items (type pt, dash or ln, each with a label), left-aligned or centered in [left, left + width]
row_keys <- function(items, left, width, y, gap = .16, centre = TRUE) {
  ww <- vapply(items, function(it) str_in(it$label) + if (it$type == "ln") .26 else .14, 0); total <- sum(ww) + gap * (length(ww) - 1)
  x <- if (centre) left + (width - total) / 2 else left
  for (i in seq_along(items)) { it <- items[[i]]; yy <- H - y - .045
    if (it$type == "ln") grid.lines(unit(c(x, x + .2), "in"), unit(rep(yy, 2), "in"), gp = gpar(col = it$col, lwd = it$lwd %||% 1.3, lineend = "butt"))
    if (it$type == "pt") grid.points(unit(x + .05, "in"), unit(yy, "in"), pch = it$pch %||% 16, size = unit(it$size %||% 1.5, "mm"), gp = gpar(col = it$col, fill = it$fill %||% it$col, lwd = it$plwd %||% .8))
    if (it$type == "dash") grid.lines(unit(c(x + .005, x + .095), "in"), unit(rep(yy, 2), "in"), gp = gpar(col = it$col, lwd = it$lwd %||% 1.3, lineend = "butt"))
    text_at(it$label, x + if (it$type == "ln") .26 else .14, y - .005)
    x <- x + ww[i] + gap }
  invisible(total)
}
# 5A's color ramp: magma, the same ends and the same stretch toward the low end as the RNA-level heatmaps, 100 steps
A_PAL <- scales::gradient_n_pal(viridisLite::magma(256, begin = .08, end = .96), seq(0, 1, length.out = 256)^1.5, "Lab")(seq(0, 1, length.out = 100))
# horizontal color bar in 5A's ramp with ticks and a title above it, drawn as an interpolated raster
colourbar <- function(x, y, w, lims, title) {
  grid.raster(matrix(A_PAL, nrow = 1), x = unit(x, "in"), y = unit(H - y, "in"), width = unit(w, "in"), height = unit(.075, "in"), just = c("left", "top"), interpolate = TRUE)
  br <- pretty(lims, 4); br <- br[br >= lims[1] & br <= lims[2]]
  for (b in br) { xx <- x + (b - lims[1]) / diff(lims) * w
    grid.lines(unit(rep(xx, 2), "in"), unit(H - y - c(.075, .105), "in"), gp = gpar(col = "black", lwd = .5))
    text_at(format(b), xx, y + .115, just = c("centre", "top")) }
  text_at(title, x, y - .03, just = c("left", "bottom"))
}

# ---- A: regulator RNA by cell type (heatmap) and SC-SIC minus SC-beta per SC-islet study (strip below it) ----------------
plot_A <- function() {
  # the 34 regulators in their display groups (the heatmap's column blocks), in drawing order
  GRP <- read.delim(file.path(REPO, "workflow", "scripts", "figure5", "plot", "regulator_groups.tsv"), comment.char = "#",
                    check.names = FALSE, stringsAsFactors = FALSE, fileEncoding = "UTF-8")
  SETS <- split(GRP$tf, factor(GRP$group, unique(GRP$group))); TFS <- unlist(SETS, use.names = FALSE)
  a <- read_table("regulator_rna/f5_5a_regulator_rna.tsv")
  lev <- a[a$kind == "level", ]; cnt <- a[a$kind == "sc_sic_above_sc_beta_studies", ]
  NSC <- unique(lev$studies[lev$column == "SC-SIC"]); NAD <- unique(lev$studies[lev$column == "adult islet beta"])
  lev <- rbind(lev, read_table("regulator_rna/f5_5a_regulator_rna__hio_mccauley_rows.tsv")[, names(lev)])   # HIO graft EC row
  LIMS <- range(lev$value)                          # color limits over every row of the table, before the two Elmentaite rows merge
  # Elmentaite 2020 + 2021 is one study row: the mean of its two parts
  ELM <- c("Elmentaite 2020 gut EC", "Elmentaite 2021 gut EC")
  e20 <- lev[lev$column == ELM[1], ]; e21 <- lev[lev$column == ELM[2], ]; e21 <- e21[match(e20$tf, e21$tf), ]
  stopifnot(nrow(e20) == length(TFS), !anyNA(e21$tf), all(is.finite(c(e20$value, e21$value))))
  elm <- e20; elm$column <- "Elmentaite 2020 + 2021 gut EC"; elm$value <- (e20$value + e21$value) / 2
  lev <- rbind(lev[!lev$column %in% ELM, ], elm)
  ROWA <- data.frame(column = c("Elmentaite 2020 + 2021 gut EC", "Egozi 2021 gut EC", "Hickey 2023 gut EC", "McCauley 2023 HIO graft EC",
                                "SC-SIC", "SC-beta", "SC-alpha", "adult islet beta", "adult islet alpha"),
                     label = c("Elmentaite 2020 + 2021", "Egozi 2021", "Hickey 2023 (nuclei)", "McCauley 2023",
                               "SC-SIC", "SC-beta", "SC-alpha", "Beta", "Alpha"),
                     rg = c(rep("gut", 3), "hio", rep("sc", 3), rep("ad", 2)))
  stopifnot(all(lev$column %in% ROWA$column))
  RGL <- c(gut = "Primary\nintestinal EC", hio = "HIO graft EC", sc = sprintf("SC-islet\n(%d studies)", NSC), ad = sprintf("Adult islet\n(%d studies)", NAD))
  SETL <- c(`EC / SC-SIC association` = "EC / SC-SIC\nassociation", `Endocrine association` = "Endocrine\nassociation",
            `Islet association` = "Islet\nassociation", `Intestinal association` = "Intestinal\nassociation")
  AL <- data.frame(tf = lev$tf, row = ROWA$label[match(lev$column, ROWA$column)], rg = ROWA$rg[match(lev$column, ROWA$column)], value = lev$value)
  AL$x <- factor(AL$tf, TFS); AL$set <- factor(SETL[rep(names(SETS), lengths(SETS))[match(AL$tf, TFS)]], SETL)
  AL$rg <- factor(RGL[AL$rg], RGL); AL$y <- factor(AL$row, rev(ROWA$label))
  # tiles in the published drawing order (neighboring tiles' white borders overlap; regulator_rna_draw_order.tsv)
  ORD <- read.delim(file.path(REPO, "workflow", "scripts", "figure5", "plot", "regulator_rna_draw_order.tsv"), comment.char = "#",
                    check.names = FALSE, stringsAsFactors = FALSE, fileEncoding = "UTF-8")
  o <- match(paste(ORD$row, ORD$tf), paste(AL$row, AL$tf)); stopifnot(nrow(ORD) == nrow(AL), !anyNA(o), !anyDuplicated(o))
  AL <- AL[o, ]
  pAh <- ggplot(AL, aes(x, y)) +
    geom_tile(aes(fill = value), color = "white", linewidth = .3) +
    scale_fill_gradientn(colours = A_PAL, limits = LIMS) +
    facet_grid(rg ~ set, scales = "free", space = "free", switch = "y") +
    scale_x_discrete(expand = c(0, 0)) + scale_y_discrete(expand = c(0, 0)) + labs(x = NULL, y = NULL) + th(FALSE) + NOXL +
    theme(axis.text.x = element_blank(),
          strip.placement = "outside", strip.text.y.left = element_text(angle = 0, hjust = 1, vjust = .5, face = "bold", size = 7, lineheight = .9, margin = margin(0, 4, 0, 0)),
          strip.text.x = element_text(size = 7, face = "bold", lineheight = .9, vjust = 0, margin = margin(0, 0, 2, 0)), strip.clip = "off",
          panel.spacing.x = unit(3, "pt"), panel.spacing.y = unit(2.5, "pt"), axis.text.y = element_text(size = 7, margin = margin(0, 1, 0, 0)),
          plot.margin = margin(2, 4, 0, 2), plot.background = element_blank())
  # strip: SC-SIC minus SC-beta in each SC-islet study (dots; samples with at least 20 cells averaged within a study) and
  # their mean (bar), which is the heatmap's SC-SIC row minus its SC-beta row
  sx <- read_table("regulator_rna/tf_expression_by_sample.tsv"); sx <- sx[sx$cohort == "sc" & as.logical(sx$min20) & sx$tf %in% TFS, ]
  sm <- aggregate(mean_log2cpm ~ family + population + tf, sx, mean)
  AS <- merge(sm[sm$population == "SC_SIC", c("family", "tf", "mean_log2cpm")], sm[sm$population == "SC_BETA", c("family", "tf", "mean_log2cpm")],
              by = c("family", "tf"), suffixes = c("_sic", "_beta")); AS$d <- AS$mean_log2cpm_sic - AS$mean_log2cpm_beta
  AS$x <- factor(AS$tf, TFS); AS$set <- factor(SETL[rep(names(SETS), lengths(SETS))[match(AS$tf, TFS)]], SETL)
  AM <- aggregate(d ~ tf + x + set, AS, mean)
  # the per-study differences give back the table's study counts (SC-SIC above SC-beta) and the heatmap's rows
  chk_n <- aggregate(d ~ tf, AS, function(v) c(n = length(v), up = sum(v > 0)))
  stopifnot(chk_n$d[, "n"] == NSC, chk_n$d[match(cnt$tf, chk_n$tf), "up"] == cnt$value)
  lv <- function(col) lev$value[lev$column == col][match(AM$tf, lev$tf[lev$column == col])]
  stopifnot(max(abs(AM$d - (lv("SC-SIC") - lv("SC-beta")))) <= 1e-9)
  AS$lab <- AM$lab <- " "                           # empty row-group strip: keeps the axis beside the strip; its label is drawn on the page
  # dot color: which cell type has more RNA in that study; an exact 0 is a study where neither has any
  AS$side <- ifelse(AS$d > 0, "sic", ifelse(AS$d < 0, "beta", "none"))
  stopifnot(!any(AS$side == "none" & AS$mean_log2cpm_sic != 0))
  pAs <- ggplot(AS, aes(x, d)) +
    geom_rect(data = data.frame(set = AS$set[match(TFS[c(TRUE, FALSE)], AS$tf)], k = ave(seq_along(TFS), rep(seq_along(SETS), lengths(SETS)), FUN = seq_along)[c(TRUE, FALSE)]),
              aes(xmin = k - .5, xmax = k + .5, ymin = -Inf, ymax = Inf), fill = "grey96", inherit.aes = FALSE) +   # alternating gene bands
    zero_h +
    geom_crossbar(data = AM, aes(ymin = d, ymax = d), width = .62, linewidth = .2, colour = "black") +       # mean bar, behind the dots
    ggbeeswarm::geom_quasirandom(aes(colour = side), width = .22, size = .7, stroke = 0) +
    scale_colour_manual(values = c(sic = POP_COL[["SC-SIC"]], beta = POP_COL[["SC-beta"]], none = "grey60")) +
    facet_grid(lab ~ set, scales = "free_x", space = "free_x", switch = "y") +
    scale_x_discrete(expand = expansion(add = .5)) + scale_y_continuous(breaks = c(-3, 0, 3)) +
    labs(x = NULL, y = NULL) + th(FALSE, TRUE) + NOXL +
    theme(axis.text.x = element_text(angle = 90, hjust = 1, vjust = .5, face = "italic", size = 7), strip.text.x = element_blank(),
          strip.placement = "outside", strip.text.y.left = element_text(angle = 0, hjust = 1, vjust = .5, face = "bold", size = 7, lineheight = .9, margin = margin(0, 4, 0, 0)),
          axis.text.y = element_text(size = 7, hjust = 1, margin = margin(0, 1, 0, 0)),
          panel.spacing.x = unit(3, "pt"), plot.margin = margin(3, 4, 2, 2), plot.background = element_blank())
  p <- (pAh / pAs + plot_layout(heights = unit(c(1, A_STRIP_H), c("null", "in")))) & theme(plot.background = element_blank())
  list(plot = p + plot_annotation(theme = theme(plot.background = element_blank())), lims = LIMS)
}

# ---- B: the regulon screen (one row per selected factor, one dot per regulon) ----------------------------------------------
# Three columns: SC-islets (SC-SIC minus SC-beta, median over SC-islet studies), intestine (EC minus other enteroendocrine,
# median over intestinal studies; Elmentaite 2020 + 2021 the mean of its parts) and adult islets (beta minus other endocrine,
# mean over adult islet studies; axis broken at GG_AD_CUT). Reproducible regulons are colored by direction, the rest grey.
GG_W <- 7.40; GG_H <- 3.26 + .50                    # slot size (in)
GG_BLOCK_TOP <- .40                                 # top of the three-column plot block in the slot
GG_PT <- 7.5                                        # key size (pt)
GG_KEY_GAP <- -.01                                  # key text's top below the block's foot (in)
GG_STUDY <- c("ELMENTAITE", "HICKEY_2023", "EGOZI_2021")
GG_HALF <- c("Targets higher\nin SC-SIC", "Targets lower\nin SC-SIC")
GG_COLS <- c(pass_ec = unname(POP_COL["SC-SIC"]), pass_beta = unname(POP_COL["SC-beta"]), fail = "grey45")
GG_PTITLE <- c(sc = "SC-islets", gut = "Intestine", ad = "Adult islets")
GG_XLAB <- c(sc = "SC-SIC minus SC-beta,\nmedian over SC-islet studies", gut = "EC minus other enteroendocrine,\nmedian over intestinal studies",
             ad = "Beta minus other endocrine,\nmedian over adult islet studies")
GG_AD_CUT <- -.075                                  # adult islets: axis break here; values below it go to a condensed stub on the left
GG_KEY_HEAD <- "Reproducible in SC-islets and intestine:"
GG_KEYS <- list(list("higher in SC-SIC", GG_COLS[["pass_ec"]]), list("lower in SC-SIC", GG_COLS[["pass_beta"]]), list("not reproducible", GG_COLS[["fail"]]))
GG_PM <- c("(+) activated targets", "(−) repressed targets")   # second key row
GG_KEY2_DY <- .14                                   # the second key row's top below the first's (in)
gg_tb <- function(x) !is.na(x) & as.character(x) %in% c("TRUE", "True", "true")
gg_th <- function() theme_minimal(base_family = "Arial", base_size = 7) + theme(
  text = element_text(color = "black", size = 7), axis.text = element_text(color = "black", size = 7),
  axis.title = element_text(color = "black", size = 7, lineheight = .95), plot.title = element_text(color = "black", size = 7.5, face = "bold", margin = margin(0, 0, 2, 0)),
  panel.grid.minor = element_blank(), panel.grid.major.y = element_blank(), panel.grid.major.x = element_line(color = "#E4E4E4", linewidth = .2),
  axis.line.x = element_line(color = "black", linewidth = .35), axis.ticks.x = element_line(color = "black", linewidth = .3),
  axis.ticks.length = unit(1.3, "pt"), plot.margin = margin(2, 6, 2, 2), legend.position = "none", plot.title.position = "panel",
  strip.text = element_blank(), strip.placement = "outside", panel.spacing.y = unit(6, "pt"))
plot_B <- function() {
  GGF <- read.delim(result_path("regulons/factors.tsv"), check.names = FALSE, stringsAsFactors = FALSE, na.strings = c("NA", ""), fileEncoding = "UTF-8")
  GGF <- GGF[gg_tb(GGF$selected), ]
  GGV <- read.delim(result_path("regulons/program_values.tsv"), check.names = FALSE, stringsAsFactors = FALSE, na.strings = c("NA", ""), fileEncoding = "UTF-8")
  GGV$pass <- gg_tb(GGV$program_pass)
  # rows: two groups by direction, factors alphabetical within a group
  GGF$half <- factor(ifelse(GGF$direction == "SC-SIC higher", GG_HALF[1], GG_HALF[2]), GG_HALF)
  GGF$tf <- sub("\\(.*$", "", GGF$factor)
  GGF <- GGF[order(GGF$half, GGF$tf, GGF$factor), ]
  GGF$expr <- sprintf('italic("%s")~"(%s)"', GGF$tf, ifelse(grepl("\\(\\+\\)$", GGF$factor), "+", "−"))   # italic symbol, roman (+)/(−)
  GGV$half <- GGF$half[match(GGV$factor, GGF$factor)]
  GGV <- GGV[order(match(GGV$factor, GGF$factor), GGV$network), ]
  # a factor's regulons spread over +-0.30 of its row
  GGV$off <- ave(seq_len(nrow(GGV)), GGV$factor, FUN = function(i) if (length(i) == 1) 0 else seq(-.30, .30, length.out = length(i)))
  GGV$y <- as.numeric(factor(GGV$factor, rev(GGF$factor))) + GGV$off
  GGV$col <- ifelse(!GGV$pass, "fail", ifelse(GGV$direction == "SC-SIC higher", "pass_ec", "pass_beta"))
  NL <- data.frame(factor = GGF$factor, half = GGF$half, y = as.numeric(factor(GGF$factor, rev(GGF$factor))), expr = GGF$expr)
  # grey band on every second row; the lowest row gets 0.2 row of padding below it
  NL$band <- (seq_len(nrow(NL)) %% 2L) == 0L
  YMIN <- min(NL$y); pad <- function(y) ifelse(y == YMIN, .2, 0)
  BANDS <- transform(NL[NL$band, ], ymin = y - .5 - pad(y), ymax = y + .5)
  YB <- rbind(transform(NL, y = y - .5 - pad(y)), transform(NL, y = y + .5))
  gg_col <- function(d, key, lab_y, lims = NULL, brk = NULL, part = "") { m <- max(abs(d$x)) * 1.07
    if (is.null(lims)) lims <- c(-m, m)
    d <- d[d$x >= lims[1] & d$x <= lims[2], ]
    br <- if (is.null(brk)) scales::breaks_pretty(4)(lims) else brk; br <- br[is.finite(br) & br >= lims[1] & br <= lims[2]]
    gl <- merge(BANDS[, c("half", "ymin", "ymax")], data.frame(x = br))          # x grid lines redrawn over the bands
    ggplot() + geom_rect(data = BANDS, aes(xmin = -Inf, xmax = Inf, ymin = ymin, ymax = ymax), inherit.aes = FALSE, fill = BAND, color = NA) +
      geom_segment(data = gl, aes(x = x, xend = x, y = ymin, yend = ymax), inherit.aes = FALSE, color = "#E4E4E4", linewidth = .2) +
      (if (lims[1] < 0 && lims[2] > 0) geom_vline(xintercept = 0, color = "#777777", linewidth = .3) else list()) +
      geom_point(data = d[!d$pass, ], aes(x = x, y = y), shape = 16, size = .85, color = GG_COLS[["fail"]]) +
      (if (any(d$pass)) geom_point(data = d[d$pass, ], aes(x = x, y = y, color = col), shape = 16, size = .95) else list()) +
      facet_grid(half ~ ., scales = "free_y", space = "free_y", switch = "y") +
      geom_blank(data = YB, aes(x = 0, y = y), inherit.aes = FALSE) +
      scale_y_continuous(breaks = NL$y, labels = if (lab_y) parse(text = NL$expr) else NULL, expand = c(0, 0)) +
      (if (any(d$pass)) scale_color_manual(values = GG_COLS) else list()) +
      scale_x_continuous(limits = lims, breaks = br, labels = function(v) sub("-", "−", format(v, drop0trailing = FALSE)), expand = c(0, 0), oob = scales::oob_keep) +
      labs(x = if (part == "stub") NULL else GG_XLAB[[key]], y = NULL, title = if (part == "main") " " else GG_PTITLE[[key]]) + gg_th() +
      # axis break: one "//" on the axis line in the gap between the stub and the main axis
      (if (part == "stub") coord_cartesian(clip = "off") else if (part == "main") list(coord_cartesian(clip = "off"),
        geom_text(data = data.frame(half = factor(GG_HALF[2], GG_HALF), x = lims[1] - diff(lims) * .012, y = YMIN - .5 - pad(YMIN)),
                  aes(x = x, y = y, label = "//"), inherit.aes = FALSE, size = 7 / .pt, family = "Arial", vjust = .5)) else list()) +
      (if (part == "stub") theme(plot.margin = margin(2, 0, 2, 2)) else if (part == "main") theme(plot.margin = margin(2, 6, 2, 0)) else theme()) +
      (if (lab_y) theme(strip.text.y.left = element_text(angle = 0, hjust = 1, vjust = .5, face = "bold", size = 7, lineheight = .95, margin = margin(0, 5, 0, 0)),
                        axis.text.y = element_text(size = 7, hjust = 1, margin = margin(0, 2, 0, 0))) else theme()) }
  d1 <- GGV; d1$x <- d1$sc_signed_median; d3 <- GGV; d3$x <- d3$adult_mean
  # intestine: one dot per regulon, the median of its intestinal study values (Elmentaite 2020 + 2021 the mean of its parts)
  d2 <- GGV; d2$x <- apply(as.matrix(GGV[, paste0("gut_mean__", GG_STUDY)]), 1, function(v) median(v[is.finite(v)]))
  # adult islets (not used for the colors): axis broken at GG_AD_CUT; the values below it in a condensed stub
  p3 <- if (any(d3$x < GG_AD_CUT, na.rm = TRUE)) {
    lo <- min(d3$x, na.rm = TRUE) * 1.05; hi <- max(d3$x, na.rm = TRUE) * 1.07
    stub <- gg_col(d3, "ad", FALSE, lims = c(lo, GG_AD_CUT), brk = -.3, part = "stub")
    main <- gg_col(d3, "ad", FALSE, lims = c(GG_AD_CUT, hi), brk = c(-.05, 0, .05), part = "main")
    wrap_plots(stub, main, nrow = 1, widths = c(.16, .84)) & theme(plot.margin = margin(2, 0, 2, 0)) }
    else gg_col(d3, "ad", FALSE)
  wrap_plots(list(gg_col(d1, "sc", TRUE), gg_col(d2, "gut", FALSE), p3), nrow = 1, widths = c(1, 1, 1))
}
# text width (in): the wider of systemfonts' and the device's measure
gg_tw <- function(s, pt = 7) {
  sf <- systemfonts::string_width(s, family = "Arial", size = pt, res = 600, weight = "normal") / 600
  dv <- vapply(s, function(x) { pushViewport(viewport(gp = gpar(fontfamily = "Arial", fontsize = pt, fontface = "plain")))
    v <- convertWidth(stringWidth(x), "in", valueOnly = TRUE); popViewport(); v }, 0, USE.NAMES = FALSE)
  pmax(sf, dv) }
# draw panel B into the current viewport (its GG_W x GG_H slot): letter, plot block, and the two key rows centred under it.
# dy: height freed above the block; extra: added block height; grow: the block's top moves up this much (in)
gg_render <- function(block, dy, extra, grow) {
  heading("B", .06, .03, vp_h = GG_H)
  hw <- gg_tw(GG_KEY_HEAD, GG_PT)
  kw <- hw + .10 + sum(vapply(GG_KEYS, function(k) .14 + gg_tw(k[[1]], GG_PT) + .16, 0))
  ky <- GG_H - .04 + extra - dy + GG_KEY_GAP
  x <- (GG_W - (kw - .16)) / 2
  text_at(GG_KEY_HEAD, x, ky, vp_h = GG_H); x <- x + hw + .10
  for (k in GG_KEYS) {
    grid.points(unit(x + .05, "in"), unit(GG_H - (ky + .005) - .045, "in"), pch = 16, size = unit(1.5, "mm"), gp = gpar(col = k[[2]]))
    text_at(k[[1]], x + .14, ky, vp_h = GG_H); x <- x + .14 + gg_tw(k[[1]], GG_PT) + .16 }
  pushViewport(viewport(x = unit(0, "in"), y = unit(GG_H - (GG_BLOCK_TOP - dy - grow), "in"), width = unit(GG_W, "in"), height = unit(GG_H - GG_BLOCK_TOP - .04 + extra + grow, "in"), just = c("left", "top")))
  grid.draw(patchworkGrob(block)); popViewport()
  pw <- gg_tw(GG_PM, GG_PT); x <- (GG_W - (sum(pw) + .3)) / 2
  for (i in seq_along(GG_PM)) { text_at(GG_PM[i], x, ky + GG_KEY2_DY, vp_h = GG_H); x <- x + pw[i] + .3 }
}

# ---- C: the CDX2 knockout within each cell type ------------------------------------------------------------------------
# full page width, 7 pt text, its key drawn under it (no ggplot legends)
C_BAR <- "grey82"                                   # random-range bars
C_REF <- "#777777"
theme_c <- function() {
  gl <- element_line(color = "#E4E4E4", linewidth = .2)
  theme_minimal(base_family = "Arial", base_size = 7) +
    theme(text = element_text(color = "black", size = 7), axis.text = element_text(color = "black", size = 7),
          axis.title = element_text(color = "black", size = 7), strip.text = element_text(color = "black", size = 8, face = "bold", hjust = 0),
          plot.title = element_text(color = "black", size = 8, face = "bold"), plot.subtitle = element_text(color = "black", size = 7),
          legend.text = element_text(color = "black", size = 7), legend.title = element_text(color = "black", size = 7),
          legend.position = "none", legend.key.size = unit(.12, "in"), panel.grid.minor = element_blank(),
          panel.grid.major.x = gl, panel.grid.major.y = element_blank(),
          axis.line.x = element_line(color = "black", linewidth = .35), axis.ticks.x = element_line(color = "black", linewidth = .3),
          axis.ticks.length = unit(1.3, "pt"), panel.spacing = unit(.16, "in"), plot.margin = margin(2, 6, 2, 2),
          plot.background = element_blank())
}
plot_C <- function() {
  # knockout minus wild type within each cell type (no SC-beta subtraction), the 12,615 nuclei of the study's GEO deposit;
  # grey bars: middle 95% of 500 expression-matched random gene sets
  E <- read_table("cdx2_knockout/group_test/group_effects.tsv")
  SETD <- c("G2 serotonin core" = "Four EC\nmarkers", "G2 other" = "Other shared\nEC genes (G2)", "G3" = "Separate SC-SIC\nprogram (G3)")
  E <- E[E$variant == "all transferred" & E$set %in% names(SETD) & E$population %in% c("SC_EC", "SC_beta", "SC_alpha"), ]
  stopifnot(nrow(E) == 18)
  k <- ifelse(grepl("per-cell", E$readout), 100, 1)                     # per-cell detection as percentage points
  E$v <- E$ko_minus_wt * k; E$lo <- E$random_q025 * k; E$hi <- E$random_q975 * k
  PCD <- "Per-cell detection\n(points of group genes)"
  E$rd <- factor(ifelse(grepl("per-cell", E$readout), PCD, "Pseudobulk\n(log2)"), c("Pseudobulk\n(log2)", PCD))
  E$pop <- c(SC_EC = "SC-SIC", SC_beta = "SC-beta", SC_alpha = "SC-alpha")[E$population]
  E$grp <- factor(SETD[E$set], SETD)               # one facet row per gene set (bold label at left)
  E$y <- c("SC-SIC" = 1, "SC-beta" = 2, "SC-alpha" = 3)[E$pop]
  # dots in each cell type's color, filled outside the random range, a ring inside it
  E$out <- factor(ifelse(E$v < E$lo | E$v > E$hi, "Outside random range", "Inside"), c("Outside random range", "Inside"))
  ggplot(E) +
    geom_vline(xintercept = 0, colour = C_REF, linewidth = .25) +
    geom_linerange(aes(y = y, xmin = lo, xmax = hi, linetype = "Middle 95% of 500 random sets"), colour = C_BAR, linewidth = 1.6) +
    geom_point(data = E[E$out == "Inside", ], aes(y = y, x = v), shape = 21, fill = "white", colour = POP_COL[E$pop[E$out == "Inside"]], size = 1.44, stroke = .8) +
    geom_point(data = E[E$out != "Inside", ], aes(y = y, x = v), shape = 21, fill = POP_COL[E$pop[E$out != "Inside"]], colour = POP_COL[E$pop[E$out != "Inside"]], size = 2.1, stroke = .3) +
    facet_grid(grp ~ rd, scales = "free_x", switch = "y") +
    scale_y_reverse(breaks = 1:3, labels = c("SC-SIC", "SC-beta", "SC-alpha"), limits = c(3.5, .5), expand = c(0, 0)) +
    scale_x_continuous(expand = expansion(mult = .09)) +
    scale_linetype_manual(values = 1, name = NULL, guide = guide_legend(order = 3, override.aes = list(linewidth = 1.6, colour = C_BAR))) +
    labs(x = expression(italic("CDX2")~"knockout minus wild type, within each cell type"), y = NULL) +
    theme_c() +
    theme(strip.text.x = element_text(hjust = .5, size = 7, face = "bold"), strip.placement = "outside",
          strip.text.y.left = element_text(angle = 0, hjust = 1, size = 7, face = "bold"), panel.spacing.y = unit(.05, "in"),
          legend.margin = margin(0, 0, 0, 0), legend.location = "plot", legend.key.spacing.x = unit(.1, "in"),
          panel.border = element_rect(fill = NA, colour = "grey80", linewidth = .3),
          axis.text.y = element_text(colour = POP_COL, face = "bold"),
          panel.grid.major.x = element_blank(), panel.grid.minor = element_blank(),
          legend.key = element_rect(fill = "white", colour = NA))
}

# ---- page layout (inches from the page's top-left) ----------------------------------------------------------------------
# Positions are the published page's exact values; where a short decimal is not exactly that value, the sum that gives it
# is written instead.
A_DOWN <- .29; A_KEY_UP <- .19; A_HM_GROW <- .30; A_STRIP_GROW <- .26
A_STRIP_H <- .45 + A_DOWN + A_STRIP_GROW           # height of 5A's per-study strip (in)
A_LAB <- "SC-SIC minus SC-beta\n(log2), per study"
FREED <- c(A = .16, G = .16); A_EXTRA <- .24; A_DROP <- 79 / 600; G_EXTRA <- .1108333
LAY <- list(y1 = .46, h1 = 2.16); KEY_DY <- .06
GROW <- c(A = 147, G = 72) / 600
Y1 <- LAY$y1 - FREED[["A"]] - GROW[["A"]]           # A's top
H1 <- LAY$h1 + GROW[["A"]] + A_EXTRA + A_DOWN - A_DROP + A_HM_GROW + A_STRIP_GROW
Y_E <- Y1 + H1 + .03                                # B's letter line: the slot's top is A's foot
Y_FG <- Y_E - .03 + GG_H + .10 + GG_KEY2_DY         # C's letter line
H_FG <- 1.49                                        # C's height

# Draw one panel ("A", "B" or "C") at its place on the page.
render_figure5 <- function(only) {
  if (only == "A") {
    A <- plot_A()
    heading("A", .10, .10)
    draw_at(patchworkGrob(A$plot), .04, Y1, 7.40, H1)
    grid.text(A_LAB, x = unit(1.50, "in"), y = unit(H - (Y1 + H1 - 1.01 - A_STRIP_GROW / 2), "in"), just = c("right", "centre"),
              gp = gpar(fontfamily = "Arial", fontsize = 7, fontface = "bold", lineheight = .95))
    colourbar(.30, Y1 + H1 - A_KEY_UP - .18, .9, A$lims, "Mean log2(1 + CPM)")
    row_keys(list(list(type = "pt", label = "SC-SIC higher", col = POP_COL[["SC-SIC"]]), list(type = "dash", label = "Mean", col = "black", lwd = .2 * ggplot2::.pt * 2.5)),
             .10, 1.7, Y1 + H1 - A_KEY_UP - .66, centre = FALSE, gap = .12)
    row_keys(list(list(type = "pt", label = "SC-beta higher", col = POP_COL[["SC-beta"]])), .10, 1.7, Y1 + H1 - A_KEY_UP - .54, centre = FALSE, gap = .12)
    row_keys(list(list(type = "pt", label = "No RNA at this depth", col = "grey60")), .10, 1.7, Y1 + H1 - A_KEY_UP - .42, centre = FALSE, gap = .12)
  }
  if (only == "B") {
    block <- (plot_B() & theme(plot.background = element_blank())) + plot_annotation(theme = theme(plot.background = element_blank()))
    pushViewport(viewport(x = unit(.04, "in"), y = unit(H - (Y_E - .03), "in"), width = unit(GG_W, "in"), height = unit(GG_H, "in"), just = c("left", "top")))
    gg_render(block, FREED[["G"]], G_EXTRA, GROW[["G"]])
    popViewport()
  }
  if (only == "C") {
    heading("C", .10, Y_FG)
    draw_at(ggplotGrob(plot_C() + theme(plot.background = element_blank())), .04, Y_FG, W - .08, H_FG)
    K1 <- .3 * ggplot2::.stroke / 2
    row_keys(list(list(type = "pt", pch = 21, label = "Outside random range", col = "grey35", fill = "grey35", size = 1.9, plwd = K1),
                  list(type = "pt", pch = 21, label = "Inside", col = "grey35", fill = "white", size = 1.57, plwd = .8 * ggplot2::.stroke / 2),
                  list(type = "ln", label = "Middle 95% of 500 random sets", col = C_BAR, lwd = 1.6 * ggplot2::.pt)), .04, W - .08, Y_FG + H_FG + KEY_DY, gap = .12)
  }
}
