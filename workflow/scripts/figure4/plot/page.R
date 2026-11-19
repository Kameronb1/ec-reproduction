# The Figure 4 page, shared by its panel scripts: the tables, the eleven rows, the measured layout of every panel and one
# drawing function per panel. Each panel script sources this file and draws its own part.
#   A  share of mature endocrine cells (published labels) of each graft against its in vitro comparison samples
#   B  per row, the share of SC-SICs and SC-beta cells with a four-marker score at or above the cutoff, start (open
#      circle) to graft or extended culture (arrowhead), samples weighted equally; the four EC markers' expression change
#   C  median expression difference from primary intestinal EC on the gene groups G1, G2, G3, start to end; pale bar =
#      middle half of the group's genes at the end
#   D  promoter accessibility of each group against the typical gene of the same cell type: Zhu 2023 (days 21, 32, 39)
#      and Augsornworawat 2023 (weeks 2-4, months 6 and 12, graft median) on one axis
#   E  distal accessibility change after grafting against 10,000 random gene sets
# Positions are in inches from the page's top-left corner. The layout is measured on a cairo device of the page's size.

# ---- page table -------------------------------------------------------------------------------------------------------
LT <- c(A = .10)                       # top of each panel letter; the others follow from the content above them
RIGHT <- 7.44
RU <- .115                             # row pitch of B and C (in per row unit)
B_H <- 2.34                            # proportions plot height
GAP_ROWS <- .10                        # from a key's last line to the next panel letter
B_PITCH <- .097                        # line pitch of the proportions strip labels
DE_H <- 1.82                           # D and E plot height, strips to the marks under D's axis
E_W <- 2.15                            # E's width
SEC_PT <- 7.5                          # section headers: 7.5 pt bold, flush left at .06 in; study names 7 pt bold, indented

# ---- text -------------------------------------------------------------------------------------------------------------
SECTION <- c(graft = "Grafts (in mice)", cult = "Extended culture (in vitro)")
GLAB <- c(G1 = "G1: EC genes low in SC-SICs", G2 = "G2: shared EC genes", G3 = "G3: separate SC-SIC program")
GLAB2 <- c(G1 = "G1: EC genes\nlow in SC-SICs", G2 = "G2: shared\nEC genes", G3 = "G3: separate\nSC-SIC program")
A_AXIS <- "Median expression difference from primary intestinal EC (log2)"
A_KEY <- c(ref = "primary intestinal EC", bar = "middle half of genes")
B_TEXT <- c(x = "Share of mature endocrine cells (%, published labels)", header = "Graft (mature endocrine cells)",
            range = "in vitro samples: range, mean", graft = "graft", low = "graft with < 50 mature endocrine cells",
            n1 = "Samples", n2 = "start → end")
TAB_HEAD <- c(study = "Study", line = "Line and end point")
C_HEAD <- c(C1 = "Four-marker score at or above the cutoff", C2 = "Four EC markers")
C_AXIS <- c(C1 = "Cells at or above the cutoff (%)", C2 = "Expression change (log2)")
MARKERS_NOTE <- c(lead = "Four EC markers: ", genes = "TPH1, SLC18A1, FEV, LMX1A")

# ---- the 11 comparisons of B and C ---------------------------------------------------------------------------------------
ROWS <- data.frame(
  id = c("A-HUES8", "A-WS4", "B-0-M1", "B-0-M3", "B-3-M6-CTRL", "B-3-M6-STZ", "H-10w", "P-TXP", "B-w0-w6", "P-w2-m6", "P-w2-m12"),
  study = c("Augsornworawat 2020", "Augsornworawat 2020", rep("Balboa 2022", 4), "Du 2022", "Augsornworawat 2023",
            "Balboa 2022", "Augsornworawat 2023", "Augsornworawat 2023"),
  short = c("HUES8, 6-mo graft", "WS4corr, 6-mo graft", "H1, 1-mo graft", "H1, 3-mo grafts (pooled)", "H1, 6-mo graft (control host)",
            "H1, 6-mo graft (STZ host)", "hCiPSC, 10-wk graft", "HUES8 (nuclei), 6-mo graft", "H1, wk 6", "HUES8 (nuclei), mo 6", "HUES8 (nuclei), mo 12"),
  y = c(1:8, 10.25, 11.25, 12.25), stringsAsFactors = FALSE)
GRAFT <- ROWS$id[1:8]; CULT <- ROWS$id[9:11]
Y_SEC <- c(graft = 0, cult = 9.125)                          # the two section-header lines, in row units
YR <- c(max(ROWS$y) + .6, -.6); BAND_ROWS <- ROWS$y[seq(1, nrow(ROWS), by = 2)]
rowy <- function(id, yy = ROWS$y) yy[match(id, ROWS$id)]
sep <- function(at = 9.55) geom_hline(yintercept = at, color = "#9A9A9A", linetype = "dotted", linewidth = .3)   # grafts | extended culture
yrows <- function(yr = YR) scale_y_reverse(limits = yr, breaks = NULL, expand = c(0, 0))   # row labels drawn on the page
YA <- list(y = ROWS$y, yr = YR, sep = 9.55); YC <- YA
YG <- list(y = ROWS$y, yr = c(YR[1], .4), sep = 9.55)        # C: "Grafts (in mice)" on the facet-header line
RS <- read_table("row_samples.tsv")                          # each comparison's start and end samples
REF <- setNames(strsplit(RS$ref_samples, ";"), RS$id); TEST <- setNames(strsplit(RS$test_samples, ";"), RS$id)
nsamp <- function(L) vapply(L, function(s) length(unlist(strsplit(s, "+", fixed = TRUE))), 0L)   # a pooled end ("a+b+c") counts its samples
NCOUNT <- setNames(sprintf("%d → %d", nsamp(REF[ROWS$id]), nsamp(TEST[ROWS$id])), ROWS$id)   # samples, start -> end

# ---- C: median difference from primary intestinal EC on G1, G2, G3 --------------------------------------------------------
EN <- read_table("E_group_gaps.tsv")
OFF3 <- c(SC_SIC = -.31, SC_BETA = 0, SC_ALPHA = .31)
PE <- expand.grid(population = names(POP_COL), id = ROWS$id, group = names(GLAB), stringsAsFactors = FALSE)[, c("group", "id", "population")]
k <- match(paste(PE$group, PE$id, PE$population), paste(EN$group, EN$id, EN$population)); PE$from <- EN$from[k]; PE$to <- EN$to[k]
PE$y <- rowy(PE$id) + OFF3[PE$population]; PE$panel <- factor(GLAB[PE$group], GLAB)
ES <- read_table("E_end_spread.tsv"); ES$y <- rowy(ES$id) + OFF3[ES$population]; ES$panel <- factor(GLAB[ES$group], GLAB)
stopifnot(nrow(ES) == sum(!is.na(PE$to)), all(paste(ES$group, ES$id, ES$population) %in% paste(PE$group, PE$id, PE$population)[!is.na(PE$to)]))

# ---- A: each graft's share of mature endocrine cells against its in vitro comparison samples -----------------------------
ab <- read_table("f4b_composition.tsv")
B_STRIP <- list(list(c(b = "Augsornworawat 2020"), c(p = "HUES8, from stage 6")), list(c(b = "Augsornworawat 2020"), c(p = "WS4corr, from stage 6")),
                list(c(b = "Balboa 2022"), c(p = "H1, from stage 7 wk 0")),
                list(c(b = "Balboa 2022"), c(p = "H1, from stage 7 wk 3")),
                list(c(b = "Du 2022", p = ", hCiPSC"), c(p = "stage 6 (vs day 2)")),
                list(c(b = "Augsornworawat 2023"), c(p = "HUES8 (nuclei)"), c(p = "from stage 6 wk 2")))
names(B_STRIP) <- vapply(B_STRIP, function(L) paste(vapply(L, paste, "", collapse = ""), collapse = "\n"), "")   # b: bold, p: plain
SRC <- data.frame(id = GRAFT,
                  study = c("Augsornworawat2020", "Augsornworawat2020", "Balboa2022", "Balboa2022", "Balboa2022", "Balboa2022", "Du2022", "Augsornworawat2023"),
                  graft = c("M6_graft_HUES8", "M6_graft_WS4corr", "M1_graft", "M3_graft_pooled", "M6_graft_CTRL", "M6_graft_STZ", "10_wpt", "M6_graft"),
                  source = c("S6_in_vitro", "WS4_baseline", "S7w0", "S7w0", "S7w3", "S7w3", "S6D2", "sc_islet_1|sc_islet_2|S6_week2"),
                  strip = names(B_STRIP)[c(1, 2, 3, 3, 4, 4, 5, 6)], stringsAsFactors = FALSE)
GC <- c("SC-SIC", "SC-alpha", "SC-beta")
PBs <- do.call(rbind, lapply(seq_len(nrow(SRC)), function(i) {
  s <- ab[ab$study == SRC$study[i] & ab$context %in% strsplit(SRC$source[i], "|", fixed = TRUE)[[1]] & ab$category %in% GC, ]
  gr <- ab[ab$study == SRC$study[i] & ab$context == SRC$graft[i] & ab$category %in% GC, ]
  do.call(rbind, lapply(GC, function(cg) { sv <- s$pct[s$category == cg]; gv <- gr[gr$category == cg, ]
    data.frame(strip = SRC$strip[i], context = SRC$graft[i], row_lab = gv$row_lab, mature = gv$mature, category = cg, pct = gv$pct,
               src_mean = mean(sv), src_lo = min(sv), src_hi = max(sv), stringsAsFactors = FALSE) })) }))
short_lab <- function(x) { x <- sub(" \\(from S7 wk [03]\\)", "", x); x <- sub(" \\((HUES8|WS4corr)\\)", "", x)
  x <- sub("6-mo graft, ", "6-mo, ", x); x <- sub("10-week graft", "10-wk graft", x); x <- sub(" grafts, pooled", ", pooled", x); x <- sub(" graft", "", x); gsub("  #", " #", x) }
PBs$glab <- sprintf("%s (%s)", short_lab(PBs$row_lab), trimws(format(PBs$mature, big.mark = ",")))
PBs$strip <- factor(PBs$strip, unique(SRC$strip)); PBs$category <- factor(PBs$category, c("SC-SIC", "SC-beta", "SC-alpha"))
PBs$glab <- factor(PBs$glab, rev(unique(PBs[order(PBs$strip, PBs$context, PBs$row_lab, method = "radix"), "glab"])))
PBs$low <- PBs$mature < 50                                   # fewer than 50 mature endocrine cells: open arrowhead

# ---- B: share of cells at or above the four-marker cutoff, start to end; the four markers' expression change -----------
CPCT <- read_table("four_marker_row_shares.tsv")             # per row end, samples weighted equally (ec_calls.R)
PC <- do.call(rbind, lapply(split(CPCT, paste(CPCT$id, CPCT$population)), function(x)
  data.frame(id = x$id[1], population = x$population[1], from = x$pct[x$end == "from"], to = x$pct[x$end == "to"], stringsAsFactors = FALSE)))
PC <- PC[order(match(PC$id, ROWS$id), match(PC$population, c("SC_SIC", "SC_BETA"))), ]; rownames(PC) <- NULL
stopifnot(nrow(PC) == 2L * nrow(ROWS))
C_OFF <- c(SC_SIC = -.18, SC_BETA = .18)
PC$y <- rowy(PC$id, YC$y) + C_OFF[PC$population]
XL1 <- c(-3, 103)                                            # 0 to 100 %, with room for a circle or head at 0 or 100
CB <- read_table("C_bars.tsv"); CB <- CB[CB$id %in% ROWS$id & CB$population %in% c("SC_SIC", "SC_BETA"), ]
CB$y <- rowy(CB$id, YC$y) + c(SC_SIC = -.18, SC_BETA = .18)[CB$population]; CB$col <- POP_COL[CB$population]
BV <- CB$four_marker
XL2 <- c(min(-3.9, floor(min(BV) * 2) / 2 - .2), max(1.2, ceiling(max(BV) * 2) / 2 + .2))

# ---- D, E -----------------------------------------------------------------------------------------------------------------
source(file.path(dirname(PLOT_SCRIPT), "panels_de.R"))

# ---- plots (marks = FALSE: the same plot without its marks, used to measure the layout) ------------------------------
colplot <- function(xlim, brk, lab, xline = NULL, minor = waiver()) ggplot() + sep(YC$sep) +
  (if (!is.null(xline)) geom_vline(xintercept = xline, color = "#777777", linewidth = .3)) +
  scale_x_continuous(limits = xlim, breaks = brk, labels = lab, minor_breaks = minor, expand = c(0, 0)) + yrows(YC$yr) + labs(x = NULL, y = NULL, title = NULL) + th()
prop_plot <- function(marks) {
  p <- ggplot(PBs, aes(y = glab)) + geom_linerange(aes(xmin = src_lo, xmax = src_hi), linewidth = 2.6, colour = "grey75") +
    geom_blank(aes(x = src_mean)) + geom_blank(aes(x = pct))
  if (marks) p <- p +                                         # shaft, head (open when < 50 cells), then the start circle on top
    geom_segment(data = PBs[!PBs$bare, ], aes(x = src_mean, xend = shaft_end, yend = glab, colour = category), linewidth = .45, lineend = "butt") +
    geom_segment(data = PBs[!PBs$low, ], aes(x = hbase, xend = pct, yend = glab, colour = category), linewidth = .45, arrow = ARW) +
    geom_segment(data = PBs[PBs$low, ], aes(x = hbase, xend = pct, yend = glab, colour = category), linewidth = .30, arrow = ARW, arrow.fill = "white") +
    geom_point(data = PBs[!PBs$bare, ], aes(x = src_mean, colour = category), shape = 21, fill = "white", size = PT_SIZE, stroke = PT_STROKE) +
    scale_colour_manual(values = GCOL) + scale_fill_manual(values = GCOL, guide = "none")
  p + facet_grid(strip ~ category, scales = "free_y", space = "free_y", switch = "y") +
    scale_x_continuous(limits = c(0, 82), breaks = c(0, 25, 50, 75), expand = c(0, 0)) +
    labs(x = B_TEXT[["x"]], y = NULL) + th() +
    theme(strip.placement = "outside", strip.text.y.left = element_blank(), strip.text.x = element_text(size = 7.5, face = "bold"),
          panel.spacing.y = unit(2.5, "pt"), panel.spacing.x = unit(10, "pt"), panel.border = element_rect(fill = NA, color = "#BDBDBD", linewidth = .25),
          axis.text.y = element_text(size = 7)) }
gg_plot <- function(marks) ggplot() + sep() + geom_vline(xintercept = 0, color = "#555555", linetype = "dashed", linewidth = .4) +
  (if (marks) geom_segment(data = ES, aes(x = q25, xend = q75, y = y, yend = y, color = population), linewidth = 1.6, alpha = .3, lineend = "butt")
   else list(geom_blank(data = ES, aes(x = q25, y = y)), geom_blank(data = ES, aes(x = q75, y = y)))) + scale_color_manual(values = POP_COL) +
  lapply(names(POP_COL), function(k) { d <- PE[PE$population == k, ]
    if (marks) arrows_ef(d, POP_COL[[k]]) else list(geom_blank(data = d, aes(x = from, y = y)), geom_blank(data = d, aes(x = to, y = y))) }) +
  facet_grid(. ~ panel, scales = "free_x", space = "free_x", labeller = as_labeller(setNames(GLAB2, GLAB))) + yrows(YG$yr) +
  labs(x = A_AXIS, y = NULL) + th() + theme(panel.spacing.x = unit(14, "pt")) + theme(plot.margin = margin(2, 4, 0, GG_LABEL_PT))

# ---- measured layout ----------------------------------------------------------------------------------------------------
MEASURE_PNG <- tempfile(fileext = ".png"); open_png(MEASURE_PNG)
SEC_W <- max(vapply(SECTION, str_in, 0, size = SEC_PT, face = "bold"))
# A: strip labels from x = .10, then the sample-count column, then the plot
B_X <- .10; B_LINES <- unlist(B_STRIP, recursive = FALSE)
run_w <- function(r) sum(vapply(seq_along(r), function(i) str_in(r[[i]], 7, if (names(r)[i] == "b") "bold" else "plain"), 0))
B_ASC <- max(unlist(lapply(B_LINES, function(r) mapply(function(s, b) asc_in(s, 7, if (b) "bold" else "plain"), r, names(r) == "b"))))
B_DSC <- max(unlist(lapply(B_LINES, function(r) mapply(function(s, b) dsc_in(s, 7, if (b) "bold" else "plain"), r, names(r) == "b"))))
B_STRIP_W <- max(vapply(B_LINES, run_w, 0))
B_N_W <- max(str_in(B_TEXT[["n1"]], 7), str_in("start  → end", 7), vapply(NCOUNT, str_in, 0, size = 7))
B_N_RIGHT <- B_X + B_STRIP_W + .16 + B_N_W
B_PLOT_X <- B_N_RIGHT + .06; B_PLOT_W <- RIGHT - B_PLOT_X
B_Y0 <- LT[["A"]] + .07 + 7.5 / 72
gB0 <- ggplotGrob(prop_plot(FALSE)); nnB <- unitType(gB0$widths) == "null"
B_IPU <- (B_PLOT_W - sum(convertWidth(gB0$widths[!nnB], "in", valueOnly = TRUE))) / 3 / diff(ggplot_build(prop_plot(FALSE))$layout$panel_params[[1]]$x.range)
B_KEY_TOP <- B_Y0 + B_H + .04; B_KEY_PITCH <- .135
LT[["B"]] <- B_KEY_TOP + B_KEY_PITCH * (1 + any(PBs$low)) + GAP_ROWS - .02
# B: label table (study at .14, line and end point, samples right-aligned), then the share arrows and the bars 3 : 2
TB_X_STUDY <- .14
STUDY_W <- max(vapply(unique(ROWS$study), str_in, 0, size = 7, face = "bold"), str_in(TAB_HEAD[["study"]], 7))
TB_X_LINE <- TB_X_STUDY + STUDY_W + .12
LINE_W <- max(vapply(ROWS$short, str_in, 0, size = 7), str_in(TAB_HEAD[["line"]], 7))
TB_N_RIGHT <- TB_X_LINE + LINE_W + .14 + B_N_W
FM_PANEL_LEFT <- TB_N_RIGHT + .10
FM_LABEL_PT <- (FM_PANEL_LEFT - .06) * 72.27
f1 <- colplot(XL1, c(0, 25, 50, 75, 100), c("0", "25", "50", "75", "100")) + theme(plot.margin = margin(2, 9, 0, FM_LABEL_PT))
f2 <- colplot(XL2, seq(ceiling(XL2[1] / 2) * 2, floor(XL2[2] / 2) * 2, 2), function(x) ifelse(x == 0, "0", sprintf("%+g", x)), xline = 0,
              minor = seq(ceiling(XL2[1]), floor(XL2[2]), 1)) + theme(plot.margin = margin(2, 4, 0, 6), panel.grid.minor.x = element_line(color = "#E4E4E4", linewidth = .2))
NP <- vapply(list(f1, f2), function(p) nonpanel_w(ggplotGrob(p)), 0)
PW <- (RIGHT - .06 - sum(NP)) * c(.6, .4)
C_W <- NP + PW; C_X <- .06 + c(0, C_W[1])
C_PANEL_TOP <- LT[["B"]] + .30 + 7.5 / 72; C_SUB_BASE <- C_PANEL_TOP - .069   # two header lines (Samples / start -> end)
T2G <- align_h(list(f1, f2)); pr2 <- T2G[[1]]$layout$t[T2G[[1]]$layout$name == "panel"]
C_PLOT_Y <- C_PANEL_TOP - sum(convertHeight(T2G[[1]]$heights[seq_len(pr2 - 1L)], "in", valueOnly = TRUE))
C_H <- nonpanel_h(T2G[[1]]) + (YR[1] - YR[2]) * RU
stopifnot(abs(panels_in(T2G[[1]], C_X[1], C_PLOT_Y, C_W[1], C_H)[[1]]$left_in - FM_PANEL_LEFT) < .002)
C_AX_H <- .17
C_KEY_TOP <- C_PLOT_Y + C_H + C_AX_H + .07
K_START_X <- .10 + CIRCLE_R
ipu_C1 <- PW[1] / diff(XL1)
LT[["C"]] <- C_KEY_TOP + .135 + GAP_ROWS
# C: short right-aligned labels; the panels start right of the widest label or section header
GG_PANEL_LEFT <- .06 + max(SEC_W, max(vapply(ROWS$short, str_in, 0, size = 7)) + .05) + .10
GG_LABEL_PT <- (GG_PANEL_LEFT - .06) * 72.27
A_BOX <- c(x = .06, w = RIGHT - .06)
pA0 <- gg_plot(FALSE); gA0 <- ggplotGrob(pA0); hA <- convertHeight(gA0$heights, "in", valueOnly = TRUE); LA <- gA0$layout
prA <- unique(LA$t[grepl("^panel", LA$name)])
A_PLOT_Y <- LT[["C"]] + .257 - sum(hA[seq_len(prA - 1L)])
A_H <- nonpanel_h(gA0) + (YG$yr[1] - YG$yr[2]) * RU
stA <- unique(LA$t[grepl("^strip-t", LA$name)])              # the facet headers' second line: "Grafts (in mice)" sits on it
A_STRIP_BASE <- A_PLOT_Y + sum(hA[seq_len(stA - 1L)]) + hA[stA] / 2 - .5 / 72.27 + (asc_in("G1: EC genes", 7.5, "bold") + 7.5 * .9 * 1.2 / 72) / 2
stopifnot(abs(panels_in(gA0, A_BOX[["x"]], A_PLOT_Y, A_BOX[["w"]], A_H)[[1]]$left_in - GG_PANEL_LEFT) < .002)
nnA <- unitType(gA0$widths) == "null"; nvA <- as.numeric(gA0$widths[nnA]); bA0 <- ggplot_build(pA0); LBA <- bA0$layout$layout
xrA <- vapply(GLAB, function(lab) bA0$layout$panel_params[[LBA$PANEL[LBA$panel == lab]]]$x.range, c(0, 0))
ipu_A <- (A_BOX[["w"]] - sum(convertWidth(gA0$widths[!nnA], "in", valueOnly = TRUE))) * nvA / sum(nvA) / (xrA[2, ] - xrA[1, ])
A_KEY_TOP <- A_PLOT_Y + A_H + .10 + c(0, 0)   # one key line
LT[["D"]] <- A_KEY_TOP[2] + .135 + GAP_ROWS
# D, E: one row; plots start .067 in below the letters
F5_DE_H <- 1.48 + 176 / 600                                  # reference plot height that fixes the marks under D's axis
F5_NP <- nonpanel_h(ggplotGrob(d_plot(1, 1) + labs(x = quote("Stage 6 culture (" * bold("Augsornworawat 2023") * ")")) +
  theme(plot.margin = margin(2, 4, 2, 2), axis.title.x = element_text(size = 7, margin = margin(t = 10.5, b = 2.5)))))
F5_PANEL_H <- F5_DE_H - F5_NP
DE_Y <- LT[["D"]] + .067
E_X <- RIGHT - E_W; D_X <- .06; D_W <- E_X - .12 - D_X
TDE0 <- align_h(list(d_plot(1, 1), pE)); D_PANEL_H <- DE_H - nonpanel_h(TDE0[[1]])
pD <- d_plot(D_PANEL_H, F5_PANEL_H)
DE_KEY_TOP <- DE_Y + DE_H + .02 + c(0, .135)
# D's study names under the axis must not touch: measure their extents in inches on the drawn facet
gD <- ggplotGrob(pD); pD_box <- panels_in(gD, D_X, DE_Y, D_W, DE_H)[[1]]; ipu_D <- pD_box$width_in / (diff(D_XLIM) * 1.1)
sw <- vapply(D_STUDY$lab, str_in, 0, size = 7, face = "bold")
D_STUDY_GAP <- (D_STUDY$x[2] - D_STUDY$x[1]) * ipu_D - sum(sw) / 2
invisible(dev.off()); unlink(MEASURE_PNG)
PAGE_END <- DE_KEY_TOP[2] + .135
stopifnot(PAGE_END < H - .05, D_STUDY_GAP > .05)

# ---- plots with their marks --------------------------------------------------------------------------------------------
PBs$s <- sign(PBs$pct - PBs$src_mean); PBs$hbase <- PBs$pct - PBs$s * STUB_IN / B_IPU
PBs$shaft_end <- PBs$pct - PBs$s * HEAD_AX0 / B_IPU; PBs$bare <- abs(PBs$pct - PBs$src_mean) <= HEAD_AX0 / B_IPU
pB <- prop_plot(TRUE)
PC$thr <- HEAD_MIN_IN / ipu_C1; PC$stub <- STUB_IN / ipu_C1
pC1 <- f1 + lapply(c("SC_SIC", "SC_BETA"), function(k) arrows_ef(PC[PC$population == k, ], POP_COL[[k]]))
pC2 <- f2 + geom_segment(data = CB, aes(x = 0, xend = four_marker, y = y, yend = y, color = I(col)), linewidth = .9)
k <- match(GLAB[PE$group], GLAB); PE$thr <- (HEAD_MIN_IN / ipu_A)[k]; PE$stub <- (STUB_IN / ipu_A)[k]
pA <- gg_plot(TRUE)

# ---- page parts --------------------------------------------------------------------------------------------------------
text_arrow <- function(s, x, y, just = "left") {             # text on the baseline y with each arrow drawn bold
  parts <- strsplit(s, "→", fixed = TRUE)[[1]]; if (endsWith(s, "→")) parts <- c(parts, "")
  runs <- c(rbind(parts, c(rep("→", length(parts) - 1L), "")))[seq_len(2L * length(parts) - 1L)]
  faces <- rep(c("plain", "bold"), length.out = length(runs)); w <- mapply(str_in, runs, 7, faces)
  if (just == "right") x <- x - sum(w)
  for (k in seq_along(runs)) { x0 <- x + sum(w[seq_len(k - 1L)])
    if (runs[k] == "→") grid.lines(x = unit(x0 + w[k] * c(.10, .90), "in"), y = unit(rep(H - y + .034, 2), "in"),
                                     arrow = arrow(angle = 25, length = unit(.032, "in"), type = "closed"),
                                     gp = gpar(col = "black", fill = "black", lwd = 1.1, linejoin = "mitre"))
    else if (nzchar(runs[k])) text_bl(runs[k], x0, y) } }
yv_of <- function(b1, Y = YA) function(v) b1$top_in + (v - Y$yr[2]) / (Y$yr[1] - Y$yr[2]) * b1$height_in
band_full <- function(b1, x0 = .06, x1 = RIGHT, Y = YA) {   # row shading: one unbroken band per shaded row, labels through data
  yv <- yv_of(b1, Y)
  for (r in BAND_ROWS) grid.rect(x = unit(x0, "in"), y = unit(H - yv(r - .5), "in"), width = unit(x1 - x0, "in"), height = unit(yv(r + .5) - yv(r - .5), "in"),
                                 just = c("left", "top"), gp = gpar(fill = "black", alpha = BAND_ALPHA, col = NA)) }
row_base <- function(yc) yc + asc_in("H", 7) / 2             # baseline that centres a 7 pt capital on a row
section_heads <- function(b1, Y = YA, graft_base = NULL) {   # section headers: own lines, 7.5 pt bold, flush left at .06 in
  yv <- yv_of(b1, Y)
  base <- c(graft = if (is.null(graft_base)) yv(Y_SEC[["graft"]]) + asc_in("H", SEC_PT, "bold") / 2 else graft_base,
            cult = yv(Y_SEC[["cult"]]) + asc_in("H", SEC_PT, "bold") / 2)
  for (k in names(base)) text_bl(SECTION[[k]], .06, base[[k]], size = SEC_PT, face = "bold") }
draw_prop <- function() {                                    # A: proportions
  B_LABEL_SPAN <<- list()
  g <- ggplotGrob(pB); draw_at(g, B_PLOT_X, B_Y0, B_PLOT_W, B_H)
  P <- panels_in(g, B_PLOT_X, B_Y0, B_PLOT_W, B_H); col1 <- vapply(P, `[[`, 0, "l") == min(vapply(P, `[[`, 0, "l"))
  P <- P[col1][order(vapply(P[col1], `[[`, 0, "t"))]
  x_r <- B_X + B_STRIP_W
  for (r in seq_along(B_STRIP)) { bx <- P[[r]]; L <- B_STRIP[[r]]; n <- length(L)
    top_ink <- bx$top_in + bx$height_in / 2 - (B_ASC + (n - 1) * B_PITCH + B_DSC) / 2
    B_LABEL_SPAN[[r]] <<- c(top_ink, top_ink + B_ASC + (n - 1) * B_PITCH + B_DSC)
    for (k in seq_len(n)) { runs <- L[[k]]; base <- top_ink + B_ASC + (k - 1) * B_PITCH; xr <- x_r
      for (j in rev(seq_along(runs))) { f <- if (names(runs)[j] == "b") "bold" else "plain"; text_bl(runs[[j]], xr, base, face = f, just = "right"); xr <- xr - str_in(runs[[j]], 7, f) } } }
  bld <- ggplot_build(pB); LB <- bld$layout$layout
  for (i in seq_len(nrow(SRC))) { r <- match(SRC$strip[i], levels(PBs$strip)); pan <- LB$PANEL[LB$ROW == r & LB$COL == 1]
    lv <- levels(droplevels(PBs$glab[PBs$strip == SRC$strip[i]])); yr <- bld$layout$panel_params[[pan]]$y.range
    stopifnot(isTRUE(all.equal(yr, c(.4, length(lv) + .6))))
    j <- match(unique(as.character(PBs$glab[PBs$strip == SRC$strip[i] & PBs$context == SRC$graft[i]])), lv)
    yin <- P[[r]]$top_in + (yr[2] - mean(j)) / diff(yr) * P[[r]]$height_in
    text_arrow(NCOUNT[[SRC$id[i]]], B_N_RIGHT, yin + (asc_in("1", 7) - dsc_in("1", 7)) / 2, just = "right") }
  la <- unique(g$layout$l[grepl("^axis-l", g$layout$name)]); st <- unique(g$layout$t[grepl("^strip-t", g$layout$name)])
  x_right <- B_PLOT_X + sum(convertWidth(g$widths[seq_len(la - 1L)], "in", valueOnly = TRUE)) + max(vapply(levels(PBs$glab), str_in, 0, size = 7))
  base <- P[[1]]$top_in - (convertHeight(g$heights[st], "in", valueOnly = TRUE) / 2 + .5 / 72.27 - asc_in("SC-SIC", 7.5, "bold") / 2)
  text_bl(B_TEXT[["header"]], x_right, base, just = "right")
  text_bl(B_TEXT[["n1"]], B_N_RIGHT, base - 2 * 7.5 / 72, just = "right"); text_arrow(B_TEXT[["n2"]], B_N_RIGHT, base - 7.5 / 72, just = "right")
  PP <- panels_in(g, B_PLOT_X, B_Y0, B_PLOT_W, B_H); kc <- (min(vapply(PP, `[[`, 0, "left_in")) + max(vapply(PP, function(p) p$left_in + p$width_in, 0))) / 2
  gl <- c(range = .20, graft = .18, low = .18); kline <- c(range = 1, graft = 1, low = 2); K <- NULL
  if (!any(PBs$low)) { gl <- gl[1:2]; kline <- kline[1:2] }
  for (ln in unique(kline)) { it <- names(gl)[kline == ln]; x <- kc - (sum(gl[it] + .06 + vapply(B_TEXT[it], str_in, 0)) + .20 * (length(it) - 1)) / 2
    for (k in it) { tl <- x + gl[[k]] + .06; K <- rbind(K, data.frame(item = k, gl = x, gr = x + gl[[k]], tl = tl, top = B_KEY_TOP + (ln - 1) * B_KEY_PITCH)); x <- tl + str_in(B_TEXT[[k]]) + .20 } }
  yc <- K$top + .045; yy <- function(i) unit(rep(H - yc[i], 2), "in")
  grid.rect(x = unit(K$gl[1], "in"), y = unit(H - yc[1], "in"), width = unit(K$gr[1] - K$gl[1], "in"), height = unit(2.6 * ggplot2::.pt / 96, "in"), just = c("left", "center"), gp = gpar(fill = "grey75", col = NA))
  grid.points(x = unit((K$gl[1] + K$gr[1]) / 2, "in"), y = unit(H - yc[1], "in"), pch = 21, size = unit(CIRCLE_PT, "pt"), gp = gpar(col = "#555555", fill = "white", lwd = CIRCLE_LWD))
  grid.lines(x = unit(c(K$gl[2], K$gr[2]), "in"), y = yy(2), arrow = ARW, gp = gpar(col = "#555555", fill = "#555555", lwd = 1.3))
  if (nrow(K) == 3) grid.lines(x = unit(c(K$gl[3], K$gr[3] - HEAD_AX0), "in"), y = yy(3), gp = gpar(col = "#555555", lwd = 1.3, lineend = "butt"))
  if (nrow(K) == 3) grid.lines(x = unit(K$gr[3] - c(STUB_IN, 0), "in"), y = yy(3), arrow = ARW, gp = gpar(col = "#555555", fill = "white", lwd = .30 * ggplot2::.pt))
  for (i in seq_len(nrow(K))) text_at(B_TEXT[[K$item[i]]], K$tl[i], K$top[i] - .005)
  gB <- vapply(seq_len(length(B_LABEL_SPAN) - 1L), function(r) B_LABEL_SPAN[[r + 1L]][1] - B_LABEL_SPAN[[r]][2], 0)
  if (min(gB) < .02) stop(sprintf("A: strip labels %.3f in apart (< .02 in)", min(gB))) }
draw_fm <- function() {                                      # B: four markers with its label table
  t2 <- align_h(lapply(list(pC1, pC2), nobg))
  bx <- lapply(1:2, function(i) panels_in(t2[[i]], C_X[i], C_PLOT_Y, C_W[i], C_H)[[1]]); names(bx) <- c("C1", "C2")
  band_full(bx$C1)
  for (i in 1:2) draw_at(t2[[i]], C_X[i], C_PLOT_Y, C_W[i], C_H)
  yv <- yv_of(bx$C1)
  first <- c(TRUE, ROWS$study[-1] != ROWS$study[-nrow(ROWS)]) | ROWS$id == CULT[1]   # a study is named on its first row in each section
  for (i in seq_len(nrow(ROWS))) { b <- row_base(yv(ROWS$y[i]))
    if (first[i]) text_bl(ROWS$study[i], TB_X_STUDY, b, face = "bold")
    text_bl(ROWS$short[i], TB_X_LINE, b); text_arrow(NCOUNT[[ROWS$id[i]]], TB_N_RIGHT, b, just = "right") }
  section_heads(bx$C1)
  # header lines: table column heads and the two plot heads on the lower line; "Samples" above "start -> end"
  text_bl(TAB_HEAD[["study"]], TB_X_STUDY, C_SUB_BASE); text_bl(TAB_HEAD[["line"]], TB_X_LINE, C_SUB_BASE)
  text_bl(B_TEXT[["n1"]], TB_N_RIGHT, C_SUB_BASE - 7.5 / 72, just = "right"); text_arrow(B_TEXT[["n2"]], TB_N_RIGHT, C_SUB_BASE, just = "right")
  y_ax <- C_PLOT_Y + C_H + convertHeight(calc_element("axis.title.x.bottom", th())$margin[1], "in", valueOnly = TRUE)
  for (p in names(bx)) { b <- bx[[p]]
    grid.text(C_AXIS[[p]], x = unit(b$left_in + b$width_in / 2, "in"), y = unit(H - y_ax, "in"), just = c("centre", "top"), gp = gpar(fontfamily = "Arial", fontsize = 7, col = "black"))
    text_bl(C_HEAD[[p]], b$left_in + b$width_in / 2, C_SUB_BASE, size = 7, just = "centre") }
  yk <- C_KEY_TOP; x <- start_graft_key(K_START_X, yk) + .30
  x <- population_key(x, yk, c("SC_SIC", "SC_BETA")) + .40
  text_at(MARKERS_NOTE[["lead"]], x, yk - .005); text_at(MARKERS_NOTE[["genes"]], x + str_in(MARKERS_NOTE[["lead"]]), yk - .005, face = "italic") }
draw_gg <- function() {                                      # C: gene groups with short right-aligned labels
  g <- ggplotGrob(nobg(pA)); bx <- panels_in(g, A_BOX[["x"]], A_PLOT_Y, A_BOX[["w"]], A_H)
  band_full(bx[[1]], Y = YG); draw_at(g, A_BOX[["x"]], A_PLOT_Y, A_BOX[["w"]], A_H)
  yv <- yv_of(bx[[1]], YG)
  for (i in seq_len(nrow(ROWS))) text_bl(ROWS$short[i], bx[[1]]$left_in - .05, row_base(yv(ROWS$y[i])), just = "right")
  section_heads(bx[[1]], YG, graft_base = A_STRIP_BASE)
  # one complete key line: start, end, the three cell types, the reference line, the pale bars; items .18 in apart
  y2 <- A_KEY_TOP[1]; G <- .18
  x <- start_graft_key(K_START_X, y2) + G; x <- population_key(x, y2) + G
  key_line(x, y2, "#555555", lwd = .4 * ggplot2::.pt, lty = "dashed", butt = TRUE); text_at(A_KEY[["ref"]], x + .26, y2 - .005)
  xb <- x + .26 + str_in(A_KEY[["ref"]]) + G
  if (xb + .26 + str_in(A_KEY[["bar"]]) > RIGHT) stop("C key passes the right edge")
  for (k in names(POP_COL)) grid.lines(unit(c(xb, xb + .2), "in"), unit(rep(H - y2 - .045 - OFF3[[k]] * .13, 2), "in"),
                                       gp = gpar(col = POP_COL[[k]], alpha = .3, lwd = 1 * ggplot2::.pt, lineend = "butt"))
  text_at(A_KEY[["bar"]], xb + .26, y2 - .005) }
DE_GROBS <- function() align_h(lapply(list(pD, pE), function(p) p + theme(plot.background = element_blank())))
draw_D <- function() draw_at(DE_GROBS()[[1]], D_X, DE_Y, D_W, DE_H)
draw_E <- function() draw_at(DE_GROBS()[[2]], E_X, DE_Y, E_W, DE_H)
draw_DE_key <- function() de_key(.10, RIGHT - .10, DE_KEY_TOP[1], DE_KEY_TOP[2])
LETTER_X <- c(A = .10, B = .06, C = .06, D = .10, E = E_X - .03)
draw_letter <- function(L) text_at(L, LETTER_X[[L]], LT[[if (L == "E") "D" else L]], 11, "bold")
