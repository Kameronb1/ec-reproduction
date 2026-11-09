# Figure S5C: share of the EC cells (primary study label) of each donor, sample, day (Singh 2024) or HIO sample group
# (McCauley 2023) classified as EC by the four-marker, TPH1-only and 20-gene classifiers, one subplot per study in
# three columns (Elmentaite 2021 | Hickey 2023 | Egozi 2021, Singh 2024, McCauley 2023), rows banded light grey, with
# each unit's EC cells on the right. Key under the shortest column; under the Hickey 2023 column, the share of each
# study's other terminal EEC cells that each classifier classifies as EC.

tissue <- read_table("calls_by_donor.tsv")                     # test studies: EC cells, and how many each classifier classifies as EC
bh <- read_table("hio_calls_by_organoid.tsv")                  # McCauley 2023's two HIO sample groups
SINGH <- SCORED[4]; EGOZI <- SCORED[3]
du <- rbind(tissue, bh[names(tissue)])
t1 <- read_table("tph1_calls_by_unit.tsv")                     # the TPH1-only classifier's counts on the same units
m <- match(paste(du$study, du$unit), paste(t1$study, t1$unit)); stopifnot(!anyNA(m), all(t1$n_rule[m] == du$n_rule))
du$called_tph1 <- t1$called_tph1[m]
oe <- read_table("other_eec_calls.tsv")                        # other terminal EEC classified as EC, per study and classifier
oe$st <- factor(study(oe$study), rev(study(SCORED)))
ob <- do.call(rbind, lapply(1:3, function(j) data.frame(st = oe$st, k = c("Four-marker classifier", "TPH1-only classifier", "20-gene classifier")[j],
  called = 100 * oe[[c("called_rule", "called_tph1", "called_cls")[j]]] / oe$n, j = j)))
ob$y <- as.integer(ob$st) + c(.27, 0, -.27)[ob$j]; OBH <- .22
pct <- function(k, n) 100 * k / n
du$rule <- pct(du$called_rule, du$n_rule); du$tph1 <- pct(du$called_tph1, du$n_rule); du$cls_pct <- pct(du$called_cls, du$n_rule)   # all three on the same EC cells
day <- suppressWarnings(as.numeric(sub("D", "", du$unit)))                           # Singh 2024 days; NA elsewhere
du <- du[order(match(du$study, SCORED), ifelse(du$study == SINGH, day, -du$n_rule)), ]   # days in order, else largest first
du$ylab <- du$unit
du$ylab[du$study == EGOZI] <- gsub("_", " ", du$unit[du$study == EGOZI], fixed = TRUE)
du$ylab[du$study == SINGH] <- sub("^D", "Day ", du$unit[du$study == SINGH])

# row pitch (in), dot size, dot offset from the row center (rows), x range, x of the EC-cell counts, gaps (in)
PT16 <- .135; DOT16 <- 1.75; OFF16 <- .30; X16 <- c(-3, 131); NX16 <- 106; GAP16 <- .07; COLGAP16 <- .22
KEY16 <- c("Four-marker classifier" = RULE, "TPH1-only classifier" = TPH1_COL, "20-gene classifier" = PANEL)
TITLE16 <- setNames(study(SCORED), SCORED); TITLE16[[2]] <- paste(TITLE16[[2]], "(nuclei)")
TYPE16 <- setNames(study_type(SCORED), SCORED)
COL16 <- list(c1 = SCORED[1], c2 = SCORED[2], c3 = SCORED[3:5])
BANDL16 <- c(c1 = .30, c2 = .30, c3 = .64)     # row-label width the bands reach over, per column (in)
XT16 <- "Called EC (% of EC cells)"            # x-axis title under each column
make16 <- function(s, pw, xt = NULL) {
  d <- du[du$study == s, ]; d$i <- rev(seq_len(nrow(d)))                              # first row on top
  pts <- rbind(data.frame(i = d$i + OFF16, x = d$rule, k = names(KEY16)[1]), data.frame(i = d$i, x = d$tph1, k = names(KEY16)[2]),
               data.frame(i = d$i - OFF16, x = d$cls_pct, k = names(KEY16)[3]))
  cn <- names(COL16)[vapply(COL16, function(k) s %in% k, TRUE)]
  top <- s %in% vapply(COL16, `[`, "", 1); bottom <- s %in% vapply(COL16, function(v) v[length(v)], "")
  upi <- diff(X16) / pw; bl <- X16[1] - BANDL16[[cn]] * upi
  bd <- d[(nrow(d) - d$i) %% 2 == 0, ]
  ggplot() + geom_rect(data = bd, aes(xmin = bl, xmax = X16[2], ymin = i - .5, ymax = i + .5), fill = "#F0F0F0", color = NA) +
    geom_vline(xintercept = c(0, 50, 100), color = GRID_COL, linewidth = .2) +
    geom_point(data = pts, aes(x, i, fill = k), size = DOT16, shape = 21, color = "black", stroke = .2) +
    geom_text(data = d, aes(x = NX16, y = i, label = n_rule), hjust = 0, size = tx, family = FONT) +
    (if (top) annotate("text", x = NX16, y = nrow(d) + .5 + .04 / PT16, label = "EC cells", hjust = 0, vjust = 0, size = tx, family = FONT, fontface = "bold") else NULL) +
    scale_fill_manual(values = KEY16, guide = "none") +
    scale_y_continuous(breaks = d$i, labels = d$ylab, expand = c(0, 0)) +
    scale_x_continuous(breaks = c(0, 50, 100), labels = function(v) paste0(v, "%"), expand = c(0, 0)) +
    guides(x = guide_axis(cap = "both")) + coord_cartesian(xlim = X16, ylim = c(.5, nrow(d) + .5), expand = FALSE, clip = "off") +
    labs(x = xt, y = NULL, title = TITLE16[[s]], subtitle = TYPE16[[s]]) + theme_kit("x") +
    theme(panel.widths = unit(pw, "in"), panel.heights = unit(nrow(d) * PT16, "in"),
          plot.title = element_text(size = 8, face = "bold", margin = margin(0, 0, 1, 0)),
          plot.subtitle = element_text(size = PT, margin = margin(0, 0, 2, 0)), plot.margin = margin(2, 2, 2, 2)) +
    (if (bottom) NULL else theme(axis.text.x = element_blank()))
}
grobs16 <- function(pw) {                       # the bottom subplot of each column carries the x-axis title
  bottoms <- vapply(COL16, function(v) v[length(v)], "")
  g <- setNames(lapply(SCORED, function(s) ggplotGrob(make16(s, pw, if (s %in% bottoms) XT16))), SCORED)
  for (k in COL16) { cw <- do.call(unit.pmax, lapply(g[k], function(x) x$widths)); for (s in k) g[[s]]$widths <- cw }   # align a column
  g
}
gin <- function(g, wh = "w") if (wh == "w") sum(convertWidth(g$widths, "in", valueOnly = TRUE)) else sum(convertHeight(g$heights, "in", valueOnly = TRUE))
# The panel is 7.4 in wide and starts 0.1 in left of the other panels; the columns fill it to 7.36 in, the "EC cells"
# header of the last column included. Panel width: the same in every column, measured from subplots 1.5 in wide.
W_C <- 7.4; RIGHT_C <- 7.36
g0 <- grobs16(1.5); fixedj <- vapply(COL16, function(k) gin(g0[[k[1]]]) - 1.5, 0)
HDR <- str_in("EC cells", face = "bold"); RN <- (X16[2] - NX16) / diff(X16)          # the header runs past the panel by HDR - RN * pw
PW16 <- floor(100 * (RIGHT_C - sum(fixedj) - 2 * COLGAP16 - HDR + 2 / 72) / (3 - RN)) / 100
gC <- grobs16(PW16); CWj <- vapply(COL16, function(k) gin(gC[[k[1]]]), 0)
XJ <- c(0, cumsum(CWj + COLGAP16))[1:3]
colh <- vapply(COL16, function(k) sum(vapply(gC[k], gin, 0, wh = "h")) + GAP16 * (length(k) - 1), 0)
SUBMAX <- max(vapply(SCORED, function(s) str_in(TYPE16[[s]]), 0))
stopifnot((NX16 - X16[1]) / diff(X16) * PW16 - SUBMAX >= .08,                       # subtitles clear the "EC cells" header
          XJ[3] + CWj[[3]] + HDR - RN * PW16 - 2 / 72 <= RIGHT_C + 1e-9)
H_C <- max(colh) + .02
H_CP <- ceiling((.47 + H_C) * 100) / 100                                              # panel height (in)
draw_c <- function() {
  vh <- convertHeight(unit(1, "npc"), "in", valueOnly = TRUE)
  for (j in seq_along(COL16)) { y <- 0; x <- XJ[j]; CW <- CWj[[j]]
    for (s in COL16[[j]]) { g <- gC[[s]]; h <- gin(g, "h")
      pushViewport(viewport(x = unit(x, "in"), y = unit(vh - y, "in"), width = unit(CW, "in"), height = unit(h, "in"), just = c("left", "top")))
      grid.draw(g); popViewport(); y <- y + h + GAP16 } }
  ky <- colh[["c1"]] + .22; x <- XJ[1] + .35                                             # key under the shortest column (Elmentaite 2021)
  for (i in seq_along(KEY16)) {
    grid.points(unit(x + .05, "in"), unit(vh - ky - (i - 1) * .16, "in"), pch = 21, size = unit(1, "char"), gp = gpar(col = "black", fill = KEY16[[i]], lwd = .2 * .pt, fontsize = DOT16 * .pt))
    grid.text(names(KEY16)[i], unit(x + .13, "in"), unit(vh - ky - (i - 1) * .16, "in"), just = c("left", "centre"), gp = gpar(fontfamily = FONT, fontsize = PT, col = "black")) }
  # under the Hickey 2023 column: other terminal EEC cells classified as EC by each classifier
  gx <- XJ[2]; gy <- colh[["c2"]] + .10; gw <- CWj[[2]] + .10; gh <- vh - gy - .08
  pushViewport(viewport(x = unit(gx, "in"), y = unit(vh - gy, "in"), width = unit(gw, "in"), height = unit(gh, "in"), just = c("left", "top")))
  grid.draw(ggplotGrob(oe_plot)); popViewport()
}

oe_plot <- ggplot(ob) +
  geom_rect(data = data.frame(i = seq(length(levels(ob$st)), 1, by = -2)), aes(xmin = 0, xmax = 8, ymin = i - .5, ymax = i + .5),
            fill = "#F0F0F0", inherit.aes = FALSE) +                 # alternate rows banded, top row first, as in the subplots above
  geom_rect(aes(xmin = 0, xmax = called, ymin = y - OBH / 2, ymax = y + OBH / 2, fill = k)) +
  scale_fill_manual(values = KEY16, guide = "none") +
  scale_y_continuous(breaks = seq_along(levels(ob$st)), labels = levels(ob$st), expand = c(0, 0), limits = c(.5, length(levels(ob$st)) + .5)) +
  scale_x_continuous(limits = c(0, 8), breaks = c(0, 2, 4, 6, 8), labels = function(v) paste0(v, "%"), expand = c(0, 0)) +
  labs(x = "Other terminal EEC called EC (% of cells)", y = NULL) + coord_cartesian(clip = "off") + theme_kit("x") + theme(plot.margin = margin(2, 6, 2, 2))

# The row bands of the first two columns reach past the panel's left edge, which cuts them off (the subplots draw
# unclipped, so the strip left of the panel is painted white again).
draw_panel <- function() {
  place(draw_c, 60, 4134, 4440, 1836, W_C, H_CP - (.30 + .17))
  grid.rect(x = unit(0, "in"), y = unit((H * DPI - 4134) / DPI, "in"), width = unit(60 / DPI, "in"), height = unit(1836 / DPI, "in"),
            just = c("left", "top"), gp = gpar(fill = "white", col = NA))
}
