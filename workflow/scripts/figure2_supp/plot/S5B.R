# Figure S5B: the 20 genes (across, in Fig. 2D's order) and their weights, as % of the classifier's total, in the 20-gene classifier
# of each scored study (down), with each classifier's score cutoff on the right. The rows above the dotted line are each
# trained on the primary study labels of the other studies above it; Singh 2024 and McCauley 2023 share one classifier.
# Key row underneath: the weight color bar and the dotted line's meaning.

ELM <- "SF_ELMENTAITE2021_GUT_ATLAS"; HIC <- "SF_HICKEY2023_ADULT_EEC"; EGO <- "SF_EGOZI2021_DEVELOPMENTAL_EEC"
SINGH <- "SF_INTESTINAL_GSE238274"; HIO <- "SF_GSE214852_HIO_GRAFT"
pred <- read_table("predictions.tsv.gz")
cutoff <- function(s) unique(pred$threshold[pred$study == s & pred$panel_size == 20])
Dt <- read_table("D_effect_heatmap.tsv", "figure2")
ORD <- unique(Dt$gene_symbol[order(Dt$panel_rank)])                                   # Fig. 2D's 20 genes, in its order
wb <- read_table("weights.tsv"); wb <- wb[wb$panel_size == 20 & wb$study %in% c(ELM, HIC, EGO, SINGH), ]
# row labels: runs of c(text, fontface) on two lines
ROWS_B <- list(list(id = ELM, l1 = list(c(study(ELM), "bold")), l2 = list(c(study_type(ELM), "plain"))),
               list(id = HIC, l1 = list(c(study(HIC), "bold")), l2 = list(c(study_type(HIC), "plain"))),
               list(id = EGO, l1 = list(c(study(EGO), "bold")), l2 = list(c(study_type(EGO), "plain"))),
               list(id = SINGH, l1 = list(c(study(SINGH), "bold"), c(paste0(", ", study_type(SINGH)), "plain")),
                    l2 = list(c(study(HIO), "bold"), c(paste0(", ", study_type(HIO)), "plain"))))
WR <- range(100 * wb$weight); PALW <- scales::seq_gradient_pal("white", "#5FA9B6", "Lab")   # intestinal classifier teal, light
fillw <- function(v) PALW(scales::rescale(v, from = WR))
TX0 <- 1.62; TW <- .245; TX1 <- TX0 + 20 * TW; CUTX <- TX1 + .42                       # tiles; cutoff column center (in)
Y0 <- .62; RHB <- .25; SEP <- .06                                                      # first tile row top; row height; gap at dotted line
row_top <- function(i) Y0 + (i - 1) * RHB + (i == 4) * SEP
B_BOT <- row_top(4) + RHB
KEY_Y <- B_BOT + .14; H_B <- KEY_Y + .26                                               # key row; content height (in)
H_BP <- ceiling((.47 + H_B) * 100) / 100                                               # panel height (in)
draw_b <- function() {
  vh <- convertHeight(unit(1, "npc"), "in", valueOnly = TRUE)
  tt <- function(s, x, y, face = "plain", just = c("left", "centre"), rot = 0)
    grid.text(s, x = unit(x, "in"), y = unit(vh - y, "in"), just = just, rot = rot, gp = gpar(fontfamily = FONT, fontsize = PT, fontface = face, col = "black"))
  text_runs <- function(p, x, y) for (s in p) { tt(s[1], x, y, s[2]); x <- x + str_in(s[1], face = s[2]) }
  for (k in seq_along(ORD)) tt(ORD[k], TX0 + (k - .5) * TW, Y0 - .05, "italic", c("left", "centre"), rot = 90)   # gene names
  tt("Score cutoff", CUTX, Y0 - .09, "bold", c("centre", "centre"))
  for (i in 1:4) { r <- ROWS_B[[i]]; yt <- row_top(i); ym <- yt + RHB / 2
    text_runs(r$l1, .04, ym - .062); text_runs(r$l2, .04, ym + .062)
    w <- wb[wb$study == r$id, ]; w <- w[match(ORD, w$gene_symbol), ]
    for (k in seq_along(ORD)) {
      grid.rect(unit(TX0 + (k - 1) * TW, "in"), unit(vh - yt, "in"), unit(TW, "in"), unit(RHB, "in"), just = c("left", "top"),
                gp = gpar(fill = fillw(100 * w$weight[k]), col = "white", lwd = .4 * .pt))
      tt(sprintf("%.1f", 100 * w$weight[k]), TX0 + (k - .5) * TW, ym, just = c("centre", "centre")) }
    tt(sprintf("%.2f", cutoff(r$id)), CUTX, ym, just = c("centre", "centre")) }
  yd <- row_top(4) - SEP / 2                                                            # dotted line between trained and untrained
  grid.lines(unit(c(.04, 7.26), "in"), unit(rep(vh - yd, 2), "in"), gp = gpar(col = "#777777", lty = "dotted", lwd = .3 * .pt))
  # key row: weight color bar (left), dotted-line meaning (right)
  tt("Weight (% of row)", .04, KEY_Y)
  bx0 <- .04 + str_in("Weight (% of row)") + .15; bw <- 1.5; by <- KEY_Y - .045
  v <- seq(WR[1], WR[2], length.out = 200)
  grid.rect(unit(bx0 + (seq_along(v) - 1) / length(v) * bw, "in"), unit(vh - by, "in"), unit(bw / length(v) * 1.5, "in"), unit(.09, "in"),
            just = c("left", "top"), gp = gpar(fill = fillw(v), col = NA))
  grid.rect(unit(bx0, "in"), unit(vh - by, "in"), unit(bw, "in"), unit(.09, "in"), just = c("left", "top"), gp = gpar(fill = NA, col = "#999999", lwd = .25 * .pt))
  for (b in c(4, 6, 8, 10)) { xb <- bx0 + (b - WR[1]) / diff(WR) * bw
    grid.lines(unit(rep(xb, 2), "in"), unit(vh - by - c(.09, .12), "in"), gp = gpar(col = "#999999", lwd = .25 * .pt))
    tt(as.character(b), xb, by + .19, just = c("centre", "centre")) }
  nx <- bx0 + bw + .45
  grid.lines(unit(nx + c(0, .3), "in"), unit(rep(vh - KEY_Y, 2), "in"), gp = gpar(col = "#777777", lty = "dotted", lwd = .3 * .pt))
  tt("Above the line: used to rank Fig. 2D's genes", nx + .38, KEY_Y - .065)
  tt("Below: not used to rank Fig. 2D's genes or to train", nx + .38, KEY_Y + .065)
}

draw_panel <- function() place(draw_b, 120, 2820, 4380, 1284, 7.3, H_BP - (.30 + .17))
