# Supplementary Figure S1C: intestinal EC fold changes followed in each population. Each point a gene judgeable in both;
# x = intestinal EC vs other EEC (median of 3 intestinal studies), y = the population vs its comparison cells; colored =
# intestinal EC signature genes; black = the 18 marker genes, left out of rho and the fitted line (five named per plot).
# Header of each plot: Spearman rho without the marker genes (concordance.tsv) and how many signature genes are also up
# 2-fold (share_two_fold.tsv).

G <- read_table("gene_signatures.tsv.gz"); names(G)[1] <- "gene"
for (b in grep("judgeable|up_in_half|excluded18", names(G), value = TRUE)) G[[b]] <- G[[b]] == "True"
CON <- read_table("concordance.tsv"); CON <- CON[CON$genes == "without them", ]
TWO <- read_table("share_two_fold.tsv"); TWO <- TWO[TWO$gene_set == "own", ]
POPC <- c(SC_SIC = "SC-SIC\nvs other SC endocrine", SC_BETA = "SC-beta vs other SC\nendocrine (not SC-SIC)",
          SC_ALPHA = "SC-alpha vs other SC\nendocrine (not SC-SIC)", HIO_GRAFT_EC = "HIO graft EC\nvs other EEC (McCauley 2023)")
LC <- do.call(rbind, lapply(names(POPC), function(p) {
  ok <- G$gut_judgeable & G[[paste0(p, "_judgeable")]] & !is.na(G[[paste0(p, "_lfc")]]) & !is.na(G$gut_ec_lfc)
  d <- data.frame(gene = G$gene, x = G$gut_ec_lfc, y = G[[paste0(p, "_lfc")]], up = G$gut_ec_up_in_half, mk = G$excluded18, pop = p)[ok, ]
  rho <- CON$spearman[CON$population == p]; n2 <- TWO$n_also_up_2fold[TWO$population == p]; nup <- TWO$n_gut_ec_up[TWO$population == p]
  k <- d[!d$mk, ]; fit <- coef(lm(y ~ x, k))                                # the fitted line drawn on the plot
  d$strip <- sprintf("%s\nρ = %s; %d of %d up 2-fold", POPC[[p]], sub("-", "−", sprintf("%.2f", rho), fixed = TRUE), n2, nup)
  d$a <- fit[[1]]; d$b <- fit[[2]]; d }))
LC$strip <- factor(LC$strip, unique(LC$strip)); LC$col <- POPCOL[LC$pop]
LC$kind <- ifelse(LC$mk, "marker", ifelse(LC$up, "up", "other"))
FITC <- unique(LC[c("strip", "pop", "a", "b", "col")])
LABC_GENES <- c("TPH1", "SLC18A1", "LMX1A", "INS", "GCG")                 # a few of the 18 marker genes named on each plot
XL_C <- c(-8.5, 7); YL_C <- c(-8, 10)
LABC <- LC[LC$mk & LC$gene %in% LABC_GENES & LC$x >= XL_C[1] & LC$x <= XL_C[2] & LC$y >= YL_C[1] & LC$y <= YL_C[2], ]
# each gene name is set in clear space beside the cloud and joined to its black point by a short leader. Positions are
# label centres in data units; the leader ends at the label box (box size from the string width at 7 pt and the plot's
# units per inch).
LABPOS <- read.table(header = TRUE, stringsAsFactors = FALSE, text = '
pop          gene     cx    cy
SC_SIC       GCG     -7.1  -0.7
SC_SIC       INS     -1.8  -4.4
SC_SIC       SLC18A1 -3.2   7.6
SC_SIC       LMX1A   -2.9   5.9
SC_SIC       TPH1     4.9   7.6
SC_BETA      GCG     -7.1  -2.4
SC_BETA      INS     -4.6   3.8
SC_BETA      SLC18A1  1.9  -5.0
SC_BETA      LMX1A    2.3   5.9
SC_BETA      TPH1     5.5  -3.6
SC_ALPHA     GCG     -5.6   5.9
SC_ALPHA     INS     -1.9  -4.2
SC_ALPHA     SLC18A1  2.4   7.6
SC_ALPHA     LMX1A    1.9  -5.6
SC_ALPHA     TPH1     5.5  -4.0
HIO_GRAFT_EC INS     -5.0  -0.4
HIO_GRAFT_EC LMX1A   -2.4   9.1
HIO_GRAFT_EC SLC18A1 -3.4   7.3
')
C_UPI <- c(x = 15.5 / 1.57, y = 18 / 1.50)                                # data units per inch of a plot (measured on the page)
local({ page <- dev.cur()                                                  # label widths measured on a 72 dpi device
  png(tempfile(fileext = ".png"), width = 2, height = 2, units = "in", res = 72, type = "cairo")
  LABPOS$hw <<- vapply(LABPOS$gene, function(g) str_in(g, PT, "italic"), 0) / 2 * C_UPI[["x"]] + .08
  dev.off(); dev.set(page) })
LABPOS$hh <- .55 * PT / 72 * C_UPI[["y"]] + .15
LABC <- merge(LABC, LABPOS, by = c("pop", "gene"))
if (nrow(LABC) != nrow(LABPOS)) stop("panel C: a placed label without its point in range")
LABC$ex <- pmin(pmax(LABC$x, LABC$cx - LABC$hw), LABC$cx + LABC$hw); LABC$ey <- pmin(pmax(LABC$y, LABC$cy - LABC$hh), LABC$cy + LABC$hh)
LABC$lead <- sqrt((LABC$ex - LABC$x)^2 + (LABC$ey - LABC$y)^2) > .5     # a leader only where the label is not right at its point
if (any(LABC$cx - LABC$hw < XL_C[1] | LABC$cx + LABC$hw > XL_C[2] | LABC$cy - LABC$hh < YL_C[1] | LABC$cy + LABC$hh > YL_C[2])) stop("panel C: a label runs past the plot edge")
pC <- ggplot(LC, aes(x, y)) +
  geom_hline(yintercept = 0, colour = "grey55", linewidth = .25) + geom_vline(xintercept = 0, colour = "grey55", linewidth = .25) +
  geom_abline(slope = 1, intercept = 0, colour = "grey45", linewidth = .3, linetype = "22") +
  geom_point(data = LC[LC$kind == "other", ], colour = "grey80", size = .2, stroke = 0) +
  geom_point(data = LC[LC$kind == "up", ], aes(fill = I(col)), shape = 21, colour = "black", size = .75, stroke = .08) +
  geom_abline(data = FITC, aes(slope = b, intercept = a, colour = I(col)), linewidth = .5) +
  geom_point(data = LC[LC$mk, ], colour = "black", size = .9, stroke = 0) +
  geom_segment(data = LABC[LABC$lead, ], aes(x = x, y = y, xend = ex, yend = ey), colour = "black", linewidth = .2) +
  geom_text(data = LABC, aes(cx, cy, label = gene), size = MM, family = FONT, fontface = "italic", colour = "black") +
  facet_wrap(~strip, nrow = 1) +
  scale_x_continuous(breaks = c(-8, -4, 0, 4), labels = minus) + scale_y_continuous(breaks = c(-5, 0, 5), labels = minus) +
  coord_cartesian(xlim = XL_C, ylim = YL_C, expand = FALSE) +
  labs(x = "Intestinal EC vs other EEC (log2 fold change, median of 3 intestinal studies)",
       y = "Population vs its comparison\ncells (log2 fold change)") +   # plain log2: a subscript would print below 7 pt
  theme_pg("both") + theme(panel.spacing.x = unit(.12, "in"))

draw_panel <- function() draw_after_72dpi(function() {
  letter("C", LX, ROW2); place(pC, X0, ROW2, XR - X0, C_H)
  key_row(list(item("Other genes", function(x, y) grid.points(unit(x, "in"), unit(H - y, "in"), pch = 16, size = unit(.04, "in"), gp = gpar(col = "grey80")), .12),
               item("Intestinal EC-up genes (2-fold in 2 of 3 intestinal studies)", function(x, y) for (k in 0:3) sym_dot(x - .09 + k * .06, y, POPCOL[c("SC_SIC", "SC_BETA", "SC_ALPHA", "HIO_GRAFT_EC")][k + 1], .05), .30),
               item("18 marker genes, left out of ρ and the fit", function(x, y) grid.points(unit(x, "in"), unit(H - y, "in"), pch = 16, size = unit(.055, "in"), gp = gpar(col = "black")), .12)),
          X0, XR - X0, KEY_C_Y)
  key_row(list(item("Fit without the marker genes", function(x, y) for (k in 0:3) sym_line(x - .105 + k * .07, y, POPCOL[c("SC_SIC", "SC_BETA", "SC_ALPHA", "HIO_GRAFT_EC")][k + 1], w = .07), .30),
               item("Equal fold change", function(x, y) sym_line(x, y, "grey45", lty = "22", lwd = .3, w = .26), .32)), X0, XR - X0, KEY_C_Y + .15)
})
