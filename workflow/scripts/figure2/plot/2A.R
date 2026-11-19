# Figure 2A: how well each population's EC-versus-rest signature (per-gene fold changes against the other endocrine
# cells of its study) agrees with intestinal EC's. For each population and study: its Spearman rho with each intestinal
# study's EC signature (open dots, joined by a line) and their median (filled dot), the dots of a row spread vertically
# so every one shows; the grey band spans the agreement among intestinal EC studies. Key below the panel.

C_HIO <- "McCauley 2023 both organoids"                 # HIO graft EC: the two organoids pooled
dots <- read_table("signature_dots.tsv"); hio <- read_table("signature_hio_dots.tsv")
hio <- hio[hio$study == C_HIO, ]; hio$population <- "HIO_EC"
cd <- rbind(dots, hio)[, c("population", "study", "spearman_median")]
pw <- read_table("signature_pairwise.tsv"); pwh <- read_table("signature_hio_pairwise.tsv")
pwh <- pwh[pwh$study == C_HIO, ]; pwh$population <- "HIO_EC"
cp <- rbind(pw, pwh)[, c("population", "study", "gut_study", "spearman")]
# each filled dot is the median of its comparisons with the intestinal studies
m <- merge(aggregate(spearman ~ population + study, data = cp, FUN = median), cd, by = c("population", "study"))
stopifnot(nrow(m) == nrow(cd), max(abs(m$spearman - m$spearman_median)) < 1e-12)
# SC-alpha's six study dots are in the tables but not drawn
stopifnot(sum(cd$population == "SC_ALPHA") == 6)
cd <- cd[cd$population != "SC_ALPHA", ]
C_ROWS <- c(GUT_EC = "Intestinal EC vs\nother EEC", SC_SIC = "SC-SIC vs other\nSC endocrine", SC_BETA = "SC-beta vs other\nSC endocrine",
            HIO_EC = "HIO graft EC\n(EC vs other EEC)")   # top to bottom
C_COLS <- c(GUT_EC = "#43A5B6", SC_SIC = "#AA3377", SC_BETA = "#4477AA", HIO_EC = "#6943AC")
stopifnot(setequal(cd$population, names(C_ROWS)))
# each dot's height = its row number (bottom row 1) plus an offset that spreads the row's dots apart
# (positions_a_dots.tsv: the offset of each population-study dot within its row)
off <- read_placement("positions_a_dots.tsv")
cd$y <- match(cd$population, rev(names(C_ROWS))) + off$offset[match(paste(cd$population, cd$study), paste(off$population, off$study))]
stopifnot(!anyNA(cd$y))
band <- range(cd$spearman_median[cd$population == "GUT_EC"])
cp <- cp[cp$population != "SC_ALPHA", ]
cp <- merge(cp, cd[, c("population", "study", "y")], by = c("population", "study"))
stopifnot(all(paste(cd$population, cd$study) %in% paste(cp$population, cp$study)))
cr <- do.call(rbind, lapply(split(cp, paste(cp$population, cp$study)), function(d)
  data.frame(population = d$population[1], y = d$y[1], lo = min(d$spearman), hi = max(d$spearman))))
C_TOP <- length(C_ROWS) + 0.6                           # y axis top: a little room above the top row
# the key, below the panel, is drawn as shapes: open dot = each intestinal study, filled dot = median, grey square = the
# range among intestinal studies
KEY <- data.frame(x = 0, y = 1, what = factor(c("each intestinal study", "median", "intestinal range"), levels = c("each intestinal study", "median", "intestinal range")))

plot_A <- ggplot(cd, aes(spearman_median, y)) +
  annotate("rect", xmin = band[1], xmax = band[2], ymin = -Inf, ymax = Inf, fill = "#E6E6E6", colour = NA) +
  geom_vline(xintercept = 0, colour = "grey70", linewidth = 0.3) +
  geom_segment(data = cr, aes(x = lo, xend = hi, y = y, yend = y, colour = population), linewidth = 0.3, inherit.aes = FALSE) +
  geom_point(data = cp, aes(spearman, y, colour = population), shape = 21, fill = "white", size = 0.85, stroke = 0.3, inherit.aes = FALSE) +
  geom_point(aes(fill = population), shape = 21, colour = "black", size = 1.3, stroke = 0.15) +
  geom_point(data = KEY, aes(x, y, shape = what), inherit.aes = FALSE, alpha = 0) +   # invisible; only feeds the key below the panel
  scale_shape_manual(name = NULL, values = c(21, 21, 22), breaks = KEY$what,
                     guide = guide_legend(nrow = 2, byrow = TRUE, override.aes = list(alpha = 1, fill = c("white", "grey35", "#E6E6E6"),
                       colour = c("grey35", "black", "#E6E6E6"), size = c(0.85, 1.3, 2.2), stroke = c(0.3, 0.15, 0)))) +
  scale_fill_manual(values = C_COLS, guide = "none") + scale_colour_manual(values = C_COLS, guide = "none") +
  scale_y_continuous(breaks = seq_along(C_ROWS), labels = rev(unname(C_ROWS)), limits = c(0.6, C_TOP), expand = c(0, 0)) +
  scale_x_continuous(breaks = c(-0.1, 0, 0.1, 0.2, 0.3, 0.4)) +
  labs(x = "Spearman ρ with intestinal EC", y = NULL) +
  theme_fig2() + theme(panel.grid.major.y = element_blank(), axis.text.y = element_text(lineheight = 0.82),
                       legend.position = "bottom", legend.margin = margin(0, 0, 0, 0), legend.box.spacing = unit(2, "pt"),
                       legend.key.size = unit(6, "pt"), legend.spacing.x = unit(1, "pt"), legend.text = element_text(size = theme_fig2()$axis.text$size))

draw_panel <- function() plot_at(plot_A, "A_plot")
