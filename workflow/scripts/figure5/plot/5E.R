# Figure 5E: identity factors. One row per selected factor (two groups: targets higher or lower in SC-SIC) and one dot
# per regulon (one per network), in three plots: SC-islets (SC-SIC minus SC-beta regulon target score, median over the
# SC-islet studies), intestine (EC minus other enteroendocrine, median over the intestinal studies) and adult islets
# (beta minus other endocrine, median over the adult islet studies, not used for the colors; axis broken at GG_AD_CUT). Regulons reproducible in SC-islets and intestine are colored by
# direction, the rest grey. Everything is drawn into the panel's slot; positions below are inches from its top-left.

G_W <- 7.40; G_H <- 3.26                            # slot size
G_TOP <- E_TOP                                      # slot top on the page (helpers.R)
BLOCK_TOP <- .40 - .16 - 72 / 600                   # top of the three-plot block in the slot
BLOCK_H <- G_H - .40 - .04 + .1108333 + 72 / 600    # its height (with room for three-line x titles)
G_STUDY <- c("ELMENTAITE", "HICKEY_2023", "EGOZI_2021")
is_true <- function(x) !is.na(x) & as.character(x) %in% c("TRUE", "True", "true")

# ---- values: selected factors and their regulons ------------------------------------------------------------------------
GGF <- read_table("identity/factors.tsv"); GGF <- GGF[is_true(GGF$selected), ]
GGV <- read_table("identity/program_values.tsv"); GGV$pass <- is_true(GGV$program_pass)

# ---- rows: two groups by direction, factors alphabetical within a group ---------------------------------------------------
GG_HALF <- c("Targets higher\nin SC-SIC", "Targets lower\nin SC-SIC")
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
GG_COLS <- c(pass_ec = unname(POP_COL["SC-SIC"]), pass_beta = unname(POP_COL["SC-beta"]), fail = "grey45")
GG_NL <- data.frame(factor = GGF$factor, half = GGF$half, y = as.numeric(factor(GGF$factor, rev(GGF$factor))), expr = GGF$expr)
# grey band on every second row; the lowest row gets 0.2 row of padding below it
GG_NL$band <- (seq_len(nrow(GG_NL)) %% 2L) == 0L
GG_YMIN <- min(GG_NL$y); gg_pad <- function(y) ifelse(y == GG_YMIN, .2, 0)
GG_BANDS <- transform(GG_NL[GG_NL$band, ], ymin = y - .5 - gg_pad(y), ymax = y + .5)
GG_YB <- rbind(transform(GG_NL, y = y - .5 - gg_pad(y)), transform(GG_NL, y = y + .5))
# ---- text ---------------------------------------------------------------------------------------------------------------
GG_PTITLE <- c(sc = "SC-islets", gut = "Intestine", ad = "Adult islets")
# x titles: two lines each; "regulon target score" and which column sets the colors go to the legend and key
GG_XLAB <- c(sc = "SC-SIC minus SC-beta,\nmedian over SC-islet studies", gut = "EC minus other enteroendocrine,\nmedian over intestinal studies",
             ad = "Beta minus other endocrine,\nmedian over adult islet studies")
GG_AD_CUT <- -.075   # adult islets: axis break here; values below it go to a condensed stub on the left
GG_KEY_HEAD <- "Reproducible in SC-islets and intestine:"
GG_KEYS <- list(list("higher in SC-SIC", GG_COLS[["pass_ec"]]), list("lower in SC-SIC", GG_COLS[["pass_beta"]]), list("not reproducible", GG_COLS[["fail"]]))
GG_PM <- c("(+) activated targets", "(−) repressed targets")   # key under the row-group labels

# ---- the three plots ----------------------------------------------------------------------------------------------------
gg_th <- function() theme_minimal(base_family = "Arial", base_size = 7) + theme(
  text = element_text(color = "black", size = 7), axis.text = element_text(color = "black", size = 7),
  axis.title = element_text(color = "black", size = 7, lineheight = .95), plot.title = element_text(color = "black", size = 7.5, face = "bold", margin = margin(0, 0, 2, 0)),
  panel.grid.minor = element_blank(), panel.grid.major.y = element_blank(), panel.grid.major.x = element_line(color = "#E4E4E4", linewidth = .2),
  axis.line.x = element_line(color = "black", linewidth = .35), axis.ticks.x = element_line(color = "black", linewidth = .3),
  axis.ticks.length = unit(1.3, "pt"), plot.margin = margin(2, 6, 2, 2), legend.position = "none", plot.title.position = "panel",
  strip.text = element_blank(), strip.placement = "outside", panel.spacing.y = unit(6, "pt"))
gg_col <- function(d, key, lab_y, lims = NULL, brk = NULL, part = "") { m <- max(abs(d$x)) * 1.07
  if (is.null(lims)) lims <- c(-m, m)
  d <- d[d$x >= lims[1] & d$x <= lims[2], ]
  br <- if (is.null(brk)) scales::breaks_pretty(4)(lims) else brk; br <- br[is.finite(br) & br >= lims[1] & br <= lims[2]]
  gl <- merge(GG_BANDS[, c("half", "ymin", "ymax")], data.frame(x = br))          # x grid lines redrawn over the bands
  ggplot() + geom_rect(data = GG_BANDS, aes(xmin = -Inf, xmax = Inf, ymin = ymin, ymax = ymax), inherit.aes = FALSE, fill = BAND, color = NA) +
    geom_segment(data = gl, aes(x = x, xend = x, y = ymin, yend = ymax), inherit.aes = FALSE, color = "#E4E4E4", linewidth = .2) +
    (if (lims[1] < 0 && lims[2] > 0) geom_vline(xintercept = 0, color = "#777777", linewidth = .3) else list()) +
    geom_point(data = d[!d$pass, ], aes(x = x, y = y), shape = 16, size = .85, color = GG_COLS[["fail"]]) +
    (if (any(d$pass)) geom_point(data = d[d$pass, ], aes(x = x, y = y, color = col), shape = 16, size = .95) else list()) +
    facet_grid(half ~ ., scales = "free_y", space = "free_y", switch = "y") +
    geom_blank(data = GG_YB, aes(x = 0, y = y), inherit.aes = FALSE) +
    scale_y_continuous(breaks = GG_NL$y, labels = if (lab_y) parse(text = GG_NL$expr) else NULL, expand = c(0, 0)) +
    (if (any(d$pass)) scale_color_manual(values = GG_COLS) else list()) +
    scale_x_continuous(limits = lims, breaks = br, labels = function(v) sub("-", "−", format(v, drop0trailing = FALSE)), expand = c(0, 0), oob = scales::oob_keep) +
    labs(x = if (part == "stub") NULL else GG_XLAB[[key]], y = NULL, title = if (part == "main") " " else GG_PTITLE[[key]]) + gg_th() +
    # axis break: one "//" on the axis line in the gap between the stub and the main axis
    (if (part == "stub") coord_cartesian(clip = "off") else if (part == "main") list(coord_cartesian(clip = "off"),
      geom_text(data = data.frame(half = factor(GG_HALF[2], GG_HALF), x = lims[1] - diff(lims) * .012, y = GG_YMIN - .5 - gg_pad(GG_YMIN)),
                aes(x = x, y = y, label = "//"), inherit.aes = FALSE, size = 7 / .pt, family = "Arial", vjust = .5)) else list()) +
    (if (part == "stub") theme(plot.margin = margin(2, 0, 2, 2)) else if (part == "main") theme(plot.margin = margin(2, 6, 2, 0)) else theme()) +
    (if (lab_y) theme(strip.text.y.left = element_text(angle = 0, hjust = 1, vjust = .5, face = "bold", size = 7, lineheight = .95, margin = margin(0, 5, 0, 0)),
                      axis.text.y = element_text(size = 7, hjust = 1, margin = margin(0, 2, 0, 0))) else theme()) }
d1 <- GGV; d1$x <- d1$sc_signed_median; d3 <- GGV; d3$x <- d3$adult_mean
# intestine: one dot per regulon, the median of its intestinal study values (Elmentaite 2020 + 2021 the
# mean of its parts), as the SC-islet column
d2 <- GGV; d2$x <- apply(as.matrix(GGV[, paste0("gut_mean__", G_STUDY)]), 1, function(v) median(v[is.finite(v)]))
# adult islets (not used for the colors): axis broken at GG_AD_CUT; the values below it in a condensed stub
p3 <- if (any(d3$x < GG_AD_CUT, na.rm = TRUE)) {
  lo <- min(d3$x, na.rm = TRUE) * 1.05; hi <- max(d3$x, na.rm = TRUE) * 1.07
  stub <- gg_col(d3, "ad", FALSE, lims = c(lo, GG_AD_CUT), brk = -.3, part = "stub")
  main <- gg_col(d3, "ad", FALSE, lims = c(GG_AD_CUT, hi), brk = c(-.05, 0, .05), part = "main")
  wrap_plots(stub, main, nrow = 1, widths = c(.16, .84)) & theme(plot.margin = margin(2, 0, 2, 0))
} else gg_col(d3, "ad", FALSE)
GG_BLOCK <- (wrap_plots(list(gg_col(d1, "sc", TRUE), gg_col(d2, "gut", FALSE), p3), nrow = 1, widths = c(1, 1, 1)) & NOBG) +
  plot_annotation(theme = NOBG)

gg_tw <- function(s, pt = 7) {                      # text width (in): the wider of systemfonts' and the device's measure
  sf <- systemfonts::string_width(s, family = "Arial", size = pt, res = 600, weight = "normal") / 600
  dv <- vapply(s, function(x) {
    pushViewport(viewport(gp = gpar(fontfamily = "Arial", fontsize = pt, fontface = "plain")))
    v <- convertWidth(stringWidth(x), "in", valueOnly = TRUE); popViewport(); v
  }, 0, USE.NAMES = FALSE)
  pmax(sf, dv)
}

draw_panel <- function() {
  pushViewport(viewport(x = unit(.04, "in"), y = unit(H - G_TOP, "in"), width = unit(G_W, "in"), height = unit(G_H, "in"), just = c("left", "top")))
  # color key in one row under the plot block, centered in the slot
  hw <- gg_tw(GG_KEY_HEAD, PT)
  kw <- hw + .10 + sum(vapply(GG_KEYS, function(k) .14 + gg_tw(k[[1]], PT) + .16, 0))
  ky <- G_H - .04 + .1108333 - .16 - .01
  x <- (G_W - (kw - .16)) / 2
  text_at(GG_KEY_HEAD, x, ky, vp_h = G_H); x <- x + hw + .10
  for (k in GG_KEYS) {
    grid.points(unit(x + .05, "in"), unit(G_H - (ky + .005) - .045, "in"), pch = 16, size = unit(1.5, "mm"), gp = gpar(col = k[[2]]))
    text_at(k[[1]], x + .14, ky, vp_h = G_H); x <- x + .14 + gg_tw(k[[1]], PT) + .16
  }
  pushViewport(viewport(x = unit(0, "in"), y = unit(G_H - BLOCK_TOP, "in"), width = unit(G_W, "in"), height = unit(BLOCK_H, "in"),
                        just = c("left", "top")))
  grid.draw(patchworkGrob(GG_BLOCK)); popViewport()
  # (+)/(−) key under the row-group labels
  for (i in seq_along(GG_PM)) text_at(GG_PM[i], .06, c(2.235, 2.365)[i] - (2.64 - G_H) - .16, vp_h = G_H)
  popViewport()
}
