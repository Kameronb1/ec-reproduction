# Figure S1A: marker-score histograms of each fitting group (top row) and of the two HIO groups the four-marker
# EC classifier scores but is not fitted on (bottom row), bars colored by how many of the four EC markers each
# cell has a read of, with the classifier's cutoff. Each column header gives the group's cells at or above the
# cutoff. Marker key in the empty slots right of the HIO row.

CUT <- jsonlite::fromJSON(result_path("intermediate/four_marker_fit.json"))$ec_core_cutoff
d <- read_table("fitting_cell_classes.tsv.gz"); hio <- read_table("hpsc_cell_classes.tsv.gz")
cells <- rbind(data.frame(arm = d$arm_id, score = d$ec_core_score, markers = d$markers_detected, fit = TRUE),
               data.frame(arm = hio$arm_id, score = hio$ec_core_score, markers = hio$markers_detected, fit = FALSE))
BW <- 0.10                                                     # bins 0.10 wide, one edge on the cutoff
span <- range(cells$score)
breaks <- c(rev(seq(CUT, span[1] - BW, by = -BW)), seq(CUT + BW, span[2] + BW, by = BW))
mids <- breaks[-length(breaks)] + BW / 2
cells$mid <- mids[as.integer(cut(cells$score, breaks = breaks, right = FALSE, include.lowest = TRUE))]
summ <- do.call(rbind, lapply(split(cells, cells$arm), function(x)
  data.frame(arm = x$arm[1], fit = x$fit[1], cells = nrow(x), above = sum(x$score >= CUT))))
summ <- summ[order(!summ$fit, -summ$above / summ$cells), ]     # fitting groups, then HIO; highest share above first
summ$study <- study(summ$arm); summ$type <- study_type(summ$arm)
summ$count <- sprintf("%s/%s ≥ cutoff", fmt(summ$above), fmt(summ$cells))
HDR <- c(fit = "Fitting groups (pooled in Fig. 2A)", hio = "HIO groups: scored by the four-marker EC classifier, not used to fit it")

one_hist <- function(a, ylab) {
  x <- cells[cells$arm == a, ]; s <- summ[summ$arm == a, ]
  tab <- aggregate(list(cells = rep(1L, nrow(x))), by = list(mid = x$mid, f = x$markers), FUN = length)
  tab$f <- factor(tab$f, levels = rev(names(MARKS)))
  ggplot(tab, aes(mid, cells, fill = f)) +
    geom_col(width = BW * .92) +
    geom_vline(xintercept = CUT, linetype = "dashed", linewidth = .35, color = "black") +
    scale_fill_manual(values = MARKS) +
    scale_x_continuous(breaks = -1:2, limits = c(-1.05, 2.05), expand = c(0, 0)) +
    scale_y_continuous(breaks = function(l) { b <- scales::breaks_pretty(4)(l); b[b == round(b)] }, expand = expansion(mult = c(0, .06))) +
    labs(x = NULL, y = if (ylab) "Cells" else NULL, title = s$study, subtitle = paste(s$type, s$count, sep = "\n")) +   # column header
    theme_kit("y") +
    theme(plot.title.position = "panel",
          plot.title = element_text(size = PT, face = "bold", hjust = 0, margin = margin(0, 0, 1.5, 0)),
          plot.subtitle = element_text(size = PT, hjust = 0, lineheight = 1, margin = margin(0, 0, 3, 0)),
          axis.line.y = element_line(color = "black", linewidth = .3), plot.margin = margin(2, 6, 2, 0),
          plot.background = element_blank())
}

W_A <- 7.3; H_A <- 3.90; X0 <- .03                             # panel width and height (in), left margin
HDR_H <- .19; XT_H <- .17; GAP <- .07; BOT <- .03              # row header, x title, gap between rows, bottom margin
ROW_H <- (H_A - .30 - 2 * HDR_H - 2 * XT_H - GAP - BOT) / 2
draw_a <- function() {
  vh <- convertHeight(unit(1, "npc"), "in", valueOnly = TRUE)
  arms <- summ$arm; first <- arms %in% c("E-MTAB-8901_paediatric", "GSE214852_WT")   # first column of each row
  gs <- lapply(seq_along(arms), function(i) ggplotGrob(one_hist(arms[i], ylab = first[i])))
  # align the histograms: same widths left of the plot (y title only in each row's first column), same heights
  col_of <- function(g, nm) unique(g$layout$l[g$layout$name == nm])
  pc <- col_of(gs[[1]], "panel"); yl <- col_of(gs[[1]], "ylab-l")
  wmax <- do.call(unit.pmax, lapply(gs, `[[`, "widths")); hmax <- do.call(unit.pmax, lapply(gs, `[[`, "heights"))
  for (i in seq_along(gs)) {
    gs[[i]]$widths[-pc] <- wmax[-pc]
    if (!first[i]) gs[[i]]$widths[yl] <- unit(0, "in")
    gs[[i]]$heights <- hmax
    gs[[i]]$layout$clip[gs[[i]]$layout$name %in% c("title", "subtitle")] <- "off"
  }
  left_all <- in_w(sum(gs[[1]]$widths[seq_len(pc - 1)])); right_all <- in_w(sum(gs[[1]]$widths[(pc + 1):length(gs[[1]]$widths)]))
  ylw <- in_w(gs[[1]]$widths[yl])
  colw <- (W_A - X0 - .01 - ylw) / 5                           # five columns; the first also holds the y title
  P <- colw - (left_all - ylw) - right_all                     # histogram width
  colx <- X0 + c(0, ylw + colw * 1:4)                          # left edge of each column's viewport
  colw_i <- c(colw + ylw, rep(colw, 4))
  panel_l <- X0 + left_all + c(0, colw * 1:4)                  # left edge of each histogram
  ytop <- function(y) unit(vh - y, "in")
  txt <- function(s, x, y, just = c("left", "top"), bold = FALSE, italic = FALSE)
    grid.text(s, x = unit(x, "in"), y = ytop(y), just = just,
              gp = gpar(fontfamily = FONT, fontsize = PT, fontface = if (bold) "bold" else if (italic) "italic" else "plain", col = "black"))
  header <- function(s, y) {                                   # bold row header with a line to the right edge
    txt(s, X0 + .02, y + .02, bold = TRUE)
    x1 <- X0 + .02 + str_in(s, bold = TRUE) + .08
    grid.lines(unit(c(x1, W_A - .03), "in"), ytop(c(y + .085, y + .085)), gp = gpar(col = "black", lwd = .5))
  }
  row <- function(idx, y) {
    for (j in seq_along(idx)) {
      pushViewport(viewport(x = unit(colx[j], "in"), y = ytop(y), width = unit(colw_i[j], "in"), height = unit(ROW_H, "in"), just = c("left", "top")))
      grid.draw(gs[[idx[j]]]); popViewport()
    }
    txt("Four-marker score", (panel_l[1] + panel_l[length(idx)] + P) / 2, y + ROW_H + .02, just = c("centre", "top"))
  }
  y <- 0
  header(HDR[["fit"]], y); y <- y + HDR_H
  row(1:5, y); y <- y + ROW_H + XT_H + GAP
  header(HDR[["hio"]], y); y <- y + HDR_H
  row(6:7, y)
  # key in the empty slots right of the HIO row
  kx <- panel_l[3]; ky <- y + .03
  txt("Markers detected (of 4)", kx, ky)
  txt("TPH1, SLC18A1, FEV, LMX1A", kx, ky + .145, italic = TRUE)
  for (m in 0:4) {
    sx <- kx + m * .40
    grid.rect(x = unit(sx, "in"), y = ytop(ky + .355), width = unit(.095, "in"), height = unit(.095, "in"), just = c("left", "centre"),
              gp = gpar(fill = MARKS[[as.character(m)]], col = NA))
    txt(as.character(m), sx + .14, ky + .355, just = c("left", "centre"))
  }
  grid.lines(unit(c(kx, kx + .16), "in"), ytop(c(ky + .53, ky + .53)), gp = gpar(col = "black", lty = 2, lwd = 1.1))
  txt(sprintf("Cutoff (%.4f)", CUT), kx + .22, ky + .53, just = c("left", "centre"))
}

draw_panel <- function() place(draw_a, 120, 48, 4380, 2196, W_A, H_A - .30)
