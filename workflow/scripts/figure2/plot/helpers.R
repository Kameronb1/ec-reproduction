# Drawing helpers shared by the Figure 2 panels: text size, theme, placement on the page, keys and study
# names. Sourced by workflow/scripts/common/draw_panel.R, which defines the page size W x H (inches).
# Every position is in inches from the page's top-left corner.

PT <- 7.5                         # body text, points
TM <- PT / ggplot2::.pt           # body text in ggplot's text size units

theme_fig2 <- function() theme_minimal(base_family = "Arial", base_size = PT) + theme(
  text = element_text(color = "black", size = PT), axis.text = element_text(color = "black", size = PT),
  axis.title = element_text(color = "black", size = PT), plot.title = element_text(color = "black", size = 8, face = "bold"),
  plot.subtitle = element_text(color = "black", size = PT), panel.grid.minor = element_blank(),
  panel.grid.major = element_line(color = "#E4E4E4", linewidth = .2), axis.line.x = element_line(color = "black", linewidth = .35),
  axis.ticks.x = element_line(color = "black", linewidth = .3), axis.ticks.length = unit(1.3, "pt"),
  plot.margin = margin(2, 3, 2, 2), legend.position = "none",
  plot.background = element_rect(fill = "white", colour = NA))

# Draw a ggplot, patchwork or gtable with its top-left corner at (x, y), w wide and h high.
draw_at <- function(p, x, y, w, h) {
  pushViewport(viewport(x = x / W, y = 1 - y / H, width = w / W, height = h / H, just = c("left", "top")))
  grid.draw(if (inherits(p, "gtable")) p else if (inherits(p, "patchwork")) patchworkGrob(p) else ggplotGrob(p))
  popViewport()
}

# Text with its top-left corner at (x, y).
text_at <- function(s, x, y, size = PT, bold = FALSE)
  grid.text(s, x = unit(x, "in"), y = unit(H - y, "in"), just = c("left", "top"),
            gp = gpar(fontfamily = "Arial", fontsize = size, fontface = if (bold) "bold" else "plain", col = "black"))

# Left and right edges (inches from the gtable's left edge) of the plot panels in gtable g drawn w inches wide.
panel_span <- function(g, w) {
  cw <- convertWidth(g$widths, "in", valueOnly = TRUE)   # fixed widths; the panels' stretchable widths count as 0
  k <- g$layout[grepl("^panel", g$layout$name), ]
  c(sum(cw[seq_len(min(k$l) - 1)]), w - sum(cw[-seq_len(max(k$r))]))
}

# An x-axis title: 7.5 pt, centered at x, baseline at y.
draw_xtitle <- function(s, x, y)
  grid.text(s, x = unit(x, "in"), y = unit(H - y, "in"), just = c("centre", "bottom"),
            gp = gpar(fontfamily = "Arial", fontsize = PT, col = "black"))

# Width of a string in inches on the open device.
str_in <- function(s, size = PT) {
  pushViewport(viewport(gp = gpar(fontfamily = "Arial", fontsize = size)))
  v <- convertWidth(stringWidth(s), "in", valueOnly = TRUE)
  popViewport()
  v
}

# A color square followed by its label.
key <- function(label, color, x, y) {
  grid.rect(x = unit(x, "in"), y = unit(H - y - .045, "in"), width = unit(.095, "in"), height = unit(.095, "in"),
            just = "left", gp = gpar(fill = color, col = NA))
  text_at(label, x + .14, y - .005)
}

# Study names: first author and year, short name and sample type (config/study_names.tsv).
SN <- read.delim(file.path(REPO, "config", "study_names.tsv"), check.names = FALSE, colClasses = "character",
                 na.strings = character(0), fileEncoding = "UTF-8")
rownames(SN) <- SN$unit_id

# The five studies with primary study labels drawn in B, C and E, in drawing order.
STUDY_ID <- c(gut_cell_atlas = "SF_ELMENTAITE2021_GUT_ATLAS", hickey = "SF_HICKEY2023_ADULT_EEC",
              egozi = "SF_EGOZI2021_DEVELOPMENTAL_EEC", gse238274 = "SF_INTESTINAL_GSE238274", hio = "SF_GSE214852_HIO_GRAFT")
FAMS <- names(STUDY_ID)
STUDY <- setNames(SN[STUDY_ID, "short_name"], FAMS)
STUDY["hickey"] <- "Hickey 2023 (nuclei)"          # single-nucleus RNA; D keeps the plain name
CONTEXT <- setNames(SN[STUDY_ID, "type_words"], FAMS)

# The four-marker classifier's cutoff on the marker score, from its fit (Figure 2A).
EC_CUTOFF <- jsonlite::fromJSON(result_path("intermediate/four_marker_fit.json", "figure2"))$ec_core_cutoff
