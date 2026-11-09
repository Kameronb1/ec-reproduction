# Figure 4A: the eleven rows that C-F compare, as a table: each row's study and cell line (bold study name), its start and
# end (the change) and how many samples each end has (from rows_called.tsv); eight grafts, then three longer in vitro
# cultures, with grey bands on alternate rows.

A_X0 <- .06; A_X1 <- .06 + 2.34                                # left and right edges of the bands
A_TEXT <- .10; A_NUM_RIGHT <- 2.38                             # left edge of the text, right edge of the sample counts
A_INDENT <- .10                                                # the rows and their bands sit this far right of the group headers
A_CHANGE_INDENT <- .10                                         # each change line sits this much further right, under its study
A_TOP <- .80 - 111 / 600; A_BOTTOM <- 3.35 + 122 / 600         # top of the first row, bottom of the last
A_HEAD_SLOT <- .17; A_HEAD_PAD <- .04                          # slot for 'Continued culture (in vitro)'; header to its row

nsamp <- function(L) vapply(L, function(s) length(unlist(strsplit(s, "+", fixed = TRUE))), 0L)   # a pooled end ("a+b+c") counts its samples
A <- ROWS; A$n <- sprintf("%d → %d", nsamp(REF[A$id]), nsamp(TEST[A$id]))
A$study_bold <- sub(",.*$", "", A$study); A$study_rest <- substring(A$study, nchar(A$study_bold) + 1L)
# a row with the same study as the row above, in the same group, shows only its change, centered in the row
A$repeat_study <- c(FALSE, A$study[-1] == A$study[-nrow(A)]) & seq_len(nrow(A)) != length(GRAFT) + 1L
# Eleven rows at one pitch below the header. Each row's two lines (study; change) are centered in the row by the tallest
# ascent and deepest descent of their strings; the sample counts are centered on the row.
L1 <- 7.5 / 72                                                 # line pitch
pitch <- (A_BOTTOM - A_TOP - A_HEAD_SLOT) / nrow(A)
row_top <- A_TOP + (seq_len(nrow(A)) - 1) * pitch + ifelse(seq_len(nrow(A)) > length(GRAFT), A_HEAD_SLOT, 0)
a1 <- max(vapply(A$study_bold, asc_in, 0, size = 7, face = "bold"), vapply(A$study_rest, asc_in, 0, size = 7))
d2 <- max(vapply(A$change, dsc_in, 0, size = 7)); pad <- (pitch - (a1 + L1 + d2)) / 2
an <- max(vapply(A$n, asc_in, 0, size = 7)); dn <- max(vapply(A$n, dsc_in, 0, size = 7))
base1 <- row_top + pad + a1; base2 <- row_top + pad + a1 + L1; num_base <- row_top + pitch / 2 + (an - dn) / 2
head_base2 <- A_TOP - A_HEAD_PAD - dsc_in("Grafts (in mice)", 7, "bold"); head_base1 <- head_base2 - L1
culture_base <- row_top[length(GRAFT) + 1L] - A_HEAD_PAD - dsc_in("Continued culture (in vitro)", 7, "bold")

# text on the baseline y with each arrow drawn bold; right-justified text is placed by its total width
text_arrow <- function(s, x, y, just = "left") {
  parts <- strsplit(s, "→", fixed = TRUE)[[1]]; if (endsWith(s, "→")) parts <- c(parts, "")
  runs <- c(rbind(parts, c(rep("→", length(parts) - 1L), "")))[seq_len(2L * length(parts) - 1L)]
  faces <- rep(c("plain", "bold"), length.out = length(runs)); w <- mapply(str_in, runs, 7, faces)
  if (just == "right") x <- x - sum(w)
  # Arial has no arrow glyph (the fallback font has no bold), so each arrow is drawn: a 1.1 pt line with a closed head,
  # in the glyph's box, at the height of the glyph's shaft
  for (k in seq_along(runs)) { x0 <- x + sum(w[seq_len(k - 1L)])
    if (runs[k] == "→") grid.lines(x = unit(x0 + w[k] * c(.10, .90), "in"), y = unit(rep(H - y + .034, 2), "in"),
                                     arrow = arrow(angle = 25, length = unit(.032, "in"), type = "closed"),
                                     gp = gpar(col = "black", fill = "black", lwd = 1.1, linejoin = "mitre"))
    else if (nzchar(runs[k])) text_bl(runs[k], x0, y) }
}
draw_panel <- function() shifted_up(ROW_UP[["ab"]], {
  for (i in which(A$y %in% BAND_ROWS)) grid.rect(x = unit(A_X0 + A_INDENT, "in"), y = unit(H - row_top[i], "in"), width = unit(A_X1 - A_X0 - A_INDENT, "in"),
                                                 height = unit(pitch, "in"), just = c("left", "top"), gp = gpar(fill = BAND_GRAY, col = NA))
  text_bl("Samples", A_NUM_RIGHT, head_base1, just = "right"); text_arrow("start → end", A_NUM_RIGHT, head_base2, just = "right")
  text_bl("Grafts (in mice)", A_TEXT, head_base2, face = "bold"); text_bl("Continued culture (in vitro)", A_TEXT, culture_base, face = "bold")
  for (i in seq_len(nrow(A))) {
    if (A$repeat_study[i]) text_arrow(A$change[i], A_TEXT + A_INDENT + A_CHANGE_INDENT, num_base[i]) else {
      text_bl(A$study_bold[i], A_TEXT + A_INDENT, base1[i], face = "bold"); text_bl(A$study_rest[i], A_TEXT + A_INDENT + str_in(A$study_bold[i], 7, "bold"), base1[i])
      text_arrow(A$change[i], A_TEXT + A_INDENT + A_CHANGE_INDENT, base2[i]) }
    text_arrow(A$n[i], A_NUM_RIGHT, num_base[i], just = "right")
  }
})
