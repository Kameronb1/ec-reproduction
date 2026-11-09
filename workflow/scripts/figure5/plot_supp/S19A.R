# Figure S19A: factors with targets higher in SC-SIC. Per factor, the pooled targets of its reproducible regulons that
# are in G2 or G3, split into G2 (left) and G3 (right) as shares of the bar, against the split expected from group size
# alone (dashed line: G2's share of all G2 + G3 genes); at the right, the counts and the q of the lean test toward the
# group the factor has more of than expected. Key below.

TL <- targets_layout(); A <- TL$A
G2COL <- "#A89A2E"; G3COL <- "#B4553F"              # gene-group colors as in Figure 3
R2 <- mean(A$e2 / A$n)                              # G2's share of all G2 + G3 genes
pA <- ggplot(A) +
  geom_rect(aes(xmin = 0, xmax = g2 / n, ymin = y - .62 / 2, ymax = y + .62 / 2), fill = G2COL) +
  geom_rect(aes(xmin = g2 / n, xmax = 1, ymin = y - .62 / 2, ymax = y + .62 / 2), fill = G3COL) +
  geom_vline(xintercept = R2, linetype = "dashed", linewidth = LW_DASH, colour = "black") +
  geom_text(aes(x = 1 + DX / TL$PW, y = y, label = txt), hjust = 0, size = PT / .pt, family = FONT, colour = "black") +
  geom_text(aes(x = 1 + TL$X2 / TL$PW, y = y, label = txt2), hjust = 0, size = PT / .pt, family = FONT, colour = "black") +
  annotate("text", x = 1 + c(DX, TL$X2) / TL$PW, y = nrow(A) + .8, label = TL$HDR$A, hjust = 0, size = PT / .pt, family = FONT, fontface = "bold") +
  targets_y(A) + scale_x_continuous(limits = c(0, 1), breaks = c(0, .5, 1), labels = c("0%", "50%", "100%"), expand = c(0, 0), oob = scales::oob_keep) +
  coord_cartesian(clip = "off") + labs(x = "Share of each factor's G2 + G3 targets", y = NULL) + theme_targets
gA <- fix_widths(pA, nrow(A), TL); hA <- convertHeight(sum(gA$heights), "in", valueOnly = TRUE)
KEY_A <- list(list(type = "sq", label = "G2: EC genes SC-SIC have", col = G2COL),
              list(type = "sq", label = "G3: SC-SIC's own program", col = G3COL),
              list(type = "ln", label = "G2's share of all G2 + G3 genes"))

# the panel's height follows its row count (2.63 and 1.91 in were sized for 9 and 5 rows)
draw_panel <- function() place(function() draw_targets(gA, hA, 2.63 + (nrow(A) - 9) * PITCH - .47, KEY_A, "Targets higher in SC-SIC", TL), 199, 2436)
