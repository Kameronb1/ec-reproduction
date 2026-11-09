# Figure S19B: factors with targets lower in SC-SIC. Per factor, its pooled targets among beta-side genes (lower in
# SC-SIC than SC-beta in at least 5 of 6 SC-islet studies), observed (point) against expected from group size alone
# (tick), on a log axis; at the right, observed (expected) and the test q. Key below.

TL <- targets_layout(); B <- TL$B
XL_B <- c(.2, 60); LX <- log10(XL_B)
pB <- ggplot(B) +
  geom_segment(aes(x = e, xend = k, y = y, yend = y), colour = "grey60", linewidth = .5) +
  geom_segment(aes(x = e, xend = e, y = y - TICK_H / PITCH / 2, yend = y + TICK_H / PITCH / 2), colour = "black", linewidth = LW_TICK, lineend = "butt") +
  geom_point(aes(x = k, y = y), shape = 21, fill = "grey25", colour = "grey25", size = PT_SZ, stroke = PT_ST) +
  geom_text(aes(x = 10^(LX[2] + DX / TL$PW * diff(LX)), y = y, label = txt), hjust = 0, size = PT / .pt, family = FONT, colour = "black") +
  geom_text(aes(x = 10^(LX[2] + TL$X2 / TL$PW * diff(LX)), y = y, label = txt2), hjust = 0, size = PT / .pt, family = FONT, colour = "black") +
  annotate("text", x = 10^(LX[2] + c(DX, TL$X2) / TL$PW * diff(LX)), y = nrow(B) + .8, label = TL$HDR$B, hjust = 0, size = PT / .pt, family = FONT, fontface = "bold") +
  targets_y(B) + scale_x_log10(limits = XL_B, breaks = c(.3, 1, 3, 10, 30), labels = c("0.3", "1", "3", "10", "30"), expand = c(0, 0), oob = scales::oob_keep) +
  coord_cartesian(clip = "off") + labs(x = "Pooled targets among beta-side genes (log scale)", y = NULL) + theme_targets
gB <- fix_widths(pB, nrow(B), TL); hB <- convertHeight(sum(gB$heights), "in", valueOnly = TRUE)
KEY_B <- list(list(type = "pt", label = "Observed", fill = "grey25"),
              list(type = "tk", label = "Expected from group size alone"))

# the panel's height follows its row count (2.63 and 1.91 in were sized for 9 and 5 rows)
draw_panel <- function() place(function() draw_targets(gB, hB, 1.91 + (nrow(B) - 5) * PITCH - .47, KEY_B, "Targets lower in SC-SIC", TL), 2858, 1140)
