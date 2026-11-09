# Figure 5B-D key, two centered rows under B, C and D: each gene group G1-G3 as B and C draw it (its point on its light
# bar); then the marks, in the order they appear from B to D: the group's median gene, the middle half of the group's
# genes (B's and C's bars), C's graft-median diamond and D's grey bar (middle 95% of the random gene sets). A dotted
# outline sets the key apart from the plots.

KEY1_Y <- C2 + H2 + .005; KEY2_Y <- KEY1_Y + .125
draw_panel <- function() {
  w1 <- row_keys(lapply(1:3, function(k) list(type = "pb", label = GL[k], fill = GCOL[[k]], bar = GCOL_LIGHT[[k]], plwd = PT_LWD)),
           .10, 7.3, KEY1_Y, gap = .2)
  w2 <- row_keys(list(list(type = "pt", pch = 21, label = "Group's median gene", col = "black", fill = "grey75", plwd = PT_LWD),
                list(type = "bars3", label = "Middle half of the group's genes", col = unname(GCOL_LIGHT), h = .09),
                list(type = "pt", pch = 23, label = "Graft median", col = "black", fill = "white", plwd = 1.2),
                list(type = "vbar", label = "Middle 95% of 10,000 random gene sets", col = BAND, w = .85 * ggplot2::.pt / 96, h = .09)),
           .10, 7.3, KEY2_Y, gap = .2)
  key_frame(.10 + 7.3 / 2, max(w1, w2) + .2, KEY1_Y - .05, KEY2_Y + .135)
}
