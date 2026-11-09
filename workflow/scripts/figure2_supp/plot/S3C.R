# Figure S3C: both classifiers on the test studies and HIO grafts with their counts thinned at random to 50%
# and 25% (three draws each; 100% is the full data): sensitivity (EC-labeled cells classified as EC), specificity
# (other terminal EEC not classified as EC) and ROC AUC, each draw a point and the mean of the draws a line. Columns on S3B's grid; a dotted
# line separates the tissue studies from Singh 2024 and McCauley 2023.

mr <- read_table("four_marker_depth.tsv"); mh <- read_table("hio_four_marker_depth.tsv")   # four-marker classifier
cl <- read_table("depth.tsv"); ch <- read_table("hio_classifier_depth.tsv")                # 20-gene classifier
mr <- rbind(mr, mh[names(mr)]); cl <- rbind(cl, ch[names(cl)])
M <- c(sensitivity = "Sensitivity", specificity = "Specificity", roc_auc = "ROC AUC")         # row order, top to bottom
ROW_TITLE <- c(Sensitivity = "EC labels called (sensitivity)", Specificity = "Other terminal EEC not called (specificity)",
               `ROC AUC` = "ROC AUC (EC labels versus other terminal EEC)")
long <- function(d, method) do.call(rbind, lapply(names(M), function(k) data.frame(method = method, study = d$study, fraction = d$fraction,
  replicate = d$replicate, metric = M[[k]], value = d[[k]])))
x <- rbind(long(mr, "Four-marker classifier"), long(cl, "20-gene classifier"))
x$method <- factor(x$method, levels = c("Four-marker classifier", "20-gene classifier"))
x$metric <- factor(x$metric, levels = unname(M)); x$st <- factor(x$study, levels = SCORED); x$pct <- 100 * x$fraction
mn <- aggregate(value ~ method + st + metric + pct, data = x, FUN = mean)                   # mean of the draws
HEAD_BOLD <- paste0(study(SCORED), ifelse(SCORED == "SF_HICKEY2023_ADULT_EEC", " (nuclei)", "")); HEAD_TYPE <- study_type(SCORED)
H_BUDGET_C <- 3.30; ROW_T_H <- .20; KEY_GAP_C <- .05; KEY_H_C <- .19   # panel height to fill, row titles, key gap and key (in)
pc3 <- ggplot(x, aes(pct, value, color = method)) +
  geom_rect(data = data.frame(metric = factor("Sensitivity", levels = unname(M))), aes(xmin = -Inf, xmax = Inf, ymin = -Inf, ymax = Inf),
            inherit.aes = FALSE, fill = "#F0F0F0", color = NA) +                            # sensitivity row shaded
  geom_line(data = mn, linewidth = .5) + geom_point(size = .9, alpha = .8) +
  facet_grid(metric ~ st) +
  scale_color_manual(values = c("Four-marker classifier" = RULE, "20-gene classifier" = PANEL), guide = "none") +
  scale_x_continuous(breaks = c(25, 50, 100), labels = function(v) paste0(v, "%"), limits = c(18, 107)) +
  scale_y_continuous(limits = c(0, 1), breaks = c(0, .5, 1), labels = c("0", "0.5", "1")) +
  labs(x = "Counts kept in the scored study", y = NULL) + theme_kit("y") +
  theme(strip.text = element_blank(), panel.spacing.x = unit(GRID_BC$S, "in"), panel.spacing.y = unit(ROW_T_H, "in"),
        plot.margin = margin(0, 0, 0, 0))
g <- ggplotGrob(pc3)
pr <- panel_rows(g)
avail <- H_BUDGET_C - TOP - HEAD_H - ROW_T_H - other_h(g) - KEY_GAP_C - KEY_H_C - PAD_B
PANEL_H <- floor(100 * avail * c(2, 1, 1) / 4) / 100                  # sensitivity row twice as tall as the flat rows
g <- fix_rows(fix_columns(g, GRID_BC$G, GRID_BC$P, GRID_BC$S, GRID_BC$W), PANEL_H)
g <- gtable::gtable_add_rows(g, unit(ROW_T_H, "in"), pos = pr[1] - 1) # a title row above the first panel row
pr <- panel_rows(g); pc <- panel_cols(g)
for (i in seq_along(M))                                               # row titles: in the added row and in the spacing rows
  g <- gtable::gtable_add_grob(g, textGrob(ROW_TITLE[[M[[i]]]], x = unit(.02, "in"), y = unit(4.5, "pt"), just = c("left", "bottom"),
                                           gp = gpar(fontfamily = FONT, fontsize = PT, fontface = "bold", col = "black")),
                               t = pr[i] - 1, l = 1, r = pc[length(pc)], clip = "off", name = paste0("row-title-", i))
g <- add_headers(g, lapply(seq_along(SCORED), function(i) list(l = i, r = i, bold = HEAD_BOLD[i], plain = HEAD_TYPE[i])), above = 1)
pc <- panel_cols(g); pr <- panel_rows(g); head_t <- unique(g$layout$t[grepl("^head-", g$layout$name)])
g <- gtable::gtable_add_grob(g, segmentsGrob(x0 = unit(.5, "npc"), x1 = unit(.5, "npc"), y0 = unit(0, "npc"), y1 = unit(1, "npc"),   # Egozi 2021 | Singh 2024
                                             gp = gpar(col = "#777777", lty = "dotted", lwd = .3 * .pt, lineend = "butt")),
                             t = head_t, b = pr[length(pr)], l = pc[3] + 1L, r = pc[3] + 1L, clip = "off", name = "study-break")
key_c <- list(list(list(label = "Four-marker classifier", glyph = "linepoint", col = RULE),
                   list(label = "20-gene classifier", glyph = "linepoint", col = PANEL),
                   list(label = "One draw (100%: full data)", glyph = "point", col = "black"),
                   list(label = "Mean of draws", glyph = "line", col = "black")))
H_G <- total_h(g)
H_C <- TOP + H_G + KEY_GAP_C + KEY_H_C + PAD_B

draw_panel <- function()
  place(stack_pieces(list(list(h = H_G, draw = g), list(h = KEY_GAP_C, draw = NULL), list(h = KEY_H_C, draw = draw_key(key_c)),
                          list(h = PAD_B, draw = NULL))),
        120, 3403, 4380, 1724, GRID_BC$W, H_C - (.30 + .17))
