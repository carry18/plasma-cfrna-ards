# 在项目根目录独立运行；输入为本次analysis目录中对应分析结果；输出为本图及源数据。
source("scripts/figures/style.R")

out_dir <- "results/supplementary_figures"
dir.create(out_dir, recursive = TRUE, showWarnings = FALSE)

save_sup <- function(plot, name, width = 183, height = 190) {
  base <- file.path(out_dir, name)
  ggsave(paste0(base, ".png"), plot,
    width = width, height = height,
    units = "mm", dpi = 300, device = ragg::agg_png, bg = "white"
  )
  grDevices::cairo_pdf(paste0(base, ".pdf"),
    width = width / 25.4,
    height = height / 25.4, family = "Arial"
  )
  print(plot)
  dev.off()
  svglite::svglite(paste0(base, ".svg"), width = width / 25.4, height = height / 25.4)
  print(plot)
  dev.off()
  ragg::agg_tiff(paste0(base, ".tiff"), width = width, height = height, units = "mm", res = 600, compression = "lzw")
  print(plot)
  dev.off()
}

# S2: network parameter choice and M02 resampling stability.
dsoft <- read.csv("results/analysis/wgcna/discovery_soft_threshold.csv")
asoft <- read.csv("results/analysis/wgcna/all_patient_soft_threshold.csv")
dsoft$network <- "50-patient discovery network"
asoft$network <- "All-patient network"
soft <- rbind(dsoft, asoft)

p2a <- ggplot(soft, aes(Power, SFT.R.sq, colour = network)) +
  geom_hline(yintercept = .8, linetype = 2, colour = "#999999", linewidth = .3) +
  geom_line(linewidth = .6) +
  geom_point(size = 1.3) +
  scale_colour_manual(values = c("50-patient discovery network" = "#087F8C", "All-patient network" = "#7B61A8"), name = NULL) +
  labs(x = "Soft-threshold power", y = expression("Scale-free topology " * R^2)) +
  theme(legend.position = "bottom")

checkpoint_dir <- Sys.getenv("CF_STABILITY_CHECKPOINT_DIR",
  "results/analysis/stability/checkpoints")
stability_files <- list.files(checkpoint_dir, pattern = "\\.rds$", full.names = TRUE)
stability_rows <- lapply(stability_files, function(f) {
  x <- readRDS(f)
  stopifnot(identical(x$status, "success"))
  m <- subset(x$matches, original == "M02")
  data.frame(design = paste0(switch(x$mode,
    fixed = "Fixed genes",
    reselected = "Gene re-selection",
    loo = "Leave-one-out"
  ), if (x$mode == "loo") "" else paste0(" (", x$fraction * 100, "% retained)")), mode = x$mode, fraction = x$fraction, iteration = x$iteration, jaccard = m$jaccard, recall = m$recall, precision = m$precision, score_correlation = m$score_cor_retained, direction_retained = as.numeric(m$direction_retained), stringsAsFactors = FALSE)
})
stopifnot(length(stability_files) == 350L)
stability <- do.call(rbind, stability_rows)
stopifnot(nrow(stability) == 350L)
stability$design <- factor(stability$design, levels = c("Fixed genes (80% retained)", "Fixed genes (90% retained)", "Gene re-selection (80% retained)", "Leave-one-out"))
write.csv(stability, file.path(out_dir, "FigureS2_resampling_source.csv"), row.names = FALSE)
summary_rows <- lapply(split(stability, stability$design), function(d) {
  z <- data.frame(design = as.character(d$design[1]), n = nrow(d))
  for (v in c("jaccard", "recall", "precision", "score_correlation")) {
    z[[paste0(v, "_median")]] <- median(d[[v]])
    z[[paste0(v, "_min")]] <- min(d[[v]])
    z[[paste0(v, "_max")]] <- max(d[[v]])
  }
  z$direction_retained <- sum(d$direction_retained)
  z
})
write.csv(do.call(rbind, summary_rows), "results/supplementary_tables/TableS3d_M02_resampling_summary.csv", row.names = FALSE)
write.csv(stability, "results/supplementary_tables/TableS3e_M02_resampling_draws.csv", row.names = FALSE)
levels(stability$design) <- c("Fixed\n80%", "Fixed\n90%", "Reselected\n80%", "Leave-one-out")
long <- rbind(data.frame(design = stability$design, metric = "Jaccard", value = stability$jaccard), data.frame(design = stability$design, metric = "Gene recall", value = stability$recall), data.frame(design = stability$design, metric = "Gene precision", value = stability$precision), data.frame(design = stability$design, metric = "Score correlation", value = stability$score_correlation))
long$metric <- factor(long$metric, levels = c("Jaccard", "Gene recall", "Gene precision", "Score correlation"))
p2b <- ggplot(long, aes(design, value, fill = design)) +
  geom_violin(width = .82, alpha = .25, colour = NA, trim = FALSE) +
  geom_boxplot(width = .18, outlier.shape = NA, linewidth = .25, colour = "#222222") +
  stat_summary(fun = median, geom = "point", size = 1.1, colour = "#222222") +
  facet_wrap(~metric, ncol = 2) +
  scale_fill_manual(values = c("#087F8C", "#3E9C8C", "#D17A22", "#7B61A8"), guide = "none") +
  scale_y_continuous(breaks = c(.5, .75, 1), labels = scales::percent_format(accuracy = 1)) +
  coord_cartesian(ylim = c(.45, 1.01)) +
  labs(x = NULL, y = "Resampling value") +
  theme(axis.text.x = element_text(angle = 25, hjust = 1, size = 6), strip.text = element_text(size = 7), panel.grid.major.y = element_line(colour = "#E5E5E5", linewidth = .25), legend.position = "none")
p2c <- ggplot(stability, aes(design, direction_retained, fill = design)) +
  stat_summary(fun = mean, geom = "col", width = .65, colour = "#222222", linewidth = .25) +
  geom_text(stat = "summary", fun = mean, aes(label = scales::percent(after_stat(y), accuracy = 1)), vjust = -0.35, size = 2.6) +
  scale_fill_manual(values = c("#087F8C", "#3E9C8C", "#D17A22", "#7B61A8"), guide = "none") +
  scale_y_continuous(limits = c(0, 1.12), breaks = c(0, .5, 1), labels = scales::percent_format(accuracy = 1)) +
  labs(x = NULL, y = "Clinical-group direction retained") +
  theme(axis.text.x = element_text(angle = 25, hjust = 1, size = 6), panel.grid.major.y = element_line(colour = "#E5E5E5", linewidth = .25))
save_sup(panelize(p2a, "A") / (panelize(p2b, "B") | panelize(p2c, "C")), "FigureS2", height = 210)
write.csv(soft, file.path(out_dir, "FigureS2_soft_threshold_source.csv"), row.names = FALSE)
