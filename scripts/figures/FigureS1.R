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

# S1: complete Hallmark context for all four patient contrasts.
h <- subset(
  read.csv("results/analysis/hallmark_all_contrasts.csv"),
  cohort == "unique_patient79"
)
labels <- c(
  AB_C = "ARDS overall vs pneumonia without ARDS", B_C = "nonintubated ARDS vs pneumonia without ARDS",
  A_C = "MV ARDS vs pneumonia without ARDS", A_B = "MV vs nonintubated ARDS"
)
h$comparison <- factor(labels[h$contrast], levels = unname(labels))
h$pathway_label <- clean_pathway(h$pathway)
ord <- aggregate(abs(NES) ~ pathway_label, h, max)
ord <- ord$pathway_label[order(ord$`abs(NES)`)]
h$pathway_label <- factor(h$pathway_label, levels = ord)
h$significant <- h$padj < 0.05
p1 <- ggplot(h, aes(NES, pathway_label, colour = NES, shape = significant)) +
  geom_vline(xintercept = 0, linetype = 2, linewidth = .25, colour = "#999999") +
  geom_point(size = 1.25) +
  facet_wrap(~comparison, nrow = 1) +
  scale_colour_gradient2(
    low = signal_cols["lower"], mid = "#D8D8D8",
    high = signal_cols["higher"], midpoint = 0, name = "NES"
  ) +
  scale_shape_manual(
    values = c(`TRUE` = 16, `FALSE` = 1), name = NULL,
    labels = c("Adjusted P >= 0.05", "Adjusted P < 0.05")
  ) +
  labs(x = "Normalized enrichment score", y = NULL) +
  theme(
    axis.text.y = element_text(size = 4.4), strip.text = element_text(size = 6.8),
    legend.position = "bottom"
  )
save_sup(panelize(p1, "A"), "FigureS1", height = 220)
write.csv(h, file.path(out_dir, "FigureS1_source.csv"), row.names = FALSE)
