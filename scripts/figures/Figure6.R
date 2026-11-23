# 在项目根目录独立运行；输入为本次analysis目录中对应分析结果；输出为本图及源数据。
source("scripts/figures/style.R")
source("scripts/internal/common.R")

a <- readRDS("results/analysis/wgcna/analysis.rds")
external <- readRDS(file.path(root, "public_context_processed.rds"))$GSE66890
focal_genes <- a$module_sets[[a$focal_module]]
matrix <- external$expr[intersect(focal_genes, rownames(external$expr)), , drop = FALSE]
keep <- rowSums(is.finite(matrix)) == ncol(matrix) & apply(matrix, 1, sd, na.rm = TRUE) > 0
matrix <- matrix[keep, , drop = FALSE]
scores <- data.frame(
  score = as.numeric(scale(colMeans(t(scale(t(matrix)))))),
  group = factor(ifelse(external$metadata$ARDS == 1, "Sepsis with ARDS", "Sepsis without ARDS"),
    levels = c("Sepsis without ARDS", "Sepsis with ARDS")
  )
)
write_source(scores, 6, "A_GSE66890_scores")
pA <- ggplot(scores, aes(group, score, colour = group)) +
  geom_boxplot(width = .48, outlier.shape = NA, linewidth = .45) +
  geom_point(position = position_jitter(width = .12, seed = 4), size = 1.05, alpha = .75) +
  scale_colour_manual(values = c(
    "Sepsis without ARDS" = unname(group_cols["Pneumonia without ARDS"]),
    "Sepsis with ARDS" = unname(group_cols["nonintubated ARDS"])
  ), guide = "none") +
  scale_x_discrete(labels = c("Sepsis without\nARDS\nn = 28", "Sepsis with\nARDS\nn = 29")) +
  labs(x = NULL, y = "Standardized M02 score")

effect <- a$external$effect
effect$label <- paste0("M02\n", effect$measured_genes, " genes")
write_source(effect, 6, "B_external_effect")
pB <- ggplot(effect, aes(estimate, label)) +
  geom_vline(xintercept = 0, linetype = 2, colour = "#999999", linewidth = .35) +
  geom_segment(aes(x = lower, xend = upper, yend = label), linewidth = .8, colour = "#B05A92") +
  geom_point(size = 2.8, colour = "#B05A92") +
  geom_text(aes(label = sprintf("P = %.3g", P)),
    x = effect$upper + .25,
    hjust = 1, size = 2.6
  ) +
  scale_x_continuous(expand = expansion(mult = c(.08, .28))) +
  labs(x = "ARDS vs sepsis score difference (95% CI)", y = NULL)

kme <- a$external$kME
write_source(kme, 6, "C_gene_module_correspondence")
pearson <- cor(kme$internal_kME, kme$external_kME)
spearman <- cor(kme$internal_kME, kme$external_kME, method = "spearman")
pC <- ggplot(kme, aes(internal_kME, external_kME)) +
  geom_hline(yintercept = 0, colour = "#D7DDE0", linewidth = .3) +
  geom_vline(xintercept = 0, colour = "#D7DDE0", linewidth = .3) +
  geom_point(size = .75, alpha = .42, colour = "#2F5D7E") +
  geom_smooth(method = "lm", se = FALSE, linewidth = .55, colour = "#B5523B") +
  annotate("text",
    x = Inf, y = Inf,
    label = sprintf("Pearson r = %.3f\nSpearman rho = %.3f", pearson, spearman),
    hjust = 1.08, vjust = 1.2, size = 2.6
  ) +
  labs(x = "Internal plasma cfRNA kME", y = "GSE66890 whole-blood kME")

fig <- (panelize(pA, "A") | panelize(pB, "B")) /
  panelize(pC, "C") + plot_layout(heights = c(.95, 1.15), widths = c(.9, 1.1))
save_figure(fig, 6, height_mm = 165)
writeLines(
  "Figure 6. External expression and structural correspondence of duplicate-resolved M02 in GSE66890.",
  file.path(out_root, "Figure6", "legend.txt")
)
