# 在项目根目录独立运行；输入为本次analysis目录中对应分析结果；输出为本图及源数据。
source("scripts/figures/common.R")
setup(5)
v <- read.csv("results/analysis/clinical/clinical_values_BC_unique50.csv")
v$M02 <- w79$scores[v$Sample, "M02"]
d <- subset(read.csv("results/analysis/clinical/clinical_associations.csv"), signature == "M02")
ep <- c(pf_ratio = "P/F ratio", platelet = "Platelet count", il6_value = "IL-6", wbc_value = "White blood cells", bilirubin_value = "Bilirubin", creatinine_value = "Creatinine")
d$label <- factor(paste0(ep[d$endpoint], " (n=", d$n, ")"), rev(paste0(ep[d$endpoint], " (n=", d$n, ")")))
wr(d, "A_clinical_effects")
pA <- ggplot(d, aes(estimate, label)) +
  geom_vline(xintercept = 0, linetype = 2, colour = "#AAB4BE") +
  geom_segment(aes(x = lower, xend = upper, yend = label), colour = module_pal["M02"], linewidth = .6) +
  geom_point(colour = module_pal["M02"], size = 2) +
  geom_text(aes(x = 1.05, label = sprintf("adj. P %.3g", adjusted_P_30)), hjust = 0, size = 2.4) +
  scale_x_continuous(limits = c(-1, 1.65), breaks = c(-1, -.5, 0, .5, 1)) +
  labs(x = "Standardized coefficient (95% CI)", y = NULL)
scatter <- function(endpoint) {
  z <- v[is.finite(v[[endpoint]]), ]
  z$score <- as.numeric(scale(z$M02))
  z$y <- as.numeric(scale(z[[endpoint]]))
  z$group <- factor(z$group)
  z$Batch <- factor(z$Batch)
  fit <- lm(y ~ score + group + Batch, z)
  beta <- coef(fit)["score"]
  z$xr <- resid(lm(score ~ group + Batch, z))
  z$yr <- resid(lm(y ~ group + Batch, z))
  line <- data.frame(xr = seq(min(z$xr), max(z$xr), length.out = 150))
  line$fit <- beta * line$xr
  delta <- qt(.975, df.residual(fit)) * abs(line$xr) * sqrt(vcov(fit)["score", "score"])
  line$lo <- line$fit - delta
  line$hi <- line$fit + delta
  wr(z, paste0(endpoint, "_points"))
  ggplot(z, aes(xr, yr)) +
    geom_ribbon(data = line, aes(x = xr, ymin = lo, ymax = hi), inherit.aes = FALSE, fill = "#F0E3EC") +
    geom_line(data = line, aes(xr, fit), colour = module_pal["M02"], linewidth = .5) +
    geom_point(aes(colour = group), size = 1.4, alpha = .8) +
    scale_color_manual(values = pal, labels = c(B = "Nonintubated ARDS", C = "Pneumonia without ARDS"), name = NULL) +
    labs(x = "M02 residual (standardized)", y = paste0(ep[endpoint], " residual"))
}
pB <- scatter("platelet")
pC <- scatter("pf_ratio")
pa <- ggplot(v[is.finite(v$platelet), ], aes(group, platelet, colour = group)) +
  geom_boxplot(width = .5, outlier.shape = NA, linewidth = .4) +
  geom_point(position = position_jitter(width = .12, seed = 7), size = 1) +
  scale_colour_manual(values = pal, guide = "none") +
  scale_x_discrete(limits = c("C", "B"), labels = c(C = "Pneumonia without ARDS", B = "Nonintubated ARDS")) +
  labs(x = NULL, y = "Platelets (10^9/L)")
cf <- read.csv("results/analysis/clinical/platelet_count_adjustment.csv")
wr(cf, "D_all_count_adjustment")
cf <- subset(cf, program %in% c("PPBP", "PF4", "HBB", "HBA1", "M02_module", "GOBP_PLATELET_AGGREGATION", "REACTOME_PLATELET_ACTIVATION_SIGNALING_AND_AGGREGATION", "GOBP_MEGAKARYOCYTE_DEVELOPMENT"))
wr(cf, "D_count_adjustment")
cf$label <- factor(wrap(pretty(cf$program), 25), rev(wrap(pretty(cf$program), 25)))
pb <- ggplot(cf, aes(y = label)) +
  geom_vline(xintercept = -log10(.05), linetype = 2, colour = "#AAB4BE") +
  geom_segment(aes(x = -log10(adjusted_P), xend = -log10(base_P), yend = label), colour = "#CCD2D6", linewidth = .6) +
  geom_point(aes(x = -log10(base_P), colour = "Group + batch"), size = 1.8) +
  geom_point(aes(x = -log10(adjusted_P), colour = "Group + batch + platelet count"), shape = 17, size = 1.8) +
  scale_colour_manual(values = c("Group + batch" = module_pal[["M02"]], "Group + batch + platelet count" = pal[["A"]]), name = NULL) +
  labs(x = "−log10 model P", y = NULL) +
  theme(axis.text.y = element_text(size = 5.4), legend.text = element_text(size = 5.2))
pD <- pa / pb + plot_layout(heights = c(.8, 1.2))
fig <- assemble(list(pA, pB, pC, pD), design = "AB\nCD", heights = c(1, 1.4))
finish(list(pA, pB, pC, pD), fig, 5, 225, "Current 79-patient results; accepted clinical panel layout retained.")
