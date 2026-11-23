# 输入：01临床端点、02表达评分、03模块结果。
# 输出：30项临床关联、完整病例记录和血小板校正结果。
# 在项目根目录运行；前置步骤：01、02、03。
#!/usr/bin/env Rscript

source("scripts/internal/common.R")
suppressPackageStartupMessages({
  library(limma)
})

out <- "results/analysis/clinical"
dir.create(out, recursive = TRUE, showWarnings = FALSE)
wr <- function(x, name) write.csv(x, file.path(out, paste0(name, ".csv")), row.names = FALSE)

wgcna <- readRDS("results/analysis/wgcna/analysis_checkpoint_before_external.rds")
inputs <- readRDS("results/analysis/analysis_inputs.rds")
patient_md <- inputs$main$metadata
patient_expression <- inputs$main$expression

# Read clinical endpoints matched to the retained samples.
clinical <- read.csv(file.path(input_dir, "clinical_endpoints.csv"), check.names = FALSE)
stopifnot(!anyDuplicated(clinical$Sample))

bc_md <- droplevels(patient_md[patient_md$group %in% c("B", "C"), ])
stopifnot(nrow(bc_md) == 50L)
clinical <- clinical[match(bc_md$Sample, clinical$Sample), ]
stopifnot(identical(as.character(clinical$Sample), as.character(bc_md$Sample)))
clinical$group <- factor(bc_md$group, levels = c("C", "B"))
clinical$Batch <- factor(bc_md$Batch)

module_scores <- wgcna$scores[bc_md$Sample, , drop = FALSE]
reference_scores <- score_expr(
  wgcna$expression[, bc_md$Sample, drop = FALSE],
  refs$rr[c("RR_predicted40", "RR_validated7")]
)$scores
signatures <- cbind(module_scores, reference_scores)

endpoints <- c(
  "pf_ratio", "platelet", "il6_value", "wbc_value",
  "bilirubin_value", "creatinine_value"
)
clinical_results <- list()
for (signature in colnames(signatures)) {
  for (endpoint in endpoints) {
    data <- clinical
    data$score <- signatures[data$Sample, signature]
    data$outcome <- data[[endpoint]]
    if (endpoint %in% c("il6_value", "bilirubin_value", "creatinine_value")) {
      data$outcome <- log1p(data$outcome)
    }
    data <- data[is.finite(data$outcome) & is.finite(data$score), ]
    if (nrow(data) < 10L) next
    data$outcome <- as.numeric(scale(data$outcome))
    data$score <- as.numeric(scale(data$score))
    fit <- lm(outcome ~ score + group + Batch, data = data)
    tab <- summary(fit)$coefficients["score", ]
    ci <- confint(fit, "score")
    clinical_results[[paste(signature, endpoint, sep = "|")]] <- data.frame(
      signature = signature,
      endpoint = endpoint,
      n = nrow(data),
      estimate = unname(tab["Estimate"]),
      lower = unname(ci[1]),
      upper = unname(ci[2]),
      P = unname(tab["Pr(>|t|)"])
    )
  }
}
clinical_results <- do.call(rbind, clinical_results)
clinical_results$adjusted_P_30 <- p.adjust(clinical_results$P, method = "BH")
wr(clinical_results, "clinical_associations")

# Platelet-count adjustment for the focal module, fixed pathways and platelet markers.
fixed <- read.csv("results/analysis/fixed_score_values.csv",
  check.names = FALSE
)
fixed <- fixed[match(bc_md$Sample, fixed$Sample), ]
focal_module <- wgcna$focal_module
program_values <- data.frame(
  Sample = bc_md$Sample,
  M02_module = module_scores[, focal_module],
  fixed[, setdiff(names(fixed), c("Sample", "group", "Batch")), drop = FALSE],
  check.names = FALSE
)

marker_genes <- intersect(
  c("PPBP", "PF4", "HBB", "HBA1", "HBA2"),
  rownames(patient_expression)
)
marker_values <- t(patient_expression[marker_genes, bc_md$Sample, drop = FALSE])
program_values <- cbind(program_values, marker_values)
program_names <- setdiff(names(program_values), "Sample")

platelet_adjustment <- list()
for (program in program_names) {
  data <- clinical
  data$score <- program_values[match(data$Sample, program_values$Sample), program]
  data <- data[is.finite(data$score) & is.finite(data$platelet), ]
  if (nrow(data) < 10L) next
  data$score <- as.numeric(scale(data$score))
  data$platelet_z <- as.numeric(scale(data$platelet))
  base <- lm(score ~ group + Batch, data = data)
  adjusted <- lm(score ~ group + Batch + platelet_z, data = data)
  base_row <- summary(base)$coefficients["groupB", ]
  adjusted_row <- summary(adjusted)$coefficients["groupB", ]
  base_ci <- confint(base, "groupB")
  adjusted_ci <- confint(adjusted, "groupB")
  platelet_adjustment[[program]] <- data.frame(
    program = program,
    n = nrow(data),
    base_estimate = unname(base_row["Estimate"]),
    base_lower = unname(base_ci[1]),
    base_upper = unname(base_ci[2]),
    base_P = unname(base_row["Pr(>|t|)"]),
    adjusted_estimate = unname(adjusted_row["Estimate"]),
    adjusted_lower = unname(adjusted_ci[1]),
    adjusted_upper = unname(adjusted_ci[2]),
    adjusted_P = unname(adjusted_row["Pr(>|t|)"])
  )
}
platelet_adjustment <- do.call(rbind, platelet_adjustment)
wr(platelet_adjustment, "platelet_count_adjustment")

clinical_values <- clinical[, c("Sample", "group", "Batch", endpoints), drop = FALSE]
wr(clinical_values, "clinical_values_BC_unique50")
wr(data.frame(
  Sample = bc_md$Sample, group = bc_md$group, Batch = bc_md$Batch,
  program_values[, setdiff(names(program_values), "Sample"), drop = FALSE],
  check.names = FALSE
), "program_values_BC_unique50")

cat("Clinical endpoint availability\n")
print(vapply(
  clinical[, endpoints, drop = FALSE],
  function(x) sum(is.finite(x)), integer(1)
))
cat("\nFocal M02 associations\n")
print(subset(clinical_results, signature == focal_module))
cat("\nSelected platelet-count adjustment\n")
print(subset(
  platelet_adjustment,
  program %in% c(
    "M02_module", "GOBP_MEGAKARYOCYTE_DEVELOPMENT",
    "GOBP_PLATELET_AGGREGATION",
    "REACTOME_PLATELET_ACTIVATION_SIGNALING_AND_AGGREGATION",
    "PPBP", "PF4"
  )
))
