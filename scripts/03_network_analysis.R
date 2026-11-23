# 输入：01整理的counts及02分析样本/表达信息。
# 输出：发现网络、全患者网络、模块评分及功能注释。
# 在项目根目录运行；前置步骤：01、02。
#!/usr/bin/env Rscript

source("scripts/internal/common.R")
suppressPackageStartupMessages({
  library(edgeR)
  library(limma)
  library(WGCNA)
})
options(stringsAsFactors = FALSE)
disableWGCNAThreads()
set.seed(20261016)

out <- "results/analysis/wgcna"
dir.create(out, recursive = TRUE, showWarnings = FALSE)
wr <- function(x, name) write.csv(x, file.path(out, paste0(name, ".csv")), row.names = FALSE)

d <- readRDS(file.path(root, "internal_processed.rds"))
md <- read.csv("results/analysis/included_samples_descriptive113.csv",
  check.names = FALSE
)
md$group <- factor(md$group, levels = c("A", "B", "C", "D"))
md$Batch <- factor(md$Batch)

dg <- DGEList(d$counts[, md$Sample, drop = FALSE])
dg <- dg[filterByExpr(dg, group = md$group), , keep.lib.sizes = FALSE]
dg <- calcNormFactors(dg)
design <- model.matrix(~ 0 + group + Batch, md)
v <- voom(dg, design, plot = FALSE)
expression <- removeBatchEffect(v$E, batch = md$Batch, design = model.matrix(~group, md))

discovery <- md$group %in% c("B", "C")
stopifnot(sum(discovery) == 50L)
gene_mad <- apply(expression[, discovery, drop = FALSE], 1, mad)
network_genes <- names(sort(gene_mad, decreasing = TRUE))[seq_len(min(5000L, sum(gene_mad > 0)))]
dat <- t(expression[network_genes, discovery, drop = FALSE])
stopifnot(goodSamplesGenes(dat, verbose = 0)$allOK)

powers <- c(1:10, 12, 14, 16, 18, 20)
sft <- pickSoftThreshold(dat, powerVector = powers, networkType = "signed", verbose = 0)
sft_table <- sft$fitIndices
eligible_power <- sft_table$Power[sft_table$SFT.R.sq >= 0.8 & sft_table$slope < 0]
power <- if (length(eligible_power)) min(eligible_power) else 6
wr(sft_table, "discovery_soft_threshold")

network <- blockwiseModules(
  dat,
  power = power, networkType = "signed", TOMType = "signed",
  minModuleSize = 40, mergeCutHeight = 0.25, numericLabels = TRUE,
  pamRespectsDendro = FALSE, maxBlockSize = 6000,
  randomSeed = 20261016, verbose = 2, nThreads = 1
)
labels <- setNames(network$colors, network_genes)
module_numbers <- sort(unique(labels[labels > 0]))
module_sets <- setNames(
  lapply(module_numbers, function(number) names(labels)[labels == number]),
  sprintf("M%02d", module_numbers)
)
wr(data.frame(module = names(module_sets), size = lengths(module_sets)), "module_sizes")

# Fixed discovery loadings projected to every retained independent profile.
scores <- matrix(NA_real_, nrow(md), length(module_sets),
  dimnames = list(md$Sample, names(module_sets))
)
membership <- list()
loadings <- list()
for (module in names(module_sets)) {
  genes <- module_sets[[module]]
  center <- rowMeans(expression[genes, discovery, drop = FALSE])
  spread <- apply(expression[genes, discovery, drop = FALSE], 1, sd)
  zmat <- t((expression[genes, , drop = FALSE] - center) / spread)
  pca <- prcomp(zmat[discovery, , drop = FALSE], center = FALSE, scale. = FALSE)
  loading <- pca$rotation[, 1]
  if (cor(pca$x[, 1], rowMeans(zmat[discovery, , drop = FALSE])) < 0) loading <- -loading
  score <- as.numeric(zmat %*% loading)
  score <- score / sd(score[discovery])
  scores[, module] <- score
  loadings[[module]] <- loading
  membership[[module]] <- data.frame(
    module = module,
    gene = genes,
    kME = as.numeric(cor(dat[, genes, drop = FALSE], score[discovery])),
    loading = loading
  )
}
membership <- do.call(rbind, membership)
wr(membership, "module_membership")
wr(data.frame(
  Sample = md$Sample, group = md$group, Batch = md$Batch,
  scores, check.names = FALSE
), "module_scores")

patients <- md$group %in% c("A", "B", "C")
patient_md <- droplevels(md[patients, ])
patient_design <- model.matrix(~ 0 + group + Batch, patient_md)
contrast_matrix <- makeContrasts(
  AB_C = (groupA + groupB) / 2 - groupC,
  B_C = groupB - groupC,
  A_C = groupA - groupC,
  A_B = groupA - groupB,
  levels = patient_design
)
fit <- eBayes(contrasts.fit(
  lmFit(t(scores[patients, , drop = FALSE]), patient_design),
  contrast_matrix
), robust = TRUE)
effects <- do.call(rbind, lapply(colnames(contrast_matrix), function(contrast) {
  tab <- topTable(fit, coef = contrast, number = Inf, sort.by = "none", confint = TRUE)
  data.frame(
    module = rownames(tab), contrast = contrast,
    estimate = tab$logFC, lower = tab$CI.L, upper = tab$CI.R,
    P = tab$P.Value, adjusted_P = tab$adj.P.Val,
    row.names = NULL
  )
}))
wr(effects, "module_group_effects")

# Functional and lung-cell reference annotations.
functional_sets <- readRDS("results/inputs/expanded_sets.rds")
cell_sets <- refs$cell
annotate_modules <- function(collection, kind) {
  collection <- lapply(collection, intersect, network_genes)
  collection <- collection[lengths(collection) >= 5 & lengths(collection) <= 500]
  rows <- list()
  for (module in names(module_sets)) {
    hits <- lapply(collection, intersect, module_sets[[module]])
    k <- lengths(hits)
    K <- lengths(collection)
    N <- length(network_genes)
    n <- length(module_sets[[module]])
    p <- phyper(k - 1, K, N - K, n, lower.tail = FALSE)
    rows[[module]] <- data.frame(
      module = module, kind = kind, pathway = names(collection), hits = k,
      P = p, adjusted_P = p.adjust(p, "BH"),
      genes = vapply(hits, paste, collapse = ";", FUN.VALUE = "")
    )
  }
  do.call(rbind, rows)
}
annotations <- rbind(
  annotate_modules(functional_sets, "functional"),
  annotate_modules(cell_sets, "cell")
)
wr(annotations, "module_annotations")

# M02 is the focal module label in this version's 50-patient discovery network.
# Its post hoc association/annotation are exploratory, not independently validated.
focal_module <- "M02"
stopifnot(focal_module %in% names(module_sets))

# Rebuild a separate network in all 79 independent patients and match the
# discovery focal module by gene overlap.
patient_expression <- expression[, patients, drop = FALSE]
patient_mad <- apply(patient_expression, 1, mad)
patient_genes <- names(sort(patient_mad, decreasing = TRUE))[seq_len(min(5000L, sum(patient_mad > 0)))]
patient_dat <- t(patient_expression[patient_genes, , drop = FALSE])
patient_sft <- pickSoftThreshold(patient_dat,
  powerVector = powers,
  networkType = "signed", verbose = 0
)
patient_sft_table <- patient_sft$fitIndices
patient_eligible <- patient_sft_table$Power[
  patient_sft_table$SFT.R.sq >= 0.8 & patient_sft_table$slope < 0
]
patient_power <- if (length(patient_eligible)) min(patient_eligible) else 6
wr(patient_sft_table, "all_patient_soft_threshold")
patient_network <- blockwiseModules(
  patient_dat,
  power = patient_power, networkType = "signed", TOMType = "signed",
  minModuleSize = 40, mergeCutHeight = 0.25, numericLabels = TRUE,
  pamRespectsDendro = FALSE, maxBlockSize = 6000,
  randomSeed = 20261016, verbose = 2, nThreads = 1
)
patient_labels <- setNames(patient_network$colors, patient_genes)
patient_numbers <- sort(unique(patient_labels[patient_labels > 0]))
patient_sets <- setNames(
  lapply(patient_numbers, function(number) names(patient_labels)[patient_labels == number]),
  sprintf("P%02d", patient_numbers)
)
focal_genes <- module_sets[[focal_module]]
patient_overlap <- do.call(rbind, lapply(names(patient_sets), function(module) {
  genes <- patient_sets[[module]]
  overlap <- length(intersect(genes, focal_genes))
  data.frame(
    module = module, size = length(genes), overlap_focal = overlap,
    recall_focal = overlap / length(focal_genes),
    precision_patient_module = overlap / length(genes),
    jaccard = overlap / length(union(genes, focal_genes))
  )
}))
patient_overlap <- patient_overlap[order(-patient_overlap$jaccard), ]
wr(patient_overlap, "all_patient_focal_overlap")
saveRDS(list(
  metadata = md, expression = expression, discovery = discovery,
  network_genes = network_genes, network = network, power = power,
  module_sets = module_sets, scores = scores, loadings = loadings,
  membership = membership, effects = effects, annotations = annotations,
  focal_module = focal_module,
  all_patient = list(
    power = patient_power, network = patient_network,
    sets = patient_sets, overlap = patient_overlap
  )
), file.path(out, "analysis_checkpoint_before_external.rds"))

capture.output(sessionInfo(), file = file.path(out, "sessionInfo_network.txt"))
