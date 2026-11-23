# 输入：最终113份原始整数count matrix、对应113行样本/临床信息及参考注释。
# 输出：inputs/内的原始counts、分析样本名单、临床端点和参考集合。
# 在项目根目录运行；前置步骤：无。
# Find the project root for Rscript, source(), or interactive execution.
local({
  args <- commandArgs(trailingOnly = FALSE)
  cli_file <- sub("^--file=", "", args[grepl("^--file=", args)])
  source_files <- unlist(lapply(sys.frames(), function(frame) frame$ofile))
  starts <- c(getwd(), dirname(c(cli_file, source_files)))
  found <- FALSE
  for (start in starts) {
    candidate <- normalizePath(start, mustWork = FALSE)
    repeat {
      if (file.exists(file.path(candidate, "scripts/config.R"))) {
        setwd(candidate)
        found <- TRUE
        break
      }
      parent <- dirname(candidate)
      if (identical(parent, candidate)) break
      candidate <- parent
    }
    if (found) break
  }
  if (!found) stop("Cannot locate project root. Set working directory to the cf project and retry.")
})
source("scripts/config.R")
missing <- unlist(input_files)[!file.exists(unlist(input_files))]
if (length(missing)) stop("Missing inputs:\n", paste(missing, collapse = "\n"))
for (d in c(input_dir, analysis_dir, file.path(output_root, "supplementary_tables"), file.path(output_root, "supplementary_figures"), file.path(output_root, "logs"))) dir.create(d, recursive = TRUE, showWarnings = FALSE)
jsonlite::write_json(input_files, file.path(input_dir, "source_paths.json"), auto_unbox = TRUE, pretty = TRUE)
write.csv(data.frame(input = names(input_files), path = unlist(input_files), md5 = unname(tools::md5sum(unlist(input_files)))), file.path(input_dir, "input_checksums.csv"), row.names = FALSE)

# Host gene symbols are collapsed by summing counts, before any cohort-specific filtering.
path <- input_files$counts
if (grepl("[.]rds$", path, ignore.case = TRUE)) {
  raw <- as.matrix(readRDS(path))
} else {
  tab <- read.table(path, header = TRUE, sep = if (grepl("[.]csv$", path, ignore.case = TRUE)) "," else "\t", check.names = FALSE, quote = '"', comment.char = "")
  raw <- as.matrix(tab[, -1, drop = FALSE])
  rownames(raw) <- tab[[1]]
}
storage.mode(raw) <- "numeric"
stopifnot(!anyDuplicated(colnames(raw)), all(is.finite(raw)), all(raw >= 0), all(raw == floor(raw)))
if (is.null(rownames(raw)) || is.null(colnames(raw))) stop("Count matrix requires gene rows and sample columns")
stopifnot(ncol(raw) == 113L)
stopifnot(setequal(colnames(raw), sprintf("CFPL%03d", seq_len(113L))))
# Accept gene-symbol rows or the original ENSG|length|symbol|type row identifiers.
ids <- rownames(raw)
if (any(grepl("|", ids, fixed = TRUE))) {
  parts <- strsplit(ids, "|", fixed = TRUE)
  if (any(lengths(parts) != 4L)) stop("Mixed/invalid compound gene identifiers")
  symbols <- vapply(parts, `[`, character(1), 3L)
} else {
  if (any(grepl("^ENSG[0-9]+([.][0-9]+)?$", ids))) stop("Map Ensembl-only row names to gene symbols using the matching genome annotation first")
  symbols <- ids
}
keep <- !is.na(symbols) & nzchar(symbols)
counts <- rowsum(raw[keep, , drop = FALSE], symbols[keep], reorder = FALSE)
saveRDS(raw, file.path(input_dir, "feature_counts.rds"))
stopifnot(!anyDuplicated(rownames(counts)), all(is.finite(counts)), all(counts >= 0))
cl <- read.csv(input_files$clinical, check.names = FALSE)
stopifnot(all(c("Sample", "group", "Batch") %in% names(cl)))
stopifnot(!anyDuplicated(cl$Sample))
stopifnot(nrow(cl) == 113L, setequal(cl$Sample, colnames(counts)),
          !anyNA(cl$group), all(cl$group %in% c("A", "B", "C", "D")),
          !anyNA(cl$Batch), all(nzchar(cl$Batch)))
cl <- cl[match(colnames(counts), cl$Sample), ]
stopifnot(identical(as.character(cl$Sample), colnames(counts)))
md <- cl[, c("Sample", "group", "Batch")]
md$group <- factor(md$group, levels = c("A", "B", "C", "D"))
md$Batch <- factor(md$Batch)
saveRDS(list(counts = counts, cl = md), file.path(input_dir, "internal_processed.rds"))
# All inputs already represent the final analysis population; do not select again.
included <- md
stopifnot(identical(as.integer(table(included$group)), c(29L, 17L, 33L, 34L)))
write.csv(included, file.path(input_dir, "analysis_samples.csv"), row.names = FALSE)
write.csv(data.frame(Sample = md$Sample, included = TRUE,
  reason = "provided_final113_input"),
  file.path(input_dir, "sample_selection_audit.csv"), row.names = FALSE)

# Fixed clinical input: no Excel parsing, identity linking or replicate transfer.
clinical <- cl[cl$group %in% c("A", "B", "C"), , drop = FALSE]
required <- c("Sample", "group", "Batch", "pf_ratio", "il6_value", "wbc_value",
              "platelet", "creatinine_value", "bilirubin_value", "hco3", "il6_censored")
stopifnot(all(required %in% names(clinical)), !anyDuplicated(clinical$Sample))
patients <- included[included$group %in% c("A", "B", "C"), , drop = FALSE]
stopifnot(setequal(clinical$Sample, patients$Sample))
clinical <- clinical[match(patients$Sample, clinical$Sample), required, drop = FALSE]
stopifnot(identical(as.character(clinical$Sample), as.character(patients$Sample)),
          all(clinical$group == as.character(patients$group)), all(clinical$Batch == as.character(patients$Batch)))
for (field in setdiff(required, c("Sample", "group", "Batch", "il6_censored"))) {
  stopifnot(is.numeric(clinical[[field]]), all(is.na(clinical[[field]]) | is.finite(clinical[[field]])))
}
stopifnot(is.logical(clinical$il6_censored), !anyNA(clinical$il6_censored))
write.csv(clinical, file.path(input_dir, "clinical_endpoints.csv"), row.names = FALSE)

# Locked, published annotation inputs; no data-dependent gene-set discovery.
annotation <- readRDS(input_files$msigdb_reference)
hall <- annotation$hallmark_definition
write.csv(hall, file.path(input_dir, "hallmark_definition.csv"), row.names = FALSE)
expanded <- annotation$expanded
cell <- annotation$cell
definitions <- annotation$definitions
rr <- read.csv(input_files$rr_signatures)
rrsets <- lapply(split(rr$gene, rr$signature), unique)
rrsets <- rrsets[c("RR_predicted40", "RR_validated7")]
stopifnot(length(rrsets[[1]]) == 40L, length(rrsets[[2]]) == 7L)
refs <- list(hallmark = lapply(split(hall$gene_symbol, hall$gs_name), unique), expanded = expanded, cell = cell, rr = rrsets)
saveRDS(refs, file.path(input_dir, "reference_sets.rds"))
saveRDS(expanded, file.path(input_dir, "expanded_sets.rds"))
write.csv(definitions, file.path(input_dir, "pathway_definitions.csv"), row.names = FALSE)
saveRDS(list(cellsets = cell, rrsets = rrsets), file.path(input_dir, "reference_signatures.rds"))
capture.output(sessionInfo(), file = file.path(input_dir, "sessionInfo.txt"))
