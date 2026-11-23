# Run from the project root. Only input data/reference paths belong here.
options(stringsAsFactors = FALSE)
output_root <- "results"
input_dir <- file.path(output_root, "inputs")
analysis_dir <- file.path(output_root, "analysis")
input_files <- list(
  counts = "data/private/count_matrix113.tsv",
  clinical = "data/private/sample_info113.csv",
  table1_clinical = "data/private/sample_info113.csv",
  msigdb_reference = "scripts/assets/msigdb_2026.1_Hs.rds",
  rr_signatures = "scripts/assets/published_RR_signatures.csv",
  geo_matrix = "data/external/GSE66890_series_matrix.txt.gz",
  geo_platform = "data/external/GPL6244.annot.gz"
)
# Optional overrides allow the same workflow to start in another checkout.
for (key in names(input_files)) input_files[[key]] <- Sys.getenv(paste0("CF_", toupper(key)), input_files[[key]])
group_labels <- c(A = "MV ARDS", B = "nonintubated ARDS", C = "Pneumonia without ARDS", D = "Healthy controls")
