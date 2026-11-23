# Optional upstream technical-replicate QC is not part of the final113 workflow.
source("scripts/config.R")
out <- file.path(output_root, "technical_qc")
dir.create(out, recursive = TRUE, showWarnings = FALSE)
writeLines("Skipped: the final113 inputs contain one profile per participant. Upstream replicate QC requires the original libraries and internal mapping.", file.path(out, "SKIPPED.txt"))
