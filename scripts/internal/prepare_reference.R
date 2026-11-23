# Optional: rebuild reference definitions from MSigDB collection files.
cache <- Sys.getenv("CF_MSIGDB_DIR")
if (!nzchar(cache)) stop("Set CF_MSIGDB_DIR to the directory containing the MSigDB collection RDS files.")
load_collection <- function(id) readRDS(file.path(cache, paste0("msigdb.2026.1.Hs.", id, ".rds")))
h <- load_collection("H")
c2 <- load_collection("C2")
c5 <- load_collection("C5")
c8 <- load_collection("C8")
x <- rbind(c2[c2$gs_subcollection %in% c("CP:REACTOME", "CP:KEGG_LEGACY", "CP:KEGG_MEDICUS"), ], c5[c5$gs_subcollection == "GO:BP", ])
x$database <- sub("CP:", "", x$gs_subcollection)
definitions <- unique(x[, c("gs_name", "database", "gs_description", "gs_exact_source", "db_version")])
names(definitions)[1] <- "pathway"
c8 <- c8[grepl("^TRAVAGLINI_LUNG_", c8$gs_name), ]
reference <- list(version = "MSigDB 2026.1.Hs", hallmark_definition = data.frame(gs_name = h$gs_name, gene_symbol = h$db_gene_symbol), expanded = lapply(split(x$db_gene_symbol, x$gs_name), unique), cell = lapply(split(c8$db_gene_symbol, c8$gs_name), unique), definitions = definitions)
saveRDS(reference, "scripts/assets/msigdb_2026.1_Hs.rds", compress = "xz")
