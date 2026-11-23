source("scripts/internal/common.R")
md <- read.csv(file.path(input_dir, "GSE66890_metadata.csv"), check.names = FALSE)
names(md)[names(md) == "!Sample_geo_accession"] <- "sample_id"
names(md)[names(md) == "!Sample_title"] <- "title"
md <- md[md[["disease state"]] %in% c("sepsis alone", "sepsis with ARDS"), ]
md$ARDS <- as.integer(md[["disease state"]] == "sepsis with ARDS")
md$age <- as.numeric(md$age_years)
md$sex <- factor(md$gender)
md$batch <- factor(md$batch_2)
d <- fread(file.path(input_dir, "GSE66890_expression.tsv.gz"), data.table = FALSE, check.names = FALSE)
map <- read.csv(file.path(input_dir, "GPL6244_unambiguous_mapping.csv"))
symbols <- map$gene[match(d[[1]], map$probe)]
valid <- !is.na(symbols) & nzchar(symbols)
mat <- as.matrix(d[valid, md$sample_id, drop = FALSE])
storage.mode(mat) <- "numeric"
if (as.numeric(quantile(mat, .99, na.rm = TRUE)) > 100) mat <- log2(pmax(mat, 0) + 1)
dt <- data.table(gene = symbols[valid], mat, check.names = FALSE)
collapsed <- dt[, lapply(.SD, median, na.rm = TRUE), by = gene]
e <- as.matrix(collapsed[, -1])
rownames(e) <- collapsed$gene
stopifnot(nrow(md) == 57L, sum(md$ARDS) == 29L, identical(colnames(e), md$sample_id))
saveRDS(list(GSE66890 = list(expr = e, metadata = md)), file.path(input_dir, "public_context_processed.rds"))
