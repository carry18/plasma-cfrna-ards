# Plasma cfRNA analysis across ARDS clinical states

Analysis code for **Plasma cell-free RNA identifies a candidate platelet-related host-response program across ARDS clinical states**.

The study includes 113 independent plasma profiles: 29 MV ARDS, 17 nonintubated ARDS, 33 pneumonia without ARDS and 34 healthy controls. The workflow covers differential expression, pathway enrichment, co-expression networks, clinical associations, external evaluation and figure/table generation.

## Data

Plasma cfRNA sequencing data are deposited in **GSA-Human, accession HRA022356**, under controlled access through the designated Data Access Committee. Analysis requires the raw count matrix and matching sample/clinical metadata described in [docs/input_schema.md](docs/input_schema.md).

The external whole-blood dataset is [GSE66890](https://www.ncbi.nlm.nih.gov/geo/query/acc.cgi?acc=GSE66890), with platform annotation [GPL6244](https://www.ncbi.nlm.nih.gov/geo/query/acc.cgi?acc=GPL6244). Download the files listed in [data/external/README.md](data/external/README.md) before running the analysis.

## Requirements

- R 4.4.x; package versions are recorded in [scripts/environment_versions.csv](scripts/environment_versions.csv).
- R packages: edgeR, limma, fgsea, WGCNA, BiocParallel, statmod, data.table, ggplot2, patchwork, svglite, ragg, scales, VennDiagram, magick, jsonlite and digest.
- Python 3, using the standard library.
- Arial for figure labels and the system libraries required by the R graphics packages.

## Inputs

Place `count_matrix113.tsv` and `sample_info113.csv` in `data/private/`. See [input_schema.md](docs/input_schema.md) and the [metadata header template](docs/sample_info113_header.csv).

Paths are configured in `scripts/config.R`. Optional overrides are `CF_COUNTS`, `CF_CLINICAL`, `CF_TABLE1_CLINICAL`, `CF_GEO_MATRIX` and `CF_GEO_PLATFORM`. Set both clinical-file overrides when using a different metadata path. `CF_PYTHON` selects the Python executable.

## Run

Run from the repository root:

```bash
Rscript scripts/01_prepare_samples.R
Rscript scripts/02_primary_analysis.R
Rscript scripts/03_network_analysis.R
Rscript scripts/04_clinical_analysis.R
Rscript scripts/05_external_analysis.R
Rscript scripts/06_resampling.R
python3 scripts/figures/Figure1.py
Rscript scripts/figures/Figure2.R
Rscript scripts/figures/Figure3.R
Rscript scripts/figures/Figure4.R
Rscript scripts/figures/Figure5.R
Rscript scripts/figures/Figure6.R
Rscript scripts/figures/FigureS1.R
Rscript scripts/figures/FigureS2.R
Rscript scripts/tables/supplementary_tables.R
python3 scripts/internal/validate.py
```

Outputs are saved to `results/`. Figure scripts export PNG, PDF, SVG and TIFF files together with source data. Step 06 reconstructs networks for 350 resampling draws and saves checkpoints for resuming the same analysis.

## Analysis

Expression models include clinical group and experimental batch. Hallmark enrichment uses Benjamini–Hochberg correction across 50 gene sets within each contrast. Clinical associations are evaluated in nonintubated ARDS and pneumonia without ARDS, with correction across 30 tests. Network reconstruction and resampling assess module stability; GSE66890 provides an independent whole-blood evaluation.

Reference-data sources are listed in [THIRD_PARTY_NOTICES.md](THIRD_PARTY_NOTICES.md).
