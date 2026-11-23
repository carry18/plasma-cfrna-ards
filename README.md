# Plasma cfRNA analysis across ARDS clinical states

Analysis and figure-generation code for **Plasma cell-free RNA identifies a candidate platelet-related host-response program across ARDS clinical states**.

The workflow starts from the final raw count matrix and frozen sample/clinical metadata: 113 independent plasma profiles, comprising 29 MV ARDS, 17 nonintubated ARDS, 33 pneumonia without ARDS and 34 healthy controls. It performs differential-expression and pathway analysis, co-expression network discovery/projection, clinical association analysis, independent whole-blood evaluation, stability reconstruction and figure/table export.

## Data access

Internal sequencing data are deposited in **GSA-Human under accession HRA022356** and are available through controlled access, subject to approval by the designated Data Access Committee. The cohort count matrix, participant-level clinical records and private identifier mappings are not included here. Authorized users must obtain and prepare the corresponding frozen inputs described in docs/input_schema.md. Access to raw sequencing reads alone does not automatically provide the clinical metadata needed to reproduce all clinical panels.

The independent external dataset **GSE66890** is public. Its source matrix and platform annotation used here are included in data/external. See THIRD_PARTY_NOTICES.md for reference-data sources.

## Environment

Use R 4.4.x with the package versions recorded in scripts/environment_versions.csv. Required packages include edgeR, limma, fgsea, WGCNA, BiocParallel, statmod, data.table, ggplot2, patchwork, svglite, ragg, scales, VennDiagram, magick, jsonlite and digest. Python 3 scripts use the standard library. Figure export uses Arial and R graphics devices; font and device versions can affect exact image rendering.

The scripts do not install dependencies automatically. Configure the external libraries needed by the R graphics packages for your platform. `CF_PYTHON` may point to a specific Python executable.

## Inputs and configuration

Place the authorized matrix and metadata in data/private, following docs/input_schema.md. Do not replace them with normalized counts, a different clinical time window or a differently ordered matrix.

Default file names are defined in scripts/config.R. They can be overridden with `CF_COUNTS`, `CF_CLINICAL`, `CF_TABLE1_CLINICAL`, `CF_GEO_MATRIX` and `CF_GEO_PLATFORM`. When replacing the clinical file, set both CF_CLINICAL and CF_TABLE1_CLINICAL.

## Run order

Run from the repository root, in separate steps:

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

Figure 1 can run after 01, Figure 2/S1 after 02, Figure 3/4/5 after the required network and clinical results, and Figure 6 after 05. FigureS2 requires all 350 successful 06 checkpoints and exports resampling tables before the final table-assembly step. Optional 07_technical_qc.R records that upstream technical-replicate QC cannot be repeated using the final one-profile-per-participant inputs.

Outputs are written to results/; inputs and outputs are not committed to the repository. Network resampling is computationally intensive and uses single-threaded WGCNA per draw. Production checkpoints are reusable only when the algorithm/input fingerprint agrees; never copy old checkpoints into a new input run.

## Analysis overview

Differential-expression models include clinical group and experimental batch. Hallmark enrichment uses Benjamini–Hochberg correction across 50 gene sets within each contrast. Clinical associations are evaluated in nonintubated ARDS and pneumonia without ARDS, with correction across 30 tests. Network reconstruction and resampling assess module stability, and GSE66890 provides an independent whole-blood evaluation.

Each figure script exports the figure and its underlying source data. The table script assembles cohort characteristics and supplementary statistical results from the analysis outputs.

## Version

Prepared 24 November 2026 from the V27 workflow. Source checksums are listed in CODE_SHA256.json. Code repository: https://github.com/carry18/plasma-cfrna-ards . An archived release DOI has not yet been assigned.
