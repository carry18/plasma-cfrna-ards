# Input schema

## count_matrix113.tsv

Tab-delimited raw, non-negative integer counts; the first column is `gene_id`, followed by 113 sample columns named CFPL001–CFPL113. Input rows may use gene symbols or `ENSG|length|symbol|type`; compound features are summed by gene symbol before population-specific filtering. Inputs are not normalized or filtered in advance. Pure Ensembl-only identifiers require mapping with the matching annotation first.

## sample_info113.csv

See sample_info113_header.csv for a header-only template. No actual participant rows are included.

- `Sample`: unique CFPL sample identifier matching the count columns.
- `group`: A = MV ARDS (29), B = nonintubated ARDS (17), C = pneumonia without ARDS (33), D = healthy controls (34).
- `Batch`: experimental batch, encoded consistently as in the frozen analysis inputs (Sep3 and Dec2). It is metadata, not inferred from sample IDs.
- Clinical variables used in models: `pf_ratio`, `il6_value`, `wbc_value`, `platelet`, `creatinine_value`, `bilirubin_value`, `hco3`; missing values are `NA`, not zero.
- `il6_censored`: TRUE/FALSE (NA for healthy controls); indicates an original inequality-qualified IL-6 value.
- Table 1 fields: `age`, `sex`, `severity` (APACHE II), `icu_days`, `steroids`, `anticoagulants`, `transfusion`, `ECMO`, `death`. The frozen sex coding is 男/女. Treatment and outcome fields use 0/1 where recorded; healthy controls have NA.

Samples are matched by identifier. The count-column order defines analysis order and must be preserved for exact seeded resampling. The controlled-access release should provide the original frozen values and order; reconstructing clinical time points or replacing missing values is outside this workflow.

## Reference and external inputs

`scripts/assets/` contains fixed annotation definitions, published reference-signature membership and Figure 1 drawing definitions. It does not contain learned modules, cohort expression or precomputed statistical results. Source terms apply; see THIRD_PARTY_NOTICES.md.

The external dataset is NCBI GEO GSE66890, with platform GPL6244. Download `GSE66890_series_matrix.txt.gz` and `GPL6244.annot.gz` into `data/external/`, or supply their paths through CF_GEO_MATRIX and CF_GEO_PLATFORM. The external data files are not included in the current repository tree. See data/external/README.md for GEO links.
