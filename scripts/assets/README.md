# Reference assets

- `msigdb_2026.1_Hs.rds`: public annotation only; MSigDB 2026.1.Hs Hallmark, Reactome, KEGG Legacy/Medicus, GO biological processes and Travaglini lung cell gene sets. Rebuild using `../internal/prepare_reference.R` and the corresponding raw MSigDB RDS collections. No cohort counts, sample annotations, differential-expression results or learned modules are stored in this asset. MSigDB usage/citation terms remain those of the source database.
- `published_RR_signatures.csv`: externally published 40-gene predicted and 7-gene validated sets from Supplementary Table ST6, Mitsuyama et al., Respiratory Research 27:277 (2026), doi:10.1186/s12931-026-03708-0. Extracted from the original published supplementary workbook; signature membership is fixed, not learned in this study.
- `figure1_defs.svg`: drawing symbols and style definitions used by Figure 1; no statistical results.
