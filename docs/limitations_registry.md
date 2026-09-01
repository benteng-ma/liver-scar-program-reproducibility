# Limitations registry

| ID | Limitation | Consequence | Frozen mitigation |
|---|---|---|---|
| L01 | GSE136103 is the source of all 19 current frozen programs | No external claim within GSE136103 | Label internal recovery only; external evidence must come from other accessions |
| L02 | E-MTAB-10553 RDS embeds GSE136103 | Donor/cell reuse can inflate apparent validation | Use only six origin-labelled normal donors; exclude embedded public cells |
| L03 | GSE181483 has 2+2 human donors | Unstable uncertainty and invalid standalone inference | Directional display only; no formal cohort-level p-value claim |
| L04 | GSE202379 contains repeated lobes | Region-level pseudoreplication | Aggregate or nest regions within donor; count 47 donors |
| L05 | GSE290642 is mostly F1/F1-2 and lacks etiology | Cannot be called a cirrhosis or cross-etiology cohort | F4-only cirrhosis contrast; all-fibrosis analysis is separately labelled sensitivity |
| L06 | GSE244832 condition is not on GEO sample pages | Risk of guessing donor groups | Exclude formal effects until processed metadata mapping is ingested and verified |
| L07 | HBV resource is CD45-enriched and raw-only | Does not test endothelial/mesenchymal programs and would require FASTQ | Secondary macrophage-only resource; not required for core benchmark |
| L08 | MERFISH panel has 317 biological genes | Frozen programs have inadequate coverage | `CONDITIONAL_GO_NO_SPATIAL`; no spatial program effect or negative conclusion |
| L09 | Full spatial h5ad was not locally downloaded | Phase 0 checked schema/sample, not full-object load | Full download is unnecessary after spatial exclusion; retain HCA file hashes and metadata |
| L10 | Etiology is incomplete in several public records | Pan-cirrhotic claims could be overbroad | Restrict formal etiology contrasts to explicitly recovered labels; otherwise `UNRESOLVED` |
| L11 | Conda R/Python executables are blocked by host Windows code-integrity policy after installation | R workflow cannot be executed or locked on this host yet | Use verified `.venv` only for Phase 0; require executable isolated R or an explicitly approved equivalent before Phase 2 |
