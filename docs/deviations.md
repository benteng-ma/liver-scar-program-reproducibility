# Deviations

## D01 — Environment materialization

The conda specification was created. Two environment-creation attempts ended with transient HTTP failures; a third download completed, but Windows code-integrity policy blocked the conda Python/R executables (`0xC0E90002`, with an additional R dependency status `0xC0000135`). Phase 0 was therefore re-verified in a repository-local Python `.venv` using the host-signed Python 3.13.11 runtime and pinned Phase 0 packages. No statistical result was generated. An executable isolated R environment and `renv.lock`, or an explicitly approved method-equivalent Python implementation, remains mandatory before Phase 2.

## D02 — Spatial smoke-test scope

The 1.322 GB HCA MERFISH h5ad was not downloaded. Instead, the complete publication spatial source-data workbook and gene panel were read, and a fixed 1 MiB HTTP range of the coordinate-bearing HCA cell-properties table was parsed to verify schema, donor, cell label and x/y fields. This is explicitly not recorded as a full h5ad load. Because zero frozen programs passes spatial coverage, the spatial module was removed rather than expanding the download.

## D03 — GSE244832 donor condition mapping

The publication reports 18 donors and group totals, but GEO sample records omit condition. Donor identifiers were retained with `pending` condition rather than guessing. This cohort cannot enter a disease-effect model until its processed cell metadata is read in Phase 2.

## D04 — Authorized method-equivalent Python route

On 2026-08-30 the user explicitly authorized continuation of the full project. Because Windows code-integrity policy blocks the materialized conda R/Python executables, Phase 2 will use the repository-local virtual environment built from the host-signed Python runtime. This is an implementation deviation, not an analytical deviation: the frozen donor unit, endpoints, eligibility gates, scores, effect-size estimands, random seed, matched-module count, meta-analysis rules and conclusion labels remain unchanged. A Phase 2 lockfile and verification report are required before disease effects are unlocked.

## D05 — GSE290642 heterogeneous deposited feature spaces

The first and later sample archives use different gene-feature spaces. A diagnostic union-based run revealed that treating an unmeasured gene as a zero would confound program ranks with feature version. The diagnostic effect files were overwritten before interpretation. Formal GSE290642 scores and matched modules use the outcome-independent intersection of 21,725 symbols present in all 24 human donors. The union count matrix remains only as an auditable raw aggregation layer.

## D06 — Watson processed raw-layer missing features

The official Dryad AnnData raw layer contains genes with missing values in only a subset of donors. Treating these values as zero would conflate an unobserved feature with absence of expression. Formal Watson scores and random modules therefore use the 25,476-gene intersection observed across all six donors. This is a feature-integrity correction and does not alter endpoints, programs, or thresholds.

## D07 — GSE181483 reconstructed broad lineages

No author processed annotation object was deposited with the four human matrices. To complete the preregistered directional display, broad lineages were reconstructed using the already frozen outcome-blind identity-marker, score-margin, and anchor rules. Ambiguous clusters were excluded. Because the cohort is 2+2, no standalone p values or confidence intervals were computed and the result cannot contribute an inferential replication claim.

## D08 — Restricted cross-cohort meta-analysis

The preregistration anticipated cross-cohort meta-analysis when comparable eligible effects were available. Only one pair met a defensible endpoint match: advanced endothelial fibrosis in GSE202379 and GSE290642. The latter uses reconstructed labels and remains sensitivity only. A primary meta-analysis, Hartung-Knapp adjustment, and leave-one-study-out analysis were therefore not run; these omissions are required by the frozen comparability and minimum-cohort rules rather than analytic convenience.
