# Cross-cohort synthesis and frozen classifications

Date: 2026-08-30

This synthesis preserves endpoint, annotation, assay, and evidence-tier boundaries. It does not combine donor-level expression matrices across studies and does not pool clinical cirrhosis, advanced fibrosis, mixed-stage MASH, and mixed F2-F4 fibrosis into one disease label.

## Evidence available

Five independent validation resources now contribute donor-level effects:

- GSE202379: author-label snRNA-seq and the only cohort with eligible formal-primary rows.
- GSE290642: reconstructed-label scRNA-seq sensitivity; only endothelial passes the donor-state gate.
- GSE244832: author cluster mapping and MASH-versus-normal sensitivity; donor-level fibrosis stage is unavailable.
- GSE210077 Watson six-donor subset: author-label snRNA-seq mixed F2/F3/F4 sensitivity.
- GSE181483: reconstructed-label scRNA-seq 2+2 directional display; no standalone p values or intervals.

The harmonized matrix contains 274 cohort/contrast/program/score rows. Endpoint families and analysis tiers remain explicit in `results/meta/cross_cohort_effect_matrix.csv`.

## Comparable advanced-fibrosis sensitivity meta-analysis

The only cross-study endpoint pair close enough for effect-only synthesis is endothelial F3-F4 versus F0 in GSE202379 and F4 versus F0 in GSE290642. This is a sensitivity meta-analysis because GSE290642 uses reconstructed broad labels. Six programs are measurable in both cohorts; SAEndo2 is unavailable in GSE202379 and is excluded.

ENDO2 has positive fixed-effect intervals for both scoring methods. Its singscore REML interval remains positive, but its standardized-mean REML interval crosses zero and the fixed standardized-mean analysis assigns 72.5% weight to one cohort, exceeding the frozen 70% rule. ENDO4 and SAEndo1 have positive fixed intervals for standardized mean only; their REML intervals cross zero. No program has positive REML intervals with both scores while also satisfying the evidence-tier and weight requirements.

Hartung-Knapp and leave-one-study-out analyses are not applied because only two cohorts are available. The primary meta-analysis is not run because no endpoint has two comparable author-label, formal-primary cohorts.

## Frozen classification result

Eighteen programs show a positive Hedges g with both score methods in at least one independent cohort and therefore meet only the weak `INDEPENDENT_DIRECTIONAL_SUPPORT` definition. `RAM2019_MAC_SIG_D_TMO` has no such support. This directional label does not require a confidence interval or cross-cohort meta-analysis and must not be read as replication.

Final higher-order counts are:

- `WITHIN_CELL_STATE_REPLICATED`: 0/19.
- `PAN_CIRRHOTIC_TRANSPORTABLE`: 0/19.
- `ASSAY_ROBUST`: 0/19; all assay-transfer labels remain `UNRESOLVED`.
- Spatial validation: 0/19 eligible because of panel coverage, not biological absence.

## Interpretation boundary

The benchmark currently supports cohort-specific and endpoint-specific directional signals, most notably selected macrophage programs in the GSE244832 MASH sensitivity analysis. It does not support a universal cirrhosis scar program, formal cross-etiology transfer, cross-assay robustness, diagnostic performance, or causal mechanism. The null higher-order classifications are an evidence-resolution result: the available external cohorts do not provide two comparable, author-labelled, adequately gated formal-primary tests.
