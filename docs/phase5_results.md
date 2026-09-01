# Phase 5 results and frozen-gate decision

Completed: 2026-09-01

Phase 5 was executed under the outcome-blind plan frozen at 2026-09-01 11:52:16 +08:00. It is a post-lock explanatory and orthogonal-validation layer. It does not modify the Phase 2 definitions or the frozen `0/19` labels for within-cell-state replication, pan-cirrhotic transportability, or assay robustness.

## Executive result

Phase 5 sharpened the positive endothelial finding without converting it into a universal program. In the two metabolic cohorts with positive broad-lineage effects, Endo2 was driven predominantly by increased expression within pre-existing endothelial states rather than by shifts in state abundance. The signal also localized to histological scar regions in an independent eight-explant spatial study. However, the positive fine states were not homologous under the frozen non-program-marker matching rule, the primary gene sets did not recur under the frozen rules in the second spatial study, and neither program separated alcohol-associated cirrhosis from alcohol-associated hepatitis with both scoring methods. Accordingly, no primary program met the state, context-specificity, or spatial-replication gate.

## Module 5A: fine-state abundance versus intensity

The analysis generated 12,084 donor-state program scores, 73 state-abundance effects, 368 state-intensity effects, 112 exact Kitagawa decompositions, and 343 cross-cohort state matches. Seven state matches passed the marker-correlation and top-marker-overlap rule, but none connected qualifying positive states for either primary program.

- Endo2 had four qualifying positive states across GSE244832 and GSE256398. SAEndo1 had five across the three cohorts. Because no qualifying state pair met the frozen cross-cohort marker-match rule, both labels are `NOT_STATE_SUPPORTED_POST_LOCK`.
- In GSE244832 MASH F2-F4, within-state intensity accounted for 93.5% of the Endo2 rank-score difference and 96.0% of its standardized-mean difference. For SAEndo1 the corresponding shares were 92.7% and 90.9%. Both intensity-component bootstrap intervals were positive for both scores.
- In GSE256398 MASH cirrhosis, within-state intensity accounted for 84.1% of the Endo2 rank-score difference and 109.7% of its standardized-mean difference; the latter exceeded 100% because the abundance component opposed the total difference. SAEndo1 showed a mixed decomposition: intensity accounted for 48.9% and 56.8%, with abundance accounting for 51.1% and 43.2%. Both components were positive for both scores.
- GSE202379 did not show a consistent positive decomposition for Endo2, and SAEndo1 remained score-discordant or uncertain.

Interpretation: the metabolic-cohort endothelial leading edge is not explained solely by enrichment of a pre-existing fine state. It includes within-state transcriptional activation. Yet the activated source clusters cannot be claimed as one portable endothelial subtype across studies.

## Module 5B: alcohol-hepatitis context specificity

Neither primary program met `CIRRHOSIS_CONTEXT_ENRICHED_POST_LOCK`.

- Endo2 alcohol-associated cirrhosis versus alcohol-associated hepatitis produced Hedges g=1.05 (HC3 95% CI -0.20 to 2.30) for rank scoring and g=2.34 (1.08 to 3.60) for standardized mean. Both estimates exceeded their matched-random 95th percentiles, but the rank-score interval crossed zero.
- SAEndo1 produced g=0.56 (-0.70 to 1.82) and g=0.58 (-0.74 to 1.90), without dual-score random specificity.
- Alcohol-associated hepatitis versus healthy liver was reported independently and did not establish a fibrosis-exclusive signal.

Interpretation: Endo2 is enriched in the cirrhosis-versus-hepatitis contrast by one effect-scale criterion, but the two-method gate fails. The result should be described as injury-context modulation rather than fibrosis specificity.

## Module 5C: independent spatial validation

### Chung et al. 2022

The public author-processed matrix contained paired fibrotic-region and parenchymal-region averages for eight independent cirrhotic explants. All eight donors had higher scar-region scores for both primary programs and both score methods.

- Endo2 primary set: median scar-minus-parenchyma effect 1.194 for standardized mean and 0.0433 for rank scoring; both exceeded the 10,000 matched-random 95th percentiles (0.988 and 0.0305).
- SAEndo1 primary set: standardized mean exceeded random expectation (1.205 versus 1.148), whereas rank scoring did not (0.0291 versus 0.0406).

Endo2 therefore passed a clearly labelled descriptive donor-region substitute analysis. The strict frozen spatial test was not evaluable because the public file was a donor-by-region average rather than a spot-by-gene matrix.

### Hammond et al. 2025

The public supplements contained author-defined scar-cluster top-100 genes for seven independent samples and spot-level cell2location output for three samples, but no complete spot-by-gene expression matrix with scar labels for all donors.

- The Endo2 primary set produced three top-100 hits in two of seven samples, with recurrence rank weight 1.42 versus a random 95th percentile of 0.97. It failed the frozen descriptive recurrence rule because the number of positive samples did not exceed its random 95th percentile.
- The SAEndo1 primary set produced two hits in one of seven samples and failed.
- The secondary 17-gene Endo2 set produced 19 hits across all seven samples and a recurrence rank weight of 12.53 versus a random 95th percentile of 2.37. This is strong secondary localization, not a primary-gate pass.
- In the two samples with evaluable author scar-cluster labels, the fraction of scar-associated endothelial states EC4-EC6 among EC1-EC6 was higher in scar clusters than other spots by 0.471 and 0.050.

Because neither study satisfied the full spot-level frozen test and neither primary set passed both independent-study rules, both primary programs are `NOT_SPATIALLY_REPLICATED_POST_LOCK`.

## Module 5D: conditional regulatory triangulation

The preregistered trigger required at least one primary program to be both `STATE_SUPPORTED_POST_LOCK` and `SPATIALLY_REPLICATED_POST_LOCK`. Neither program met either combined prerequisite. E-MTAB-13131 snATAC analysis is therefore `NOT_TRIGGERED_BY_FROZEN_RULE`; no ATAC outcome was searched.

## Manuscript-level conclusion

The strongest defensible positive claim is now more precise: an Endo2-centered endothelial transcriptional leading edge is activated within states in aligned metabolic fibrosis and is localized to histological scar, with broader Endo2 genes recurring across a second spatial cohort. The same evidence does not identify a conserved fine endothelial subtype, a fibrosis-exclusive signal, or a primary-set spatial replication. This hierarchy increases biological meaning while preserving the central clinical boundary: the published programs are not transferable staging, diagnostic, prognostic, or treatment-selection scores.

## Auditable outputs

- Primary conclusion matrix: `results/phase5/phase5_primary_conclusion_matrix.csv`.
- Regulatory trigger record: `results/phase5/phase5_regulatory_trigger.json`.
- State outputs: `results/phase5/state_*` and `results/phase5/primary_state_*`.
- Alcohol-context outputs: `results/phase5/alcohol_context_*`.
- Spatial outputs: `results/phase5/spatial_2022_*` and `results/phase5/spatial_2025_*`.
- Main figure: `results/figures/figure_9_state_context_spatial_boundary.{png,pdf}`.
- Supplementary figure: `results/figures/supplementary_figure_10_phase5_audit.{png,pdf}`.
