# Phase 4 post-lock biological-structure plan

Frozen: 2026-08-31, before inspection of any Phase 4 numerical output.

## Purpose

Phase 4 is intended to strengthen the biological narrative without changing any frozen Phase 2 classification or the Phase 3 post-lock labels. It asks why the aligned metabolic comparisons preserve endothelial program ordering and whether that conditional transfer can be localized to reproducible member genes, donor-level cross-lineage coordination, or a shared cirrhosis component rather than technical composition.

All Phase 4 results are exploratory biological-structure analyses. Negative, null, heterogeneous, or unstable findings will be retained. No threshold will be relaxed after inspection, no program will be redefined from a held-out outcome, and no Phase 4 result can establish a causal cell-cell interaction, clinical biomarker, or therapeutic mechanism.

Random seed: `20260831`.

## Analysis A: three-cohort endothelial member-gene coherence

### Frozen cohorts and endpoints

1. GSE244832 endothelial MASH F2-F4 group versus normal sensitivity.
2. GSE256398 endothelial MASH cirrhosis versus healthy post-lock sensitivity.
3. GSE202379 author-labelled endothelial clinical NASH cirrhosis versus healthy.

The preselected programs are ENDO2 (`RAM2019_ENDO_2`) and SAEndo1 (`RAM2019_ENDO_6_SAENDO1`) because Phase 3 identified program-level aligned MASH support before this analysis. No additional program may be promoted based on Phase 4 output.

### Gene universe and estimands

- Use gene-level donor-effect rows already generated for the Reactome analysis.
- Restrict to genes present in all three cohort contrasts with detection fraction at least 0.20 in each.
- For every eligible member gene, report the three Hedges-g values, three-of-three sign agreement, the minimum and median effect, a fixed-effect estimate, a REML random-effects estimate where estimable, 95% intervals, Cochran Q, and I-squared.
- A `coherent_positive_member` requires positive Hedges g in all three cohorts. A `meta_supported_member` additionally requires a positive fixed-effect 95% interval and Benjamini-Hochberg FDR below 0.05 within the tested program. Random-effects results remain the higher-uncertainty sensitivity boundary.

### Frozen controls

- Perform all three leave-one-cohort-out rotations: select positive genes in the two discovery cohorts, then report sign retention in the untouched third cohort.
- Compare each program with 10,000 same-size random modules drawn from the shared gene universe with per-gene matching on average detection-decile. The one-sided empirical P value is `(1 + random modules at least as coherent as the real program) / 10001`.
- Report complete results even if neither program exceeds the random control.

## Analysis B: replicated donor-level cross-lineage coupling

### Frozen cohorts

- GSE244832 and GSE256398, because both contain all three target lineages for each included human donor and donor-level program scores for two scoring methods.

### Estimands and adjustment

- Evaluate every endothelial-macrophage, endothelial-mesenchymal, and macrophage-mesenchymal program pair; no pair is selected in advance from its correlation.
- For each cohort and score method, residualize each donor score on the deposited disease group. If donor-lineage cell count and library size are available, add log10(cell count + 1) and log10(library size + 1) as technical covariates in the primary residualization; disease-group-only residualization is sensitivity.
- Estimate Spearman correlation on the residuals and obtain two-sided P values from 10,000 within-cohort donor-label permutations.
- Combine cohort-specific correlations separately for each score method by Fisher-z inverse-variance weighting.

### Frozen coupling label

A `stable_cross_lineage_coupling` requires the same correlation sign in both cohorts and both score methods, absolute rho at least 0.30 in at least three of the four cohort-method estimates, and Benjamini-Hochberg FDR below 0.05 for the two-cohort meta-correlation with both scoring methods. This is coordinated covariation, not evidence of ligand-receptor causality.

## Analysis C: shared-cirrhosis versus etiology-divergent program geometry

Use only the frozen GSE256398 contrasts MASH cirrhosis versus healthy, alcohol-associated cirrhosis versus healthy, and MASH versus alcohol-associated cirrhosis.

For each program and score method:

- define the shared-cirrhosis component as the mean of the two cirrhosis-versus-healthy Hedges-g values;
- define etiology divergence as half their difference and cross-check its direction against the direct MASH-versus-alcohol contrast;
- use 10,000 disease-group-stratified donor bootstraps to obtain intervals for both components;
- retain the existing matched-random 95th-percentile indicators without re-estimating a more favorable null.

A `shared_directional_backbone` requires positive disease-versus-healthy effects for both etiologies and both score methods, with the shared component larger in magnitude than etiology divergence for both methods. A `shared_random_specific_backbone` additionally requires both disease-versus-healthy effects to exceed their matched-random 95th percentile with both scores. An `etiology_divergent` label requires opposite-signed etiology deviations with a bootstrap interval excluding zero for both score methods. These labels are descriptive and lineage-specific.

## Analysis D: composition independence in the two metabolic cohorts

- In GSE244832 and GSE256398, fit donor-level program-score models with disease group alone and with disease group plus log10 target-lineage cell count and log10 library size.
- For each prespecified disease comparison, report unadjusted and adjusted standardized coefficients, robust HC3 intervals, sign retention, and magnitude-retention ratio.
- A result is `composition_stable` only if its sign is unchanged and the adjusted coefficient retains at least 70% of the unadjusted magnitude in both scoring methods. This analysis diagnoses technical attenuation and cannot prove cell-intrinsic biology.

## Analysis E: program-versus-context variance and transport topology

- Use the complete Phase 3 effect matrix without pooling incompatible clinical endpoints.
- Within each lineage and score method, fit the balanced two-way effect decomposition `Hedges g = grand mean + program effect + context effect + residual` and report non-negative method-of-moments variance components and the proportion attributable to program identity, context, and residual interaction.
- Build the complete context-by-context Spearman rank-correlation matrix within each lineage and score method. Report whether endpoint/etiology alignment, same assay, and same annotation provenance are associated with higher rank concordance using a pair-level regression with 10,000 context-label permutations. Because descriptors are confounded, estimates are localization statistics rather than causal effects.

## Reporting and manuscript gate

- Produce complete machine-readable outputs and a standalone Phase 4 report before manuscript editing.
- A single six-panel main Figure 8 may summarize the analysis only if at least one prespecified positive structure survives its frozen control. Otherwise all Phase 4 results remain supplementary or are reported as a negative boundary.
- Supplementary figures may show the full member-gene, coupling, composition, etiology, and topology audits regardless of outcome.
- Phase 4 cannot modify the frozen 0/19 formal replication, pan-cirrhotic, or assay-robust counts.
