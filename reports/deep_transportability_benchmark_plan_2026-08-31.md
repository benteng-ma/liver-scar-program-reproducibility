# Frozen deep transportability benchmark amendment

Freeze date: 2026-08-31

## Status and non-negotiable boundary

This is a second, explicitly exploratory amendment requested after the frozen primary and first exploratory analyses were complete. It cannot change the 19 frozen program definitions, the donor as the inferential unit, any primary endpoint or gate, or any prespecified `WITHIN_CELL_STATE_REPLICATED`, `PAN_CIRRHOTIC_TRANSPORTABLE`, or `ASSAY_ROBUST` classification. The purpose is to determine whether the negative higher-order result can be converted into an interpretable benchmark with calibrated positive controls, quantitative failure decomposition, a held-out minimal-core experiment, and a reusable program report card.

No threshold below may be changed after new outcomes are inspected. All new findings will be labelled exploratory. Failure to recover a held-out core is an allowed outcome.

## Analysis A: positive-negative control ladder

### A1. Lineage-identity positive controls

The three target-lineage identity modules are frozen in `config/deep_benchmark_control_programs.csv`. Their genes are copied from the pre-existing, outcome-blind identity panels in `config/gse290642_identity_markers.yaml`, which were used before disease-effect analysis.

For every cohort, identity modules will be scored in every eligible donor-lineage pseudobulk using a standardized mean and rank-based singscore. The target is the author or frozen broad lineage, not disease status. Report:

- module coverage;
- matched-lineage versus off-lineage Hedges g;
- donor-stratified label-permutation P value when donor overlap permits;
- top-score lineage accuracy and macro one-versus-rest AUROC.

GSE202379 and Watson author labels are the non-circular positive-control evidence. GSE244832, GSE290642 and GSE181483 used overlapping markers during cluster reconstruction and are reported only as pipeline/QC controls.

An identity or canonical disease-response control is evaluable in a cohort only when at least five genes and at least 60% of the frozen module are measured. A five-gene minimal core requires at least four measured genes (80%) in a held-out cohort.

The lineage-identity ladder passes its non-circular positive-control check in a cohort only when top-score lineage accuracy is at least 0.80 and macro one-versus-rest AUROC is at least 0.90 with both score methods. This threshold assesses technical sensitivity only and is not a disease-effect threshold.

### A2. Canonical disease-response intermediate controls

Three frozen, direction-aware modules are used: scar-associated macrophage, LSEC capillarization, and HSC activation. They are not proposed as novel programs and may overlap published scar-state genes. They test whether each cohort contains an interpretable canonical disease-response axis at the donor level. Use the representative endpoint already frozen for each cohort, the same 30/20-cell rules, both score methods, and the same effect estimator. No intermediate control is assumed to be positive in every etiology.

### A3. Negative controls

Retain the existing 1,000 expression/detection-matched random modules for the 19 frozen programs. Add 1,000 donor-label permutations for each eligible canonical disease-response module and 1,000 matched random modules for any frozen five-gene minimal core. Report empirical percentile and two-sided permutation P values. A method is considered capable of separating controls only when non-circular identity controls outperform label-shuffled and size-matched random controls.

## Analysis B: quantitative decomposition of transfer failure

Use one representative endpoint per cohort, exactly as frozen in `reports/exploratory_value_rescue_plan.md`. Create one row per lineage, program, score method and cohort pair.

Outcomes:

- sign discordance;
- absolute difference in Hedges g;
- rank discordance after standardizing effects within cohort-lineage-score strata.

Prespecified explanatory descriptors:

- same versus different assay;
- both author-labelled versus any reconstructed mapping;
- comparable advanced-fibrosis endpoint versus endpoint mismatch;
- absolute and minimum program coverage;
- harmonic mean and minimum group donor count;
- Watson-pair indicator;
- lineage and score method.

Report median contrasts with 10,000 program-cluster bootstrap intervals and Spearman associations for continuous descriptors. These are descriptive associations, not causal variance attribution, because cohort, assay, annotation, etiology and endpoint are partially confounded.

Additionally, estimate leave-one-cohort-out prediction of standardized program effects. For each held-out cohort-lineage-score stratum, use the mean standardized effect of the same program in the remaining cohorts of that lineage to predict the held-out ordering. Report pooled out-of-cohort Spearman correlation and predictive R-squared. Non-positive predictive R-squared is interpreted as absence of useful rank transport, not as proof of biological absence.

## Analysis C: frozen discovery and held-out minimal-core experiment

### Discovery cohorts

- GSE202379: non-end-stage donor-level F0-F4 gene association within the 30-cell lineage gate.
- GSE244832: MASH versus normal donor-level gene effect within the 30-cell lineage gate.

### Candidate universe

For each lineage, use the union of genes in the frozen 19 programs. A gene is eligible only if:

1. it occurs in at least two frozen programs of that lineage;
2. it is measured in both discovery cohorts;
3. it is detected in at least 20% of eligible donor-lineage pseudobulks in both discovery cohorts;
4. its GSE202379 Spearman rho with F0-F4 stage is at least 0.10;
5. its GSE244832 MASH-versus-normal Hedges g is at least 0.20.

All frozen program genes are `UP_IN_STATE`; therefore both discovery effects must be positive. Candidate genes are ranked by the mean of their within-cohort percentile ranks, breaking ties by gene symbol. A lineage core is frozen as the top five genes only when at least five candidates pass. Otherwise no core is declared for that lineage.

### Held-out validation cohorts

- GSE290642 human scRNA-seq;
- Watson six-donor snRNA-seq;
- GSE181483 human scRNA-seq.

The five-gene core is scored without modification using both methods. A `HELD_OUT_DIRECTIONAL_CORE` requires both scores to be positive in at least two held-out cohorts and the observed effect to exceed the matched-random 95th percentile with both scores in at least one held-out cohort. This is an exploratory core label and cannot upgrade the original 19 programs.

No leave-one-gene substitution, threshold relaxation, cohort exchange, or etiology pooling is allowed after discovery results are seen.

## Analysis D: transportability report card

Generate one auditable report card for every frozen program. Five equally weighted domains, each scaled 0-20, are frozen:

1. **Measurement:** 10 × median representative-endpoint coverage + 10 × minimum coverage, each coverage clipped to 0-1.
2. **Score-method reliability:** 20 × fraction of evaluable cohort contexts in which singscore and standardized mean have the same effect direction.
3. **Directional transfer:** 20 × fraction of representative cohort effects in the expected positive direction, averaged across the two score methods.
4. **Matched-random specificity:** 20 × fraction of evaluable cohort contexts in which both score methods exceed their matched-random 95th percentile.
5. **Endpoint evidence:** 5 points each for a positive formal-primary interval with both scores, a dual-score FDR-positive F0-F4 trend, eligible advanced-fibrosis meta-analysis with positive interval under both scores, and original-program assay robustness. Frozen failures remain zero.

The total is the sum of the five domains (0-100). It is a prioritization/readiness index, not a clinical score, probability, validated classifier, or new replication definition. Domain values must always be shown beside the total. No readiness categories or cut-points will be inferred from these 19 programs.

## Multiplicity, seeds and reporting

- Base seed: `20260831`.
- Bootstrap/permutation replicates: 10,000 unless computationally prohibitive; any reduction must be documented before the affected result is interpreted.
- New gene-level discovery P values are descriptive and will not be used for candidate selection; effect thresholds and held-out performance are the selection/validation rules.
- Benjamini-Hochberg adjustment is applied to inferential control-module tests within analysis family.
- All numerical source tables, logs, hashes and scripts are saved under the project repository.

## Interpretation rule

The manuscript may be upgraded from a simple negative result only if the completed analysis demonstrates all of the following:

1. non-circular positive controls establish that the pipeline can detect stable biological identity;
2. transfer loss persists after excluding obvious coverage and score-method failure;
3. failure descriptors or held-out prediction localize where portability is lost; and
4. either a held-out minimal core is recovered or the report card provides a transparent, reusable program-level resource.

The work remains non-mechanistic and non-clinical unless independent experimental or outcome evidence is added.
