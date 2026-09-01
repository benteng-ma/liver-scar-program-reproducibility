# Phase 2 scoring and effect implementation addendum

Date: 2026-08-30. This document fixes implementation details that were not numerically specified in the Phase 1 plan. It does not alter program membership, endpoints, cell gates, donor gates, lineage labels, coverage thresholds, or the conclusion boundary.

## Donor-level expression scale

For each cohort and harmonized lineage independently, raw donor pseudobulk counts are converted to `log2(CPM + 1)` using the library size of that donor-lineage pseudobulk. Genes remain in the author feature space. Missing program genes are excluded from the score denominator and are never set to zero.

## Scores

All 19 frozen programs contain `UP_IN_STATE` genes. Direction-aware singscore ranks all measured genes within each donor, uses average ranks for ties, and rescales the mean program rank to the theoretical `[-1, 1]` range for a gene set of the observed size. The standardized-mean score z-standardizes each measured gene across the cohort-lineage reference donors that pass the 20-cell sensitivity gate, then averages the signed gene z values. The reference set is outcome-blind and includes all disease strata.

Programs with at least 80% measured genes are primary. Coverage of 60–79% is flagged, 40–59% is sensitivity only, and below 40% is not evaluated. A threshold boundary is inclusive: exactly 80% is primary.

## Effects and intervals

Positive effects mean a higher program score in the case group. Hedges g uses the pooled donor-level standard deviation and the small-sample correction `J = 1 - 3/(4df - 1)`. The reported robust standard error comes from an HC3 donor-level linear model for the raw mean difference and is transformed onto the Hedges-g scale. The primary 95% interval is `g ± 1.96 × HC3 SE`; the conventional standardized-mean-difference interval is retained as a secondary audit column.

Donor-label permutation tests are two-sided. All label allocations are enumerated when there are at most 100,000 combinations; otherwise 10,000 allocations are sampled without interpreting the resulting p value as stronger evidence than the effect and interval. Random seed is 20260830 with deterministic per-program offsets.

## Matched random modules

Each measured program is compared with 1,000 modules of the same measured size. Candidate genes are matched per program gene on cohort-lineage mean-expression decile and donor detection-rate decile, excluding real program genes and avoiding replacement within each random module when possible. If an exact two-decile bin is exhausted, the nearest Manhattan-distance bin is used and the expansion is recorded. The empirical one-sided p value is `(1 + number of random-module g values at least as large as the real g) / 1001`; performance above the random-module 95th percentile is reported separately for both score types.

## GSE202379 endpoint handling

`NASH with cirrhosis` versus `Healthy control` is the clinical cirrhosis endpoint. F3–F4 versus F0 excludes the `end stage` disease-status stratum in the primary comparison. End-stage versus healthy is reported separately; pooling end-stage F4 donors into F3–F4 is sensitivity only. No result that fails the frozen donor-count gate is estimated merely because the underlying object contains cells.
