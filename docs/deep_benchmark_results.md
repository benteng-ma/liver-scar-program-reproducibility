# Deep transportability benchmark results

Date: 2026-08-31

All results below are exploratory and preserve the frozen primary classifications.

## A. Positive-negative control ladder

- Non-circular author-label identity cohorts evaluated: 2; cohorts passing the frozen dual-method accuracy/AUROC threshold: 1/2.
- Across the non-circular identity rows, median top-score accuracy was 0.833 and median macro AUROC was 0.856.
- Canonical disease-response cohort-module comparisons positive with both scores after FDR control: 0/12.
- Matched random modules, direction randomization and donor-label permutations remain negative-control references; a positive identity result demonstrates technical sensitivity, not disease-program transportability.

## B. Quantitative failure decomposition

- Leave-one-cohort-out prediction of held-out program ordering: Spearman rho 0.364; predictive R-squared 0.024 across 152 predictions.
- Descriptor contrasts are descriptive because assay, endpoint, annotation and etiology are partly confounded. Bootstrap intervals are provided in `results/deep_benchmark/transfer_failure_decomposition_summary.csv`.
- After fully excluding Watson from both training and testing, leave-one-cohort-out Spearman rho was 0.419 and predictive R-squared was 0.010 across 90 predictions.
- In program pairs with at least 80% coverage after excluding Watson, sign agreement was 0.689 and median program-rank rho was 0.475.
- Triggered Watson audit: the standardized-mean matched-lineage margin changed by -2.045 in fibrotic versus healthy donors (exact P=0.100); disease and technical context are inseparable in this 3+3 cohort.
- Largest prespecified descriptor association by absolute estimate: includes_watson with absolute_effect_difference (5.532).

## C. Frozen discovery and held-out minimal cores

- macrophage_monocyte: 0 genes passed the frozen discovery gates; fewer than five were available, so no core was declared.
- endothelial: frozen Core5 = TFF3, PLPP1, FTL, CPE, FTH1; it did not pass the held-out directional/random-specific rule.
- mesenchymal_hsc_myofibroblast: frozen Core5 = MDK, CST3, TM4SF1, SERPINF1, TMSB10; it did not pass the held-out directional/random-specific rule.

## D. Program transportability report card

The report card is a transparent prioritization resource, not a clinical score. The five highest totals were:

- RAM2019_MAC_SIG_B_SAM: 64.6/100 (measurement 18.0, method 20.0, direction 13.3, random specificity 13.3, endpoint evidence 0.0).
- RAM2019_MAC_SIG_A_SAM: 63.2/100 (measurement 16.5, method 20.0, direction 13.3, random specificity 13.3, endpoint evidence 0.0).
- RAM2019_MAC_SIG_F_CDC1: 57.5/100 (measurement 17.5, method 20.0, direction 13.3, random specificity 6.7, endpoint evidence 0.0).
- RAM2019_MES_MESOTHELIAL: 56.9/100 (measurement 16.9, method 20.0, direction 15.0, random specificity 5.0, endpoint evidence 0.0).
- RAM2019_ENDO_6_SAENDO1: 56.0/100 (measurement 20.0, method 20.0, direction 16.0, random specificity 0.0, endpoint evidence 0.0).

## Interpretation boundary

These analyses can convert the project from a simple null report into a calibrated transportability benchmark only to the extent that technical positive controls pass and failure persists in held-out prediction. They do not establish a causal mechanism, universal scar state, diagnostic biomarker or treatment target.
