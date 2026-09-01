# Phase 4 biological-structure analysis: results and manuscript interpretation

Date: 2026-08-31  
Status: completed post-lock exploratory/validation extension  
Frozen plan: `reports/phase4_biological_structure_plan_2026-08-31.md`  
Run log: `results/logs/phase4_biological_structure_run.json`

## Purpose

Phase 4 was added after the formal Phase 2/3 benchmark had been locked. Its purpose was not to change any cohort classification or rescue failed transportability claims. It asked whether the most reproducible signals contain interpretable biological structure at the gene, multicellular-network, etiology, and composition levels, and whether the overall benchmark geometry quantitatively supports an endpoint-conditioned rather than universal-scar interpretation.

All decision rules, 10,000-permutation/random-module counts, 10,000 donor bootstraps, and seed 20260831 were fixed before inspecting Phase 4 outputs. Formal Phase 2/3 labels remain unchanged.

## 1. Endothelial member-gene coherence

The two endothelial programs selected a priori from the preceding benchmark, Endo2 and SAEndo1, showed gene-level coherence across three independent metabolic-fibrosis/cirrhosis cohorts.

- Endo2: 40 of 45 published genes were eligible in the shared universe; 17/40 (42.5%) had positive Hedges g in all three cohorts. Seven genes were supported by within-program fixed-effect meta-analysis at FDR < 0.05: **TFF3, TSPAN5, PPDPF, EFEMP1, NTS, ADIRF, and LGALS3**.
- SAEndo1: 11 of 12 published genes were eligible; 6/11 (54.5%) were positive in all three cohorts. Five were meta-supported: **GSN, RBP7, PLPP1, PLVAP, and VWA1**.
- Against 10,000 detection-matched random gene modules, Endo2 had a null median of 4 coherent genes (95th percentile 8; 99th percentile 10) versus 17 observed, empirical P = 0.0001. SAEndo1 had a null median of 1 (95th percentile 3; 99th percentile 4) versus 6 observed, empirical P = 0.0013.
- In leave-one-cohort-out analyses, held-out positive-sign retention was 70.8%, 85.0%, and 89.5% for Endo2 and 75.0%, 100%, and 100% for SAEndo1.

Interpretation: the positive program behavior is not driven by a single dominant marker or a generic detection advantage. It reflects a reproducible endothelial leading edge. This remains a transcriptional association and does not establish an endothelial cell state as causal.

## 2. Replicated cross-lineage coupling

Donor-level endothelial, macrophage/monocyte, and mesenchymal program scores were residualized for disease group and prespecified QC covariates in GSE244832 and GSE256398. A pair was called stable only if all four cohort-by-score correlations had the same sign, at least three of four had |Spearman rho| >= 0.30, and both score-specific meta-analyses had FDR < 0.05.

Sixteen cross-lineage pairs met this rule:

- 7 macrophage-mesenchymal edges;
- 6 endothelial-mesenchymal edges;
- 3 endothelial-macrophage edges.

The most robust representative pairs were:

- **SAEndo2-SAM-B:** minimum |rho| = 0.500 across four estimates; meta rho = 0.604 for rank scoring (FDR = 0.00086) and 0.514 for z-mean scoring (FDR = 0.0104).
- **SAEndo1-SAM-A:** minimum |rho| = 0.456; meta rho = 0.528 (FDR = 0.0051) and 0.493 (FDR = 0.0119).
- **TMo-E-Mesothelial:** minimum |rho| = 0.439; meta rho = 0.519 (FDR = 0.0060) and 0.497 (FDR = 0.0119).

SAEndo2 was the highest-degree hub (degree 5); SAM-B, TMo-E, and SAMes-B each had degree 4.

Interpretation: the transferable endothelial signal is embedded in a coordinated vascular-immune-stromal transcriptional structure rather than behaving as an isolated marker module. Residual correlation does not prove ligand-receptor signaling, cellular contact, or causal direction.

## 3. Shared cirrhosis versus etiology divergence

Within GSE256398, MASH cirrhosis-versus-healthy and alcohol cirrhosis-versus-healthy effects were decomposed into a shared cirrhosis component and an etiology-divergence component with 10,000 donor bootstraps.

- Eight programs formed a shared directional backbone across both scoring methods: **Endo2, Endo5, SAEndo1, SAM-A, HSC, VSMC, SAMes-A, and SAMes-B**.
- Only **Endo2** and **SAM-A** also exceeded detection-matched random-module thresholds in both etiologic contrasts and both score methods.
- Endo2 had median shared component 2.452 versus median absolute divergence 0.248.
- SAM-A had median shared component 0.972 versus median absolute divergence 0.422.
- Mesothelial was the sole program classified as etiology-divergent under the frozen rule.

Interpretation: Endo2 and SAM-A define the strongest shared MASH-alcohol cirrhosis backbone, while other programs contribute directional but less random-specific structure. This creates a biologically useful separation between common advanced-disease architecture and etiology-sensitive remodeling.

## 4. Composition and technical-depth sensitivity

Program effects in two metabolic cohorts were refit after adjustment for target-lineage cell count and library depth. Forty-one of 76 context-program pairs retained direction and at least 70% of the unadjusted magnitude across both scoring methods. Five also retained positive HC3 intervals for both adjusted score models:

- GSE244832 MASH F2-F4 versus normal: SAEndo1;
- GSE256398 alcohol cirrhosis versus healthy: Endo2 and SAEndo1;
- GSE256398 MASH cirrhosis versus healthy: SAEndo1 and SAM-A.

Interpretation: these selected signals are not explained solely by broad lineage abundance or library depth. The analysis does not establish cell-intrinsic regulation and cannot remove all compositional or sampling confounding.

## 5. Program-versus-context variance and topology

Across the complete benchmark effect matrix, context accounted for 46.4%-67.1% of variance across lineage-score strata (median 55.2%), compared with 0.1%-15.1% for program identity (median 4.7%); the median interaction/residual fraction was 38.0%.

Permutation regressions of cross-context rank concordance on same endpoint, etiology, assay, and annotation status yielded no FDR-significant coefficient (0/24). Some endpoint and annotation coefficients were nominally positive, but the recorded descriptors were correlated and did not individually explain the topology after correction.

Interpretation: the negative transportability result is now quantitatively reframed. It is not an empty failure to reproduce; it identifies **context dominance** as a measurable property of the scar-program benchmark. The available descriptors are insufficient to assign that context dependence to a single technical or biological cause.

## Revised central claim

Human cirrhosis scar programs are not universally portable as intact disease modules. Instead, transferability is hierarchically organized: a coherent endothelial leading edge, particularly Endo2 and SAEndo1, recurs across metabolically aligned fibrotic/cirrhotic cohorts; Endo2 and SAM-A define a shared MASH-alcohol cirrhosis backbone; and this signal is embedded in a replicated vascular-immune-stromal coupling network. At the same time, context explains substantially more benchmark variance than program identity, setting an explicit boundary on universal marker use.

## Clinical and basic-science relevance

- Clinical: the work identifies which components are plausible candidates for cross-cohort biomarker development (Endo2/SAEndo1 member genes) and which require indication-specific calibration. It argues against transferring a single scar score across endpoints without revalidation.
- Basic science: the work proposes a testable multicellular architecture linking scar-associated endothelial, macrophage, and mesenchymal programs, while separating shared advanced-cirrhosis structure from etiology divergence.
- Translational boundary: none of these analyses establishes diagnostic accuracy, prognosis, treatment response, or causal signaling. Prospective tissue and spatial/perturbational validation remain required.

## Figure allocation

- Main Figure 8: coherent member-gene leading edge, matched-random calibration, replicated cross-lineage network, etiology geometry, composition adjustment, and variance partition.
- Supplementary Figure S7: full endothelial member-gene audit, meta-supported core, leave-one-cohort-out retention, and random-module quantiles.
- Supplementary Figure S8: complete cross-lineage meta-maps, score-method concordance, and network degree.
- Supplementary Figure S9: full etiology decomposition, bootstrap intervals, composition sensitivity, and negative descriptor-regression boundary.

## Auditable outputs

All numerical results are in `results/phase4/`; figure source-data files are in `results/source_data/`; finalized PNG/PDF figures are in `results/figures/`; the reproducible implementation is `scripts/phase4_biological_structure.py` plus `scripts/make_phase4_figures.py`.
