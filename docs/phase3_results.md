# Phase 3 manuscript-enrichment results

Completed: 2026-08-31

## Executive interpretation

Phase 3 materially enriches the manuscript without converting it into a universal-signature paper. The sixth cohort supplies an independent, donor-resolved, same-platform comparison of MASH and alcohol-associated cirrhosis and a metabolic F0-to-fibrosis-to-cirrhosis trajectory. The strongest new positive result is **conditional rank transfer in metabolic disease**, especially for endothelial programs. The strongest new boundary result is that neither stage trend, random-effects synthesis, nor pathway ranking becomes universally transferable.

The paper can now make three linked claims:

1. published programs retain reproducible relative ordering when disease context is aligned;
2. selective programs, especially ENDO2 and SAEndo1, receive post-lock sensitivity support in MASH cirrhosis, but not formal replication;
3. program, pathway, and effect-size behavior remain strongly conditioned by stage, etiology, cohort, and lineage.

## Sixth external cohort: GSE256398

- Official public human snRNA-seq resource with 26 independent donors: six healthy, three MASLD F0, four MASH fibrosis F2-3/F3, four MASH cirrhosis, four alcohol-associated cirrhosis, and five alcohol-associated hepatitis.
- The fixed QC retained 159,784 nuclei. The frozen broad-lineage rule retained 28,895 endothelial, macrophage/monocyte, or mesenchymal/HSC nuclei for donor pseudobulk analysis.
- All 26 donors contributed all three target lineages. All five fixed binary contrasts passed the 30-cell and three-donor-per-group gates in all three lineages.
- All 19 programs met primary feature coverage (range 0.800-1.000). No identity or coverage threshold was relaxed.
- A metadata audit caught and corrected healthy controls that had initially been assigned the metabolic-order value zero. Healthy controls are now excluded from the metabolic-only ordinal trajectory; binary contrasts and cell assignments were unchanged.

### Binary program effects

The five contrasts generated 190 program-method effect rows. Thirty-three exceeded the 95th percentile of 1,000 expression- and detection-matched random modules.

Programs above the matched-random 95th percentile with both scores were:

- alcohol cirrhosis versus healthy: ENDO2, SAEndo1, SAM-A, TMo-E, and VSMC;
- MASH cirrhosis versus healthy: ENDO2, SAM-A, and cDC1-F;
- MASH fibrosis versus MASLD F0: TMo-D and TMo-E;
- MASH cirrhosis versus MASLD F0: SAM-A, TMo-D, and TMo-E;
- MASH versus alcohol cirrhosis: mesothelial.

Requiring expected direction, positive HC3 interval, and matched-random specificity in both scores left five program-context results:

- alcohol cirrhosis versus healthy: ENDO2 and SAEndo1;
- MASH cirrhosis versus healthy: ENDO2 and SAM-A;
- MASH fibrosis versus MASLD F0: TMo-D.

This pattern is not pan-etiologic. It identifies a small set of testable program-context pairs.

### Metabolic ordinal trend

The metabolic-only trajectory contained 3 MASLD F0, 4 MASH fibrosis, and 4 MASH cirrhosis donors. ENDO2 and ENDO3 singscores showed positive exact-permutation FDR-controlled trends, but their standardized-mean results did not pass FDR. No program therefore met a dual-score ordinal trend criterion.

## Cross-cohort synthesis after adding GSE256398

The phase-3 matrix contains 464 effect rows across six independent datasets. Frozen Phase 2 labels remain unchanged.

### Conditional MASH rank transfer

GSE244832 mixed F2-F4 MASH versus normal and GSE256398 MASH cirrhosis versus healthy showed high endothelial program-order concordance:

- singscore rho = 0.857, sign agreement = 1.00;
- standardized mean rho = 1.000, sign agreement = 1.00.

Macrophage/monocyte rank correlations were 0.600 and 0.657; mesenchymal/HSC rank correlations were negative or near zero. The result is therefore lineage-specific rather than a global MASH signature.

Within GSE256398, MASH and alcohol cirrhosis effects retained high rank correlation in endothelial (rho 0.964 for both scores) and macrophage/monocyte programs (rho 0.943 for both scores), although sign agreement was only 0.57-0.83. MASH fibrosis and MASH cirrhosis effects also retained strong within-cohort program ordering, particularly for macrophage and mesenchymal programs.

### Post-lock MASH cirrhosis sensitivity meta-analysis

GSE202379 clinical NASH cirrhosis versus healthy and GSE256398 MASH cirrhosis versus healthy form a post-lock, stage-and-etiology-aligned two-cohort synthesis. ENDO2 and SAEndo1 had fixed-effect intervals above zero with both scores and satisfied the 70% fixed-weight rule. Neither retained a positive random-effects interval with both scores. Because one cohort uses reconstructed broad labels, the analysis was added after the formal freeze, and k=2, the result is `POST_LOCK_SENSITIVITY_SUPPORT`, not formal replication.

## Clinical-histology specificity

GSE202379 contributed 288 univariate histology rows and 576 age/sex or joint-histology HC3 model rows.

- Five univariate rows passed within-family FDR.
- In the non-end-stage primary exploratory set, standardized-mean SAEndo1 associated with ballooning, and standardized-mean mesothelial associated with fibrosis.
- Minimal age/sex-adjusted models retained standardized-mean ENDO2-steatosis and SAEndo1-ballooning associations after FDR.
- No joint-histology model passed the same FDR threshold.

The programs therefore do not behave as a single fibrosis axis. Some reflect steatosis or ballooning, and correlated histology dimensions erode specificity in joint models.

## Independent metabolic progression in GSE244832

The NORMAL-to-MASL-to-MASH analysis produced 38 program-method rows. Five singscore rows passed within-lineage/method FDR: SAM-A, TMo-E, mesothelial, SAMes, and SAMes-B. None had a corresponding FDR-positive standardized-mean trend. This supports biological prioritization but not a robust stage or progression marker.

## Program redundancy and ranking stability

- The largest off-diagonal within-lineage gene-set Jaccard coefficient was only 0.120.
- Despite low direct gene overlap, donor-score correlations reduced the median effective dimensionality to 2.63 of 7 endothelial programs, 2.22 of 6 macrophage programs, and 1.99 of 6 mesenchymal programs.
- Ten-thousand donor bootstraps produced 264 program-level stability rows and 702 pairwise ordering rows. Several context-specific leaders had high positive-direction and top-five probabilities, but rank intervals commonly spanned three to six positions in small donor groups.
- Under 100,000 symmetric Dirichlet report-card weight draws, top-five probabilities were 95.4% for SAM-B, 78.5% for SAM-A, 74.9% for mesothelial, 65.3% for SAEndo1, and 65.2% for cDC1-F. SAM-B is the most weight-robust prioritization lead, not a validated biomarker.

## Threshold sensitivity

At the frozen 30-cell, 80%-coverage, and 95th matched-random specification, 160 program-contexts were evaluable across the six-cohort evidence network. Nineteen were positive and random-specific with both scores; seven additionally had positive HC3 intervals with both scores, representing nine unique programs with at least one dual-positive random-specific context.

At the stricter 99th random percentile, only eight contexts remained dual-positive/random-specific and four retained dual positive intervals. Raising coverage to 90% further reduced counts. The conclusion is not created by one threshold, but the number of apparent successes is threshold-dependent.

## Reactome pathway-level transport

The official Reactome June 2026 snapshot was fixed before enrichment results were inspected. Fourteen eligible cohort-contrast-lineage analyses yielded 302,163 gene-effect rows and 18,688 pathway-enrichment rows using 1,000 preranked permutations per analysis.

- Same-cohort MASH versus alcohol cirrhosis pathway-rank correlations were 0.631 in endothelial, 0.785 in macrophage/monocyte, and 0.649 in mesenchymal/HSC lineages.
- Cross-cohort MASH pathway-rank correlations were much lower: 0.372, 0.247, and 0.173.
- The matched GSE202379 versus GSE256398 cirrhosis comparison was near zero in endothelial and negative in mesenchymal/HSC; GSE202379 macrophage was ineligible.
- Two hundred forty-two lineage-pathway-direction records were FDR-significant in at least two contexts, but 216 were negative recurrences. Only 26 were positive recurrences.
- Positive same-platform, cross-etiology recurrences included endothelial extracellular-matrix organization, laminin interactions, and ECM proteoglycans, plus macrophage phagocytic and second-messenger pathways.

Pathways do not rescue a universal program. They add a biologically interpretable layer showing within-assay etiologic convergence and cross-cohort attenuation.

## Precision and future design

Across 426 eligible program-effect rows, median minimum detectable standardized effects at 80% power were approximately 1.88 in GSE202379, 1.56 in GSE244832, 1.98 in GSE256398, 2.01 in GSE290642, and 2.29 in Watson6. These public cohorts are powered mainly for large effects. Prospective validation should treat direction confirmation, rank stability, and absolute calibration as separate endpoints and should use substantially larger balanced donor groups for moderate effects.

## Figure architecture

- New main Figure 7 contains six panels: sixth-cohort geometry; full five-contrast effect landscape; selected pairwise concordance; two ordinal trajectories; Reactome transfer; and weight sensitivity.
- New Supplementary Figure S1: GSE256398 identity, cell-count, coverage, and donor-gate audits.
- New Supplementary Figure S2: program gene overlap, effective dimensionality, and report-card weight uncertainty.
- New Supplementary Figure S3: histology specificity and two metabolic progression analyses.
- New Supplementary Figure S4: Reactome transfer matrices and recurrent positive pathways.
- New Supplementary Figure S5: threshold grid and prospective precision/power.
- New Supplementary Figure S6: donor-bootstrap direction and rank stability.

## Conclusion boundary

Phase 3 changes the paper from a five-cohort benchmark dominated by failure analysis into a six-cohort benchmark with a clearer conditional-positive biological result. It does **not** justify universal, diagnostic, causal, mechanistic, therapeutic, or pan-etiologic claims. The most defensible innovation is now the donor-level demonstration that **relative ordering is reproducible in aligned metabolic contexts, while effect calibration, stage specificity, and pathway transfer degrade as biological and technical context diverges**.

## Final manuscript and package

- The manuscript was rebuilt as a 38-page, 6,090-word main document with an exactly 275-word structured abstract, 18 references, three tables, and seven multi-panel main figures.
- The 16-page supplement contains the executed methods, Tables S1-S23 index, Tables S1, S2, and S10 in-document, and six supplementary figures; complete machine-readable tables and figure source data are in `Supplementary_Data.zip`.
- The cover letter and author checklist are one page each. Bing Chen is the sole first author; Fen Wang and Xiao-ming Liu are co-corresponding authors; Acknowledgements is `None.` and the conflict statement is `The authors have no conflicts of interest to declare.`
- All 56 Word-rendered pages were visually inspected; the four final DOCX accessibility audits returned zero findings; 49/49 automated tests passed.
- Final archive: `submission/hepatology_communications/v6_2026-08-31/HepComm_INITIAL_SUBMISSION_PACKAGE_v6_2026-08-31.zip`; SHA-256 `037722147303B3126F6C1F4BAA3A80F0593EBDE10BEF36E9E74DB8A45603D2E7`.
