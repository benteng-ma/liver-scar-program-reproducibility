# Phase 3 manuscript-enrichment amendment

Frozen: 2026-08-31, before inspection of any GSE256398 program score or any new Phase 3 numerical outcome.

## Purpose

Increase biological interpretation, external evidence, and audit depth without changing the Phase 1/2 program inventory, primary thresholds, or the already-frozen higher-order classifications. Phase 3 is a post-lock external-enrichment analysis. A new result may be reported as independent directional or endpoint-specific support, but it cannot silently redefine an earlier success criterion.

## A. Sixth-cohort external validation: GSE256398

- Include only the 26 human donor libraries deposited as separate GEO samples; exclude four mouse libraries.
- Freeze GEO sample titles as the group authority: healthy control (n=6), alcohol-associated cirrhosis (n=4), alcohol-associated hepatitis (n=5), MASLD F0 (n=3), MASH fibrosis F2-3/F3 (n=4), and MASH cirrhosis (n=4).
- Treat every GEO sample as one donor. Never pool cells across donors before inference.
- Apply the existing QC envelope: 200-6,000 detected genes and at most 20% mitochondrial counts. Retain a donor-lineage only when it contains at least 30 cells; 20 cells is sensitivity only.
- Reconstruct broad endothelial, macrophage/monocyte, and mesenchymal/HSC lineages without using disease labels, fibrosis labels, program genes, or program scores. Use the frozen identity-marker panel and the existing two-anchor plus score-margin logic. Author fine-state labels, if recoverable from a publication object, may be used only after source verification.
- Aggregate raw counts to donor-lineage pseudobulks and normalize as log2(CPM+1). Use the same frozen singscore and standardized-mean definitions, primary program-coverage gate, HC3 intervals, donor-label permutations, and 1,000 expression/detection-matched random modules.
- Fixed contrasts, in this order:
  1. MASH cirrhosis versus healthy control;
  2. alcohol-associated cirrhosis versus healthy control;
  3. MASH fibrosis F2-3/F3 versus MASLD F0;
  4. MASH cirrhosis versus MASLD F0;
  5. MASH cirrhosis versus alcohol-associated cirrhosis, as an etiology contrast rather than a disease-presence contrast;
  6. ordinal MASLD F0 -> MASH fibrosis -> MASH cirrhosis.
- The within-dataset MASH and alcohol contrasts will be used to separate endpoint alignment from assay and laboratory context. No result will be described as pan-etiologic unless it satisfies the already-frozen Phase 2 definition.

## B. Clinical-histology specificity in GSE202379

- Analyze donor-level program scores against fibrosis, steatosis, ballooning, and inflammation as separate ordinal axes.
- Report Spearman effects with 10,000 donor-label permutations and Benjamini-Hochberg correction within each lineage, score method, and histology family.
- Fit two descriptive HC3 models: a minimal age/sex-adjusted model and a joint histology model containing all four axes plus age and sex. Report standardized slopes, incremental R-squared, variance-inflation factors, and donor counts.
- Histology variables are correlated and can be mediators; adjusted associations are specificity diagnostics, not causal effects.

## C. Metabolic progression in GSE244832

- Use all 18 donors in the fixed order NORMAL=0, MASL=1, MASH=2.
- Report Spearman trend, Theil-Sen slope, 10,000 label permutations, and FDR within lineage and score method.
- Retain the existing MASH-versus-NORMAL analysis unchanged. The group-level F2-F4 label remains insufficient for a fibrosis-stage claim.

## D. Stability and redundancy analyses

- Quantify within-lineage program overlap using Jaccard, overlap coefficient, shared-gene counts, and a fixed hierarchical clustering of 1-Jaccard distance.
- Estimate the effective number of program tests from the donor-score correlation eigenvalues; this is descriptive and will not replace frozen FDR.
- Bootstrap donors within comparison groups 10,000 times to estimate effect-sign probability, rank intervals, top-five probability, and pairwise ordering stability.
- Perturb report-card weights with 100,000 symmetric Dirichlet draws. Report rank acceptability and top-five probability; no weight draw creates a clinical score or replication label.
- Repeat classification-relevant summaries over the fixed sensitivity grid: 20/30-cell gates, 60/80/90% coverage, and 90/95/99th matched-random percentiles. Frozen primary results remain the 30-cell, 80%-coverage, 95th-percentile specification.

## E. Pathway-level transportability

- Compute donor-level gene-wise standardized effects only for contrasts passing the donor-lineage gate.
- Use a versioned public Reactome human gene-set snapshot fixed before inspecting enrichment results. Perform preranked enrichment separately by cohort, lineage, and contrast.
- Compare pathway direction/rank transfer with the published-program transfer results. Pathway convergence cannot rescue or relabel a failed published program.
- Do not use pathway results to select new Core5 genes.

## F. Precision and future-study design

- Report robust interval width, minimum detectable standardized effect at 80% power, and achieved donor allocation for every eligible contrast.
- For SAM-A, SAM-B, TMo-E, and the endpoint-aligned endothelial programs, simulate prespecified balanced donor counts using observed variance with bootstrap uncertainty. Provide planning curves for direction confirmation and effect-size calibration.
- These are prospective design estimates, not evidence that the programs are clinically validated.

## Locked boundaries

- The Ramachandran program inventory remains 19 programs and 2,893 program-gene rows.
- GSE136103 remains discovery/internal recovery only.
- No integrated UMAP, cell-level disease p value, diagnostic classifier, docking, causal network, or therapeutic claim will be added.
- Phase 2 formal labels remain frozen unless a separately stated, previously defined rule is genuinely met by the new independent cohort. Any such event must be reported alongside the post-lock status of GSE256398.
- All Phase 3 outputs must be versioned, checksum-audited, tested, and explicitly separated from formal Phase 2 results.
