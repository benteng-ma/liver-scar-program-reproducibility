# Frozen analysis plan (Phase 1 preregistration)

Freeze date: 2026-08-30. Formal disease effects are prohibited until this file, the Phase 0 decision and the program checksum are committed.

## Question and unit

The primary question is whether published macrophage/monocyte, endothelial and mesenchymal/HSC/myofibroblast programs transfer across independent human liver cohorts and between scRNA-seq and snRNA-seq. The biological and inferential unit is the donor. Cells, nuclei, tissue regions, sorting fractions and captures are never independent replicates.

## Frozen programs

The 19 program IDs in `literature/program_lineage_matrix.csv` are frozen. Gene membership, direction and provenance are exactly those in `literature/program_inventory.csv`. They may not be modified using validation data. GSE136103 is `REUSED/internal recovery only` for every current program. New published programs may enter only through a dated protocol amendment made before their outcomes are inspected and must remain a separate program family.

## Frozen lineages and labels

Primary lineages are `macrophage_monocyte`, `endothelial`, and `mesenchymal_hsc_myofibroblast`. Author annotations are primary; the hierarchy in `metadata/cell_state_mapping.csv` is frozen. Reference mapping is sensitivity only. Program expression, UMAP position or desired replication may not define a cell label.

## Dataset roles

- **Discovery/internal recovery:** GSE136103. No external-validation label.
- **Primary independent snRNA candidates:** GSE202379; GSE244832 after verified condition mapping.
- **Primary independent scRNA candidate:** GSE290642 human NPC, with F4-only cirrhosis analysis separated from all-fibrosis sensitivity.
- **Small directional support:** GSE181483 human (2+2); GSE210077 Watson snRNA (3+3).
- **Healthy assay baselines:** six unique E-MTAB-10553 donors after excluding embedded GSE136103 cells; GSE185477 paired healthy sc/sn donors.
- **Macrophage-only secondary HBV:** PRJNA833766 only if a processed donor-resolved object is recovered without FASTQ; otherwise omitted.
- **Excluded as an undifferentiated cohort:** the evolving full GSE210077 series. Only publication-defined subsets may be used.
- **Spatial:** excluded from formal program validation because 0/19 programs meets the frozen panel threshold.

## Disease endpoints and etiologies

Primary endpoints are (1) histologic/clinical cirrhosis versus healthy and (2) advanced fibrosis F3-F4 versus F0, analyzed separately. F1-F4 versus F0 is sensitivity only. F2/F3 will not be relabelled cirrhosis. GSE202379 end-stage donors remain a separate stratum unless their cirrhosis status is explicitly recoverable. Etiology labels are analyzed only when explicit; missing etiology is `UNRESOLVED`, never inferred from study context.

## Eligibility gates

- Primary pseudobulk: at least 30 cells/nuclei per donor×harmonized state.
- Cell-count sensitivity: at least 20.
- Formal within-cohort effect: at least 3 eligible donors per comparison group; 2+2 cohorts are directional only.
- Primary program coverage: ≥80%; flagged evaluation 60-79%; sensitivity 40-59%; <40% not evaluated.
- No missing program gene is imputed and no score denominator treats a missing gene as zero expression.

## Per-cohort processing

Each dataset is processed independently from author-provided counts and annotations. Counts from technical fractions, captures or regions belonging to one donor are summed or modelled within that donor before inference. A raw-count matrix is aggregated to `donor × harmonized_state`; cell number, detected genes and library size are retained as QC. No cross-study ComBat, scVI expression, integrated expression matrix or cell-level disease test is allowed.

## Scores and effects

Primary scores are direction-aware singscore and a standardized mean of measured signed genes on donor-level normalized pseudobulk. Camera/roast or an equivalent competitive/self-contained gene-set test is sensitivity. Within each cohort/state, report donor points, standardized mean difference, Hedges g, 95% interval, robust standard error and donor-label permutation where exchangeability permits. P values never replace effect sizes.

## Independent support and meta-analysis

`INTERNAL_RECOVERY` means expected direction in GSE136103 and is not validation. `INDEPENDENT_DIRECTIONAL_SUPPORT` requires expected-direction Hedges g in a cohort not used to define the program. `WITHIN_CELL_STATE_REPLICATED` requires expected direction in at least two independent cohorts, an effect-only meta-analysis whose 95% interval excludes zero, and no single cohort carrying >70% of inverse-variance weight. Fixed-effect and REML random-effect estimates are reported; Hartung-Knapp, I², tau² and leave-one-study-out are mandatory when at least three comparable cohorts exist.

`PAN_CIRRHOTIC_TRANSPORTABLE` additionally requires support in at least two explicitly different etiologies, both scRNA and snRNA compatibility on shared measurable genes, empirical performance above the 95th percentile of matched random modules, and evidence that the result is not entirely composition-associated. Otherwise classifications are limited to the prespecified weaker labels.

## Random modules and perturbation controls

For every dataset/state/program, generate 1,000 random modules matched on measured module size, mean expression decile and detection-rate decile, sampling without replacement where possible. Seed is `20260830`. Additional controls are direction randomization, leave-one-gene-out, leave-one-study-out, leave-one-etiology-out, independent-only, scRNA-only, snRNA-only, shared-gene-only, 20-versus-30-cell gates, author labels versus reference mapping, exclusion of CD45 enrichment, exclusion of non-cirrhotic MASH fibrosis, exclusion of cohorts with <3 donors/group, and within-cohort donor-label permutation.

## Composition versus within-state expression

Donor-level cell proportions and donor×state pseudobulk scores are separate outcomes. Proportions are descriptive with intervals or beta-binomial/Dirichlet-multinomial models only when donor counts support them. A program is `COMPOSITION_ASSOCIATED` if its apparent whole-lineage signal is not reproduced within the harmonized state and tracks donor-level state abundance. No composition change is described as intracellular activation.

## Cross-assay transfer

For every program, report full-program and shared-gene effects separately in scRNA and snRNA. Classify as `ASSAY_ROBUST`, `SNRNA_ATTENUATED`, `SCRNA_ATTENUATED`, `GENE_COVERAGE_DEPENDENT`, `ASSAY_DISCORDANT` or `UNRESOLVED` from effect direction and interval overlap, not from one significant/one nonsignificant comparisons.

## Spatial rule

Formal spatial eligibility was frozen at at least five detected program genes and ≥20% program coverage. No program passed. Therefore no spatial disease effect, location test or negative program call will be run. The MERFISH audit is reported only as a platform-coverage limitation.

## Conclusion boundary

No result will be called causal, mechanistic, diagnostic, therapeutic, universal or a novel cell type. Failure to transfer is reportable and will not trigger module redefinition, threshold lowering, etiology pooling, additional machine learning, WGCNA, pseudotime or cell communication analyses.

