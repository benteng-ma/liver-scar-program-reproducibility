# Frozen Phase 5 analysis plan

Frozen: 2026-09-01 11:52:16 +08:00

Status at freeze: no Phase 5 outcome had been computed or inspected. The current v12 archive is preserved at SHA-256 `B8E050B1ACB007428B20FE0EEDD2F2422EEAA6EB441832A34F0694B234CED0B1`.

## Scope and immutable boundary

Phase 5 is a post-lock biological-explanation and orthogonal-validation extension. It cannot change the frozen 19 programs, the Phase 2 donor gates, the two score definitions, the matched-random principle, or any formal Phase 2 classification. In particular, `WITHIN_CELL_STATE_REPLICATED`, `PAN_CIRRHOTIC_TRANSPORTABLE`, and `ASSAY_ROBUST` remain 0/19 regardless of Phase 5 outcomes.

The two primary endothelial programs are `RAM2019_ENDO_2` (Endo2) and `RAM2019_ENDO_6_SAENDO1` (SAEndo1). The primary spatial gene sets are the already frozen within-program meta-supported members:

- Endo2: `TFF3, TSPAN5, PPDPF, EFEMP1, NTS, ADIRF, LGALS3`.
- SAEndo1: `GSN, RBP7, PLPP1, PLVAP, VWA1`.

The secondary spatial gene sets are the previously frozen coherent-positive members:

- Endo2: `TFF3, TSPAN5, PPDPF, EFEMP1, NTS, ADIRF, LGALS3, LAPTM5, TMSB10, S100A6, VIM, S100A10, CALD1, ANXA2, GUK1, C4ORF48, SNCG`.
- SAEndo1: `GSN, RBP7, PLPP1, PLVAP, VWA1, TMEM88`.

## Module 5A: fine-state abundance-versus-intensity decomposition

### Cohorts and comparisons

- GSE202379: author-labelled clinical cirrhosis versus healthy liver.
- GSE244832: MASH versus normal liver; retained as mixed-stage metabolic disease, not a cirrhosis endpoint.
- GSE256398: MASH cirrhosis versus healthy liver.

Endothelial is primary. Macrophage/monocyte and mesenchymal/HSC/myofibroblast results for their lineage-matched programs are secondary.

### State definition and anti-circularity

1. Use source-generated clusters already present before Phase 5: GSE202379 `SCT_snn_res.0.8`, GSE244832 `seurat_clusters`, and GSE256398 `cluster`.
2. Restrict clusters to cells within the previously harmonized parent lineage. A source cluster split across lineages is treated as separate lineage-specific states.
3. Define state marker profiles and cross-cohort state matches after excluding every gene in the frozen 19-program inventory. Disease labels and Phase 5 program scores are not used to create or merge states.
4. Match states across cohorts using Spearman correlation of non-program marker effects and top-50-marker overlap. A match is considered supported when correlation is positive, marker Jaccard is at least 0.10, and at least three shared non-program markers are present. All pairwise matches are reported.

### Gates and estimands

- Primary donor-state gate: at least 30 cells; sensitivity gate: at least 20 cells.
- A state-specific disease effect requires at least three eligible donors per comparison group.
- Program intensity is calculated from donor-state raw-count pseudobulk using the existing standardized-mean and singscore methods.
- State abundance is the state cell fraction within the parent lineage for each donor. Primary abundance inference uses an arcsine-square-root transformed fraction with HC3 robust intervals; untransformed fractions and donor bootstrap intervals are reported for interpretability.
- The broad-lineage mean difference is decomposed exactly with the symmetric Kitagawa identity:
  - abundance component: sum over states of `(p_case - p_control) * (mu_case + mu_control) / 2`;
  - within-state intensity component: sum over states of `(mu_case - mu_control) * (p_case + p_control) / 2`.
- Ten-thousand donor-stratified bootstrap resamples estimate component intervals and the component share. States failing the donor gate remain in abundance denominators but do not receive state-specific intensity inference.

### Primary interpretation rule

A primary endothelial program is called `STATE_SUPPORTED_POST_LOCK` only when both score methods show a positive within-state intensity effect in an eligible state in at least two of the three cohorts and the contributing states meet the frozen cross-cohort marker-match rule. This is a post-lock biological-structure label and not formal replication.

All programs, states, directions, and failed gates are retained. No state is removed because it weakens the result.

## Module 5B: alcohol-hepatitis context specificity

### Prespecified comparisons

- Alcohol-associated hepatitis versus healthy liver.
- Alcohol-associated cirrhosis versus alcohol-associated hepatitis.

The existing 30-cell donor gate is primary and the 20-cell gate is sensitivity. All three target lineages have at least three eligible donors per group under the primary gate.

### Statistics

- Reuse the existing GSE256398 donor program scores and recompute contrast effects with Hedges g and HC3 robust intervals.
- For each program and score method, compare the observed effect with 1,000 frozen size/detection-matched random modules generated using the existing GSE256398 random-module procedure.
- Primary interpretation concerns Endo2 and SAEndo1. All 19 programs are reported as secondary.
- `CIRRHOSIS_CONTEXT_ENRICHED_POST_LOCK` requires a positive alcohol-cirrhosis-versus-hepatitis HC3 lower interval and matched-random 95th-percentile exceedance for both score methods. Hepatitis-versus-healthy is reported independently and is not required to be null.

Alcohol-associated hepatitis may coexist with fibrosis. Results are described as cirrhosis-versus-hepatitis context specificity, not pure fibrosis-versus-inflammation specificity.

## Module 5C: independent spatial scar validation

### Frozen resources

1. The 2022 full-transcriptome human cirrhosis spatial study of eight explants with histological fibrotic and parenchymal regions.
2. The 2025 E-MTAB-13132 experimental Visium cohort (three independent early-cirrhosis donors; one technical replicate retained only within donor) plus E-MTAB-14960 validation cohort (four independent fibrotic biopsies).

The 2025 linked E-MTAB-13130 snRNA-seq and E-MTAB-13131 snATAC-seq are not part of the primary spatial test.

### Audit and coverage gates

- Donor identity, tissue section identity, technology, expression matrix provenance, and scar/parenchyma or author-defined scar-cluster labels must be auditable.
- Technical sections from one donor are combined before inference.
- Primary gene-set coverage must be at least 80% in each spatial study: at least 6/7 Endo2 genes and 4/5 SAEndo1 genes. Secondary coherent-set coverage requires at least 14/17 and 5/6 genes.
- Stop and report `SPATIAL_RESOURCE_NOT_EVALUABLE` if expression matrices or region labels cannot be recovered without inventing annotations.

### Scoring and inference

- Score each spatial spot by within-section standardized mean and rank-based singscore using only measured frozen genes.
- Spots are measurements, not independent biological replicates.
- For every donor, compute the median scar-minus-parenchyma score and a within-donor standardized effect. Combine technical sections before generating the donor effect.
- Compare the observed donor-median effect with 10,000 detection-, abundance-, and set-size-matched random gene sets.
- A spatial study passes for a program when both score methods have a positive median donor effect, at least 75% of independent donors have positive effects, and both observed median effects exceed the matched-random 95th percentile.
- `SPATIALLY_REPLICATED_POST_LOCK` requires the pass rule in both independent spatial studies. No spot-level P value is used as replication evidence.
- Secondary analyses test spatial association with author-defined scar endothelial abundance and macrophage/mesenchymal scar abundance. These are localization/co-occurrence results only.

## Module 5D: conditional regulatory triangulation

E-MTAB-13131 snATAC-seq is analyzed only if at least one primary endothelial program is both `STATE_SUPPORTED_POST_LOCK` and `SPATIALLY_REPLICATED_POST_LOCK`.

If triggered, the frozen tests are:

1. promoter accessibility for the primary leading-edge genes in author-defined scar-associated versus sinusoidal endothelial states;
2. linked enhancer/peak-to-gene accessibility where author-provided peak-to-gene links are auditable;
3. transcription-factor motif activity associated with the accessible leading-edge loci;
4. agreement in direction with transcript-derived TF activity in at least two of GSE202379, GSE244832, and GSE256398.

All inference remains donor-aware. Same-donor RNA/ATAC agreement is orthogonal molecular support, not independent replication or causality. If the trigger is not met, the module is recorded as `NOT_TRIGGERED_BY_FROZEN_RULE` and no ATAC outcome is searched.

## Multiplicity, reporting, and stopping

- Primary tests are the two endothelial programs. All-program analyses are secondary and FDR-corrected within module and lineage.
- Both score methods must be shown; neither is selected after results are known.
- Every prespecified cohort, state, donor, and contrast is retained, including negative or discordant outcomes.
- New figures may be promoted to the main text only if they answer a distinct biological question. Otherwise they remain supplementary.
- Phase 5 stops after Modules 5A-5C and the conditional 5D decision. No additional enrichment, PPI, ROC/ML, unrestricted ligand-receptor, drug-prediction, or cross-sectional pseudotime analysis will be added.
