# Phase 0 decision

## Decision: `CONDITIONAL_GO_NO_SPATIAL`

Proceed with a three-lineage donor-level transcriptomic benchmark. Remove spatial validation from the main results and title. Spatial data may be mentioned only as a coverage limitation; absence of MERFISH signal cannot be interpreted as biological absence.

## Prespecified criteria

| # | Criterion | Result | Evidence |
|---:|---|---|---|
| 1 | No exact duplicate | PASS | No paper met all eight exact-duplicate elements; closest 2026 study lacks donor benchmark/transfer/random modules/spatial validation |
| 2 | ≥3 independent human liver datasets | PASS | GSE181483, GSE202379, GSE210077-Watson, GSE244832 and GSE290642 are independent of GSE136103 for Ramachandran programs |
| 3 | ≥2 etiologies | PASS WITH RESTRICTION | metabolic/MASH and MetALD are available independently; HBV is macrophage-only secondary; unspecified cohorts cannot be used for etiologic claims |
| 4 | scRNA and snRNA | PASS | GSE136103/GSE181483/GSE290642 scRNA; GSE202379/GSE210077/GSE244832 snRNA |
| 5 | Usable human liver spatial resource | RESOURCE PASS | Watson MERFISH has 3 healthy and F2/F3/F4 donors with coordinates and labels |
| 6 | Recoverable donor IDs | PASS | Explicit donor IDs or sample-as-donor mappings recovered; repeated fractions/lobes identified |
| 7 | All three lineages in ≥2 independent datasets | PASS AT REPORTED-LINEAGE LEVEL | Whole-liver/NPC studies report all three; per-donor ≥30-cell eligibility is a frozen Phase 2 gate |
| 8 | About 8 healthy and 8 disease donors per lineage | PASS AT FEASIBILITY LEVEL | Approximate totals exceed threshold after donor collapsing; exact eligible totals will be reported before effects |
| 9 | At least one truly independent queue | PASS | Multiple validation accessions did not participate in GSE136103 program discovery |
| 10 | ≥15 exact direction-defined programs | PASS | 19 programs, 2,893 gene rows, exact supplementary provenance |
| 11 | Adequate validation gene coverage | PASS | All 19 reach ≥80% in tested scRNA and snRNA feature spaces |
| 12 | Main result does not require cell integration | PASS | Frozen donor×state pseudobulk and effect-only meta-analysis |
| 13 | Processed counts avoid FASTQ | PASS FOR CORE | Core scRNA/snRNA candidates provide processed counts; raw-only HBV resource is optional |
| 14 | Spatial panel adequately covers some main programs | **FAIL** | 0/19 satisfies both ≥5 detected genes and ≥20% coverage on the 317-gene MERFISH panel |

## Consequences

- Working title: **Cross-cohort transportability of macrophage, endothelial and mesenchymal scar programs in human liver fibrosis and cirrhosis**.
- Do not use “spatial”, “scar niche”, “pan-cirrhotic” or “universal” in the title before final classification supports it.
- Three transcriptomic lineages remain frozen. A lineage automatically drops from formal inference if the full processed data do not yield the preregistered donor×state cell counts; this is a feasibility gate, not a post-result choice.
- The analysis plan is frozen in `reports/analysis_plan.md` before any disease-effect analysis.

