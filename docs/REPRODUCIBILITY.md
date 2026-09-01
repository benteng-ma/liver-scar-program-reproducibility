# Reproducibility map

## Frozen specification

Start with `config/`, `literature/program_inventory.csv`, `metadata/dataset_manifest.csv`, and `docs/analysis_plan.md`. Do not change program membership, cohort roles, scoring methods, or formal gates to match later outcomes.

## Analysis order

1. Run the cohort-specific extraction or aggregation scripts for GSE202379, GSE290642, GSE244832, GSE210077/Watson6, GSE181483, and GSE256398.
2. Run the matching `audit_*_gates.py` and lineage-audit scripts.
3. Run the cohort-specific `analyze_*_programs.py` scripts and matched random-module analyses.
4. Run `cross_cohort_synthesis.py` and `deep_transportability_benchmark.py`.
5. Run the Phase 3, Phase 4, and Phase 5 scripts in numerical order. The frozen Phase 2 classifications must not be overwritten by post-lock results.
6. Rebuild Figures 1-9 and Supplementary Figures S1-S10 with `make_benchmark_figures.py`, `make_deep_benchmark_figures.py`, and the Phase 3-5 figure scripts.
7. Run `python -m pytest -q` and compare generated source tables with `results/source_data/`.

Random seeds and iteration counts are recorded in the scripts and `results/logs/`. The public repository contains derived tables required to verify reported values but not the third-party input matrices.
