# Liver scar-program reproducibility benchmark

[![DOI](https://zenodo.org/badge/DOI/10.5281/zenodo.22227547.svg)](https://doi.org/10.5281/zenodo.22227547)

Code, frozen analysis specifications, and derived results supporting the manuscript **Scar cell programs show limited reproducibility across human liver fibrosis cohorts despite a recurrent endothelial response**.

Release v1.0.0 (2026-09-01) is archived at Zenodo under [doi:10.5281/zenodo.22227547](https://doi.org/10.5281/zenodo.22227547). The GitHub repository is public at https://github.com/benteng-ma/liver-scar-program-reproducibility.

## What this repository shows

Nineteen published human liver scar-cell programs were evaluated with independent donors as the inferential unit across six single-cell or single-nucleus cohorts and two spatial resources. Intact programs showed limited cross-cohort reproducibility, while a recurrent endothelial response was detectable at the member-gene, within-state, multicellular-network, and scar-localization levels. The repository does not support a universal diagnostic, prognostic, causal, or treatment-selection claim.

## Repository contents

- `config/`: frozen programs, thresholds, cohort roles, and state definitions.
- `literature/`: frozen program inventory and provenance-oriented literature records.
- `metadata/`: public dataset manifests, donor mappings, coverage, and eligibility records.
- `scripts/`: cohort processing, scoring, robustness analyses, synthesis, and figure generation.
- `results/`: derived numerical results, checksums, source data, and publication figures.
- `docs/`: data-access, reproducibility, analysis-plan, and claim-boundary documentation.
- `tests/`: public integrity tests that do not depend on manuscript submission files.

## Raw data

Raw expression matrices and third-party processed objects are not redistributed. Download them from the source repositories listed in `metadata/dataset_manifest.csv` and place them under the local `data/` hierarchy described in `docs/DATA_ACCESS.md`. The `data/` directory is intentionally ignored by Git.

## Environment

```text
conda env create -f environment.yml
conda activate cirrhosis-scar-transportability
python -m pip install -r requirements-release.txt
python scripts/verify_environment.py
```

## Verification

After the public datasets have been placed at the documented paths, run the cohort scripts in the order listed in `docs/REPRODUCIBILITY.md`. Repository-integrity tests can be run with:

```text
python -m pytest -q
```

The frozen Phase 2 labels remain 0/19 for within-cell-state replication, pan-cirrhotic transportability, and assay robustness. Post-lock analyses provide biological localization and prioritization without altering those labels.

## Citation

Use `CITATION.cff` and cite the version-specific archive at https://doi.org/10.5281/zenodo.22227547. The all-versions concept DOI is https://doi.org/10.5281/zenodo.22227546.

## License

Repository code is released under the MIT License; see `LICENSE`. This license does not relicense third-party source data, which remain governed by their original repositories and terms.
