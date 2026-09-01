#!/usr/bin/env python3
"""Build a privacy-screened, analysis-focused GitHub release candidate."""

from __future__ import annotations

import csv
import hashlib
import json
import shutil
from pathlib import Path


ROOT = Path(__file__).resolve().parents[1]
RELEASE_ROOT = ROOT / "release"
OUT = RELEASE_ROOT / "liver-scar-program-reproducibility"
VERSION = "1.0.0"
TITLE = "Scar cell programs show limited reproducibility across human liver fibrosis cohorts despite a recurrent endothelial response"


def sha256(path: Path) -> str:
    digest = hashlib.sha256()
    with path.open("rb") as handle:
        for chunk in iter(lambda: handle.read(1024 * 1024), b""):
            digest.update(chunk)
    return digest.hexdigest()


def safe_reset() -> None:
    release_root = RELEASE_ROOT.resolve()
    target = OUT.resolve()
    if target.parent != release_root or target.name != "liver-scar-program-reproducibility":
        raise RuntimeError(f"Refusing to reset unexpected path: {target}")
    RELEASE_ROOT.mkdir(parents=True, exist_ok=True)
    if OUT.exists():
        shutil.rmtree(OUT)
    OUT.mkdir(parents=True)


def copy_file(relative: str, destination: str | None = None) -> None:
    source = ROOT / relative
    target = OUT / (destination or relative)
    target.parent.mkdir(parents=True, exist_ok=True)
    shutil.copy2(source, target)


def copy_tree(relative: str, exclude_names: set[str] | None = None) -> None:
    exclude_names = exclude_names or set()
    source_root = ROOT / relative
    for source in source_root.rglob("*"):
        if not source.is_file() or source.name in exclude_names:
            continue
        target = OUT / source.relative_to(ROOT)
        target.parent.mkdir(parents=True, exist_ok=True)
        shutil.copy2(source, target)


def is_ephemeral_sync_file(path: Path) -> bool:
    """Return True for transient cloud-sync sidecars that are not project files."""
    name = path.name.lower()
    return ".baiduyun.uploading" in name or name.endswith(".tmp")


def write_release_documents() -> None:
    (OUT / ".gitattributes").write_text(
        "* text=auto eol=lf\n"
        "*.png binary\n*.pdf binary\n*.parquet binary\n*.gz binary\n*.zip binary\n"
        "*.xlsx binary\n*.xls binary\n*.tif binary\n*.tiff binary\n*.h5ad binary\n",
        encoding="utf-8",
    )
    (OUT / ".gitignore").write_text(
        "# Raw or controlled inputs\n"
        "data/\n*.h5ad\n*.h5\n*.hdf5\n*.rds\n*.RData\n*.mtx\n*.mtx.gz\n"
        "# Local environments and caches\n"
        ".venv/\n.venv_phase2/\n.conda_env/\n.tools/\n__pycache__/\n.pytest_cache/\n*.py[cod]\n"
        "# Cloud-sync and temporary sidecars\n"
        "*.baiduyun.*\n*.tmp\n"
        "# Manuscript, submission, and local QA artifacts\n"
        "submission/\nvalidation/\n*.docx\n*.zip\n~$*\n"
        "# Secrets and machine-local configuration\n"
        ".env\n.env.*\n*.pem\n*.key\n",
        encoding="utf-8",
    )
    (OUT / "README.md").write_text(
        "# Liver scar-program reproducibility benchmark\n\n"
        f"Code, frozen analysis specifications, and derived results supporting the manuscript **{TITLE}**.\n\n"
        f"Release candidate: v{VERSION} (2026-09-01). The repository is staged privately at https://github.com/benteng-ma/liver-scar-program-reproducibility pending the final public-release and Zenodo checks.\n\n"
        "## What this repository shows\n\n"
        "Nineteen published human liver scar-cell programs were evaluated with independent donors as the inferential unit across six single-cell or single-nucleus cohorts and two spatial resources. Intact programs showed limited cross-cohort reproducibility, while a recurrent endothelial response was detectable at the member-gene, within-state, multicellular-network, and scar-localization levels. The repository does not support a universal diagnostic, prognostic, causal, or treatment-selection claim.\n\n"
        "## Repository contents\n\n"
        "- `config/`: frozen programs, thresholds, cohort roles, and state definitions.\n"
        "- `literature/`: frozen program inventory and provenance-oriented literature records.\n"
        "- `metadata/`: public dataset manifests, donor mappings, coverage, and eligibility records.\n"
        "- `scripts/`: cohort processing, scoring, robustness analyses, synthesis, and figure generation.\n"
        "- `results/`: derived numerical results, checksums, source data, and publication figures.\n"
        "- `docs/`: data-access, reproducibility, analysis-plan, and claim-boundary documentation.\n"
        "- `tests/`: public integrity tests that do not depend on manuscript submission files.\n\n"
        "## Raw data\n\n"
        "Raw expression matrices and third-party processed objects are not redistributed. Download them from the source repositories listed in `metadata/dataset_manifest.csv` and place them under the local `data/` hierarchy described in `docs/DATA_ACCESS.md`. The `data/` directory is intentionally ignored by Git.\n\n"
        "## Environment\n\n"
        "```text\n"
        "conda env create -f environment.yml\n"
        "conda activate cirrhosis-scar-transportability\n"
        "python -m pip install -r requirements-release.txt\n"
        "python scripts/verify_environment.py\n"
        "```\n\n"
        "## Verification\n\n"
        "After the public datasets have been placed at the documented paths, run the cohort scripts in the order listed in `docs/REPRODUCIBILITY.md`. Repository-integrity tests can be run with:\n\n"
        "```text\npython -m pytest -q\n```\n\n"
        "The frozen Phase 2 labels remain 0/19 for within-cell-state replication, pan-cirrhotic transportability, and assay robustness. Post-lock analyses provide biological localization and prioritization without altering those labels.\n\n"
        "## Citation\n\n"
        "Use `CITATION.cff`. Cite the version-specific Zenodo DOI after the v1.0.0 release has been archived.\n\n"
        "## License\n\n"
        "Repository code is released under the MIT License; see `LICENSE`. This license does not relicense third-party source data, which remain governed by their original repositories and terms.\n",
        encoding="utf-8",
    )
    (OUT / "LICENSE").write_text(
        "MIT License\n\n"
        "Copyright (c) 2026 Bing Chen, Benteng Ma, Ting Cai, Xiao-ming Liu, and Fen Wang\n\n"
        "Permission is hereby granted, free of charge, to any person obtaining a copy\n"
        "of this software and associated documentation files (the \"Software\"), to deal\n"
        "in the Software without restriction, including without limitation the rights\n"
        "to use, copy, modify, merge, publish, distribute, sublicense, and/or sell\n"
        "copies of the Software, and to permit persons to whom the Software is\n"
        "furnished to do so, subject to the following conditions:\n\n"
        "The above copyright notice and this permission notice shall be included in all\n"
        "copies or substantial portions of the Software.\n\n"
        "THE SOFTWARE IS PROVIDED \"AS IS\", WITHOUT WARRANTY OF ANY KIND, EXPRESS OR\n"
        "IMPLIED, INCLUDING BUT NOT LIMITED TO THE WARRANTIES OF MERCHANTABILITY,\n"
        "FITNESS FOR A PARTICULAR PURPOSE AND NONINFRINGEMENT. IN NO EVENT SHALL THE\n"
        "AUTHORS OR COPYRIGHT HOLDERS BE LIABLE FOR ANY CLAIM, DAMAGES OR OTHER\n"
        "LIABILITY, WHETHER IN AN ACTION OF CONTRACT, TORT OR OTHERWISE, ARISING FROM,\n"
        "OUT OF OR IN CONNECTION WITH THE SOFTWARE OR THE USE OR OTHER DEALINGS IN THE\n"
        "SOFTWARE.\n",
        encoding="utf-8",
    )
    (OUT / "CITATION.cff").write_text(
        "cff-version: 1.2.0\n"
        "message: \"If you use this code, please cite the archived software release and the associated article.\"\n"
        "type: software\n"
        "title: \"Liver scar-program reproducibility benchmark\"\n"
        f"version: \"{VERSION}\"\n"
        "date-released: 2026-09-01\n"
        "license: MIT\n"
        "repository-code: \"https://github.com/benteng-ma/liver-scar-program-reproducibility\"\n"
        "authors:\n"
        "  - family-names: Chen\n    given-names: Bing\n    orcid: https://orcid.org/0000-0001-8828-9678\n"
        "  - family-names: Ma\n    given-names: Benteng\n    orcid: https://orcid.org/0000-0001-8795-7291\n"
        "  - family-names: Cai\n    given-names: Ting\n    orcid: https://orcid.org/0000-0002-8910-3289\n"
        "  - family-names: Liu\n    given-names: Xiao-ming\n"
        "  - family-names: Wang\n    given-names: Fen\n"
        "keywords:\n"
        "  - liver fibrosis\n  - cirrhosis\n  - single-cell RNA sequencing\n  - reproducibility\n  - endothelial cells\n",
        encoding="utf-8",
    )
    (OUT / "requirements-release.txt").write_text(
        "-r requirements-phase2.txt\npytest==9.1.1\n",
        encoding="utf-8",
    )
    docs = OUT / "docs"
    docs.mkdir(exist_ok=True)
    (docs / "DATA_ACCESS.md").write_text(
        "# Data access and redistribution boundary\n\n"
        "The analyses use public human liver datasets. Accession, publication, assay, donor-count, endpoint, annotation, and limitation fields are recorded in `metadata/dataset_manifest.csv`.\n\n"
        "Raw matrices and third-party processed objects are not included in this repository. Retrieve them from GEO, ArrayExpress/BioStudies, Dryad, or the linked study repository using the source accession and citation. Preserve downloaded files as read-only inputs under `data/raw/` or `data/external/`; generated intermediates belong under `data/interim/` and `data/processed/`. These directories are ignored by Git.\n\n"
        "Do not infer biological zero from missing values in derived tables. Missing fields indicate not applicable, unavailable, or not estimated according to the table schema. Public donor identifiers are retained only as source-study codes; no direct identifiers or re-identification keys are included.\n",
        encoding="utf-8",
    )
    (docs / "REPRODUCIBILITY.md").write_text(
        "# Reproducibility map\n\n"
        "## Frozen specification\n\n"
        "Start with `config/`, `literature/program_inventory.csv`, `metadata/dataset_manifest.csv`, and `docs/analysis_plan.md`. Do not change program membership, cohort roles, scoring methods, or formal gates to match later outcomes.\n\n"
        "## Analysis order\n\n"
        "1. Run the cohort-specific extraction or aggregation scripts for GSE202379, GSE290642, GSE244832, GSE210077/Watson6, GSE181483, and GSE256398.\n"
        "2. Run the matching `audit_*_gates.py` and lineage-audit scripts.\n"
        "3. Run the cohort-specific `analyze_*_programs.py` scripts and matched random-module analyses.\n"
        "4. Run `cross_cohort_synthesis.py` and `deep_transportability_benchmark.py`.\n"
        "5. Run the Phase 3, Phase 4, and Phase 5 scripts in numerical order. The frozen Phase 2 classifications must not be overwritten by post-lock results.\n"
        "6. Rebuild Figures 1-9 and Supplementary Figures S1-S10 with `make_benchmark_figures.py`, `make_deep_benchmark_figures.py`, and the Phase 3-5 figure scripts.\n"
        "7. Run `python -m pytest -q` and compare generated source tables with `results/source_data/`.\n\n"
        "Random seeds and iteration counts are recorded in the scripts and `results/logs/`. The public repository contains derived tables required to verify reported values but not the third-party input matrices.\n",
        encoding="utf-8",
    )
    (OUT / ".zenodo.json").write_text(
        json.dumps(
            {
                "title": "Liver scar-program reproducibility benchmark",
                "upload_type": "software",
                "version": VERSION,
                "publication_date": "2026-09-01",
                "description": "Code, frozen specifications, and derived results for a donor-level benchmark of published human liver scar-cell programs.",
                "access_right": "open",
                "creators": [
                    {"name": "Chen, Bing", "orcid": "0000-0001-8828-9678"},
                    {"name": "Ma, Benteng", "orcid": "0000-0001-8795-7291"},
                    {"name": "Cai, Ting", "orcid": "0000-0002-8910-3289"},
                    {"name": "Liu, Xiao-ming"},
                    {"name": "Wang, Fen"},
                ],
                "keywords": ["liver fibrosis", "cirrhosis", "single-cell RNA sequencing", "reproducibility", "endothelial cells"],
                "license": "mit",
            },
            indent=2,
        )
        + "\n",
        encoding="utf-8",
    )


def copy_public_content() -> None:
    for relative in ("environment.yml", "requirements-phase0.txt", "requirements-phase2.txt"):
        copy_file(relative)
    copy_tree("config")
    copy_tree("literature")

    for source in (ROOT / "metadata").glob("*"):
        if not source.is_file() or source.suffix.lower() == ".docx" or source.name == "author_funding_provenance.md":
            continue
        target = OUT / "metadata" / source.name
        target.parent.mkdir(parents=True, exist_ok=True)
        shutil.copy2(source, target)

    excluded_scripts = {
        "audit_references_v8.py", "compare_hepcomm_v14_v15_renders.py",
        "make_docx_qa_contact_sheets.py", "make_docx_qa_contact_sheets_v14.py",
        "make_v12_contact_sheets.py", "make_v13_contact_sheets.py", "rasterize_qa_pdfs.py",
    }
    for source in (ROOT / "scripts").glob("*"):
        if not source.is_file() or source.name in excluded_scripts or is_ephemeral_sync_file(source):
            continue
        if source.name.startswith("build_hepcomm_submission") or source.name.startswith("validate_hepcomm_submission"):
            continue
        target = OUT / "scripts" / source.name
        target.parent.mkdir(parents=True, exist_ok=True)
        shutil.copy2(source, target)

    for name in ("test_phase0_integrity.py", "test_phase2_integrity.py", "test_deep_benchmark.py"):
        copy_file(f"tests/{name}")
    phase3_test = (ROOT / "tests" / "test_phase3_integrity.py").read_text(encoding="utf-8")
    phase3_test = phase3_test.split("\ndef test_manuscript_preserves_author_and_declaration_requirements", 1)[0].rstrip() + "\n"
    (OUT / "tests" / "test_phase3_integrity.py").write_text(phase3_test, encoding="utf-8")

    report_map = {
        "analysis_plan.md": "docs/analysis_plan.md",
        "deviations.md": "docs/deviations.md",
        "limitations_registry.md": "docs/limitations_registry.md",
        "phase0_decision.md": "docs/phase0_decision.md",
        "phase2_scoring_implementation.md": "docs/phase2_scoring_implementation.md",
        "deep_transportability_benchmark_plan_2026-08-31.md": "docs/deep_benchmark_plan.md",
        "deep_transportability_benchmark_results.md": "docs/deep_benchmark_results.md",
        "phase3_manuscript_enrichment_plan_2026-08-31.md": "docs/phase3_plan.md",
        "phase3_manuscript_enrichment_results_2026-08-31.md": "docs/phase3_results.md",
        "phase4_biological_structure_plan_2026-08-31.md": "docs/phase4_plan.md",
        "phase4_biological_structure_results_2026-08-31.md": "docs/phase4_results.md",
        "phase5_analysis_plan_frozen_2026-09-01.md": "docs/phase5_plan_frozen.md",
        "phase5_results_and_gate_decision_2026-09-01.md": "docs/phase5_results.md",
        "cross_cohort_synthesis.md": "docs/cross_cohort_synthesis.md",
    }
    for source, destination in report_map.items():
        copy_file(f"reports/{source}", destination)
    copy_file("reports/deep_transportability_benchmark_plan_2026-08-31.md")
    copy_file("reports/phase0_decision.json")

    excluded_results = {"r_4.6.1_install.log", "phase2_environment_verification.json"}
    for source in (ROOT / "results").rglob("*"):
        if not source.is_file() or source.name in excluded_results:
            continue
        target = OUT / source.relative_to(ROOT)
        target.parent.mkdir(parents=True, exist_ok=True)
        shutil.copy2(source, target)

    environment_record = json.loads(
        (ROOT / "results" / "logs" / "phase2_environment_verification.json").read_text(encoding="utf-8")
    )
    environment_record["executable"] = "<LOCAL_ISOLATED_ENVIRONMENT>/python"
    (OUT / "results" / "logs" / "phase2_environment_verification.json").write_text(
        json.dumps(environment_record, indent=2) + "\n", encoding="utf-8"
    )


def write_manifest() -> Path:
    rows = []
    for path in sorted(OUT.rglob("*")):
        if (
            not path.is_file()
            or any(part in {".git", ".pytest_cache", "__pycache__"} for part in path.parts)
            or path.suffix.lower() in {".pyc", ".pyo"}
            or is_ephemeral_sync_file(path)
            or path.name == "SHA256_manifest.csv"
        ):
            continue
        rows.append((path.relative_to(OUT).as_posix(), path.stat().st_size, sha256(path)))
    manifest = OUT / "SHA256_manifest.csv"
    with manifest.open("w", newline="", encoding="utf-8") as handle:
        writer = csv.writer(handle)
        writer.writerow(["path", "bytes", "sha256"])
        writer.writerows(rows)
    return manifest


def main() -> None:
    safe_reset()
    write_release_documents()
    copy_public_content()
    manifest = write_manifest()
    files = [path for path in OUT.rglob("*") if path.is_file()]
    summary = {
        "path": str(OUT),
        "version": VERSION,
        "files": len(files),
        "bytes": sum(path.stat().st_size for path in files),
        "largest_file_bytes": max(path.stat().st_size for path in files),
        "manifest": str(manifest),
        "license_status": "MIT_CONFIRMED_PRIVATE_STAGING",
    }
    (RELEASE_ROOT / "public_repository_candidate_summary.json").write_text(json.dumps(summary, indent=2) + "\n", encoding="utf-8")
    print(json.dumps(summary, indent=2))


if __name__ == "__main__":
    main()
