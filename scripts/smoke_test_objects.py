"""Read one real scRNA, snRNA and spatial object without disease testing."""

from __future__ import annotations

import gzip
import hashlib
import json
import platform
from pathlib import Path

import numpy as np
import pandas as pd
import scipy
from scipy.io import mmread


ROOT = Path(__file__).resolve().parents[1]
OUT = ROOT / "results" / "qc"
OUT.mkdir(parents=True, exist_ok=True)


def sha256(path: Path) -> str:
    h = hashlib.sha256()
    with path.open("rb") as handle:
        while chunk := handle.read(1024 * 1024):
            h.update(chunk)
    return h.hexdigest()


def read_10x(folder: Path, prefix: str, feature_token: str) -> tuple[dict, set[str]]:
    matrix_path = folder / f"{prefix}_matrix.mtx.gz"
    barcode_path = folder / f"{prefix}_barcodes.tsv.gz"
    feature_path = folder / f"{prefix}_{feature_token}.tsv.gz"
    matrix = mmread(matrix_path).tocsr()
    with gzip.open(barcode_path, "rt", encoding="utf-8") as handle:
        barcodes = [line.rstrip("\n") for line in handle]
    features = pd.read_csv(feature_path, sep="\t", header=None, compression="gzip")
    symbol_col = 1 if features.shape[1] > 1 else 0
    genes = set(features.iloc[:, symbol_col].dropna().astype(str))
    assert matrix.shape == (len(features), len(barcodes))
    return {
        "matrix_shape_genes_by_cells": list(matrix.shape),
        "nonzero_entries": int(matrix.nnz),
        "sparse_dtype": str(matrix.dtype),
        "n_barcodes": len(barcodes),
        "n_features": len(features),
        "matrix_sha256": sha256(matrix_path),
        "read_status": "PASS",
    }, genes


scrna_dir = ROOT / "data" / "external" / "smoke_scrna_GSE136103"
snrna_dir = ROOT / "data" / "external" / "smoke_snrna_GSE210077"
spatial_dir = ROOT / "data" / "external" / "smoke_spatial_GSE210077"

scrna, scrna_genes = read_10x(scrna_dir, "GSM4041150_healthy1_cd45+", "genes")
snrna, snrna_genes = read_10x(snrna_dir, "GSM6416560_6854-9", "features")

spatial_path = spatial_dir / "cell_properties_healthy_diseased_merfish.head_1MiB.csv"
spatial = pd.read_csv(spatial_path, on_bad_lines="skip")
required = {"donor_id", "condition", "Cell_Type_final", "x", "y", "n_counts"}
assert required.issubset(spatial.columns)
spatial_genes = set(pd.read_csv(
    spatial_dir / "gene_names_healthy_diseased_merfish.csv", header=None
).iloc[:, 0].dropna().astype(str))
source_data_path = ROOT / "data" / "raw" / "phase0_metadata" / "PMC11697218_supplementary" / "41467_2024_55325_MOESM11_ESM.xlsx"
spatial_source = pd.read_excel(source_data_path, sheet_name="Figure 4")
source_required = {"cell_id", "cell_name", "nuclei_num", "area", "n_counts", "count_density"}
assert source_required.issubset(spatial_source.columns)
spatial_result = {
    "object": "complete publication spatial source-data table plus HCA cell-properties bounded 1 MiB coordinate range and complete MERFISH gene panel",
    "complete_source_table_rows": int(len(spatial_source)),
    "complete_source_table_sha256": sha256(source_data_path),
    "rows_read": int(len(spatial)),
    "columns": list(spatial.columns),
    "coordinate_finite_fraction": float(np.isfinite(spatial[["x", "y"]].to_numpy()).mean()),
    "donors_in_bounded_range": sorted(spatial["donor_id"].dropna().astype(str).unique().tolist()),
    "cell_labels_in_bounded_range": sorted(spatial["Cell_Type_final"].dropna().astype(str).unique().tolist()),
    "gene_panel_size": len(spatial_genes),
    "range_file_sha256": sha256(spatial_path),
    "read_status": "PASS_COMPLETE_SOURCE_TABLE_PLUS_BOUNDED_COORDINATE_RANGE_NOT_FULL_H5AD",
}

programs = pd.read_csv(ROOT / "literature" / "program_inventory.csv")
coverage_rows = []
for program_id, group in programs.groupby("program_id", sort=False):
    genes = set(group["gene_symbol"].dropna().astype(str))
    for assay, detected in [
        ("scRNA_GSE136103_feature_space", scrna_genes),
        ("snRNA_GSE210077_feature_space", snrna_genes),
        ("MERFISH_GSE210077_panel", spatial_genes),
    ]:
        n = len(genes & detected)
        fraction = n / len(genes) if genes else 0.0
        if assay.startswith("MERFISH"):
            eligible = n >= 5 and fraction >= 0.20
            tier = "spatial_evaluable" if eligible else "spatial_ineligible"
        elif fraction >= 0.80:
            tier = "primary"
        elif fraction >= 0.60:
            tier = "evaluable_flagged"
        elif fraction >= 0.40:
            tier = "sensitivity_only"
        else:
            tier = "not_evaluable"
        coverage_rows.append({
            "program_id": program_id,
            "cell_lineage": group.iloc[0]["cell_lineage"],
            "assay_resource": assay,
            "n_program_genes": len(genes),
            "n_detected": n,
            "coverage": round(fraction, 6),
            "coverage_tier": tier,
        })
pd.DataFrame(coverage_rows).to_csv(ROOT / "metadata" / "program_assay_coverage.csv", index=False, encoding="utf-8-sig")

result = {
    "generated": "2026-08-30",
    "prohibition": "No disease effect, clustering, integration or differential expression was performed.",
    "runtime": {
        "python": platform.python_version(),
        "pandas": pd.__version__,
        "scipy": scipy.__version__,
        "numpy": np.__version__,
    },
    "scRNA_GSE136103_GSM4041150": scrna,
    "snRNA_GSE210077_GSM6416560": snrna,
    "spatial_GSE210077_HCA": spatial_result,
    "coverage_summary": {
        "n_frozen_programs": int(programs["program_id"].nunique()),
        "scRNA_primary_coverage_programs": sum(r["coverage_tier"] == "primary" for r in coverage_rows if r["assay_resource"].startswith("scRNA")),
        "snRNA_primary_coverage_programs": sum(r["coverage_tier"] == "primary" for r in coverage_rows if r["assay_resource"].startswith("snRNA")),
        "spatial_evaluable_programs": sum(r["coverage_tier"] == "spatial_evaluable" for r in coverage_rows),
    },
}
(OUT / "phase0_smoke_test.json").write_text(json.dumps(result, indent=2, ensure_ascii=False), encoding="utf-8")
print(json.dumps(result["coverage_summary"], indent=2))
