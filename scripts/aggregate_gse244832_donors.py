from __future__ import annotations

import json
from pathlib import Path

import numpy as np
import pandas as pd
from scipy.io import mmwrite
from scipy.sparse import csc_matrix

from stream_matrix_market import aggregate_by_cell_group, read_header


def read_single_column(path: Path) -> list[str]:
    return [line.strip().strip('"') for line in path.read_text().splitlines() if line.strip()]


def canonical_donor_id(author_id: str) -> str:
    if not author_id.startswith("JB"):
        raise ValueError(f"unexpected donor ID: {author_id}")
    return f"JB_{author_id[2:]}"


def main() -> None:
    repo = Path(__file__).resolve().parents[1]
    source = repo / "data" / "raw" / "GSE244832" / "processed_files"
    output = repo / "data" / "interim" / "GSE244832"
    output.mkdir(parents=True, exist_ok=True)
    matrix_path = source / "hLIVER_counts.mtx.gz"
    genes = read_single_column(source / "hLIVER_genes.csv")
    cells = read_single_column(source / "hLIVER_cells.csv")
    metadata = pd.read_csv(source / "hLIVER_metadata.csv", index_col=0)
    mapping = pd.read_csv(repo / "metadata" / "gse244832_cluster_mapping.csv")
    if cells != metadata.index.tolist():
        raise RuntimeError("cell order does not match author metadata")
    n_genes, n_cells, n_nonzero, _ = read_header(matrix_path)
    if n_genes != len(genes) or n_cells != len(metadata):
        raise RuntimeError("matrix dimensions do not match genes/metadata")

    cluster_to_lineage = {
        int(row.seurat_cluster): row.harmonized_lineage
        for row in mapping.itertuples()
        if row.target_included == "yes"
    }
    metadata["harmonized_lineage"] = metadata["seurat_clusters"].map(cluster_to_lineage)
    target = metadata["harmonized_lineage"].notna()
    target_keys = (
        metadata.loc[target, "orig.ident"].astype(str)
        + "|||"
        + metadata.loc[target, "harmonized_lineage"].astype(str)
    )
    group_levels = sorted(
        target_keys.unique(),
        key=lambda key: (int(key.split("|||", 1)[0][2:]), key.split("|||", 1)[1]),
    )
    all_levels = [*group_levels, "__EXCLUDED__"]
    level_to_code = {value: index for index, value in enumerate(all_levels)}
    cell_keys = pd.Series("__EXCLUDED__", index=metadata.index)
    cell_keys.loc[target] = target_keys
    cell_codes = cell_keys.map(level_to_code).to_numpy(np.int64)

    aggregated, observed_nonzero = aggregate_by_cell_group(
        matrix_path, cell_codes, len(all_levels)
    )
    if observed_nonzero != n_nonzero:
        raise RuntimeError("Matrix Market nonzero count mismatch")
    pseudobulk = aggregated[:, : len(group_levels)]
    group_ids = [f"G{index:04d}" for index in range(1, len(group_levels) + 1)]
    mmwrite(output / "donor_lineage_raw_counts.mtx", csc_matrix(pseudobulk))
    pd.DataFrame({"gene": genes}).to_csv(output / "genes.csv", index=False)

    manifest_rows = []
    for group_id, key in zip(group_ids, group_levels):
        author_donor, lineage = key.split("|||", 1)
        selected = target & metadata["orig.ident"].eq(author_donor) & metadata[
            "harmonized_lineage"
        ].eq(lineage)
        conditions = sorted(metadata.loc[selected, "condition"].astype(str).unique())
        if len(conditions) != 1:
            raise RuntimeError(f"condition is not constant for {key}")
        condition = conditions[0]
        manifest_rows.append(
            {
                "group_id": group_id,
                "donor_id": canonical_donor_id(author_donor),
                "author_donor_id": author_donor,
                "harmonized_lineage": lineage,
                "n_cells": int(selected.sum()),
                "condition": condition,
                "disease_group": {
                    "NORMAL": "normal",
                    "NAFL": "MASL",
                    "NASH": "MASH",
                }[condition],
                "etiology": "none" if condition == "NORMAL" else "metabolic",
                "fibrosis_stage": (
                    "group-level F2-F4 only" if condition == "NASH" else "not reported"
                ),
            }
        )
    group_manifest = pd.DataFrame(manifest_rows)
    if len(group_manifest) != 54:
        raise RuntimeError("expected 18 donors x 3 target lineages")
    group_manifest.to_csv(output / "donor_lineage_manifest.csv", index=False)

    cluster_counts = (
        metadata.groupby(["orig.ident", "condition", "seurat_clusters"])
        .size()
        .reset_index(name="n_cells")
    )
    cluster_counts.to_csv(output / "donor_cluster_cell_counts.csv", index=False)
    summary = {
        "dataset_id": "GSE244832",
        "matrix_declared_nonzero": n_nonzero,
        "matrix_observed_nonzero": observed_nonzero,
        "genes": len(genes),
        "cells": len(metadata),
        "donors": metadata["orig.ident"].nunique(),
        "target_cells": int(target.sum()),
        "target_donor_lineage_groups": len(group_manifest),
        "condition_donors": (
            metadata.drop_duplicates("orig.ident").groupby("condition").size().to_dict()
        ),
        "mapping_file": "metadata/gse244832_cluster_mapping.csv",
        "endpoint_status": "MASH-vs-NORMAL sensitivity only; no donor-level fibrosis stage",
    }
    (output / "donor_aggregation_summary.json").write_text(
        json.dumps(summary, indent=2) + "\n", encoding="utf-8"
    )
    print(group_manifest.groupby(["harmonized_lineage", "disease_group"])["n_cells"].agg(["count", "min", "median", "max"]).to_string())
    print(json.dumps(summary, indent=2))


if __name__ == "__main__":
    main()
