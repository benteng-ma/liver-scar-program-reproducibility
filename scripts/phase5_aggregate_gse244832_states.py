from __future__ import annotations

from pathlib import Path

import numpy as np
import pandas as pd
from scipy.io import mmwrite
from scipy.sparse import csc_matrix

from aggregate_gse244832_donors import canonical_donor_id, read_single_column
from stream_matrix_market import read_header


def aggregate_by_cell_group_fast(
    path: Path,
    cell_groups: np.ndarray,
    n_groups: int,
    chunksize: int = 5_000_000,
) -> tuple[np.ndarray, int]:
    """Aggregate an uncompressed coordinate Matrix Market file with pandas' C parser."""
    n_genes, n_cells, n_nonzero, _ = read_header(path)
    if cell_groups.shape != (n_cells,):
        raise ValueError("cell group length does not match matrix columns")
    flat = np.zeros(n_genes * n_groups, dtype=np.int64)
    observed = 0
    reader = pd.read_csv(
        path,
        sep=" ",
        skiprows=2,
        names=["row", "column", "value"],
        dtype={"row": "int32", "column": "int32", "value": "int32"},
        chunksize=chunksize,
        engine="c",
        compression=None,
    )
    for chunk in reader:
        rows = chunk["row"].to_numpy(np.int64) - 1
        columns = chunk["column"].to_numpy(np.int64) - 1
        values = chunk["value"].to_numpy(np.int64)
        indices = rows * n_groups + cell_groups[columns]
        flat += np.bincount(indices, weights=values, minlength=flat.size).astype(np.int64)
        observed += len(chunk)
    if observed != n_nonzero:
        raise RuntimeError(f"observed {observed} coordinates; header declares {n_nonzero}")
    return flat.reshape(n_genes, n_groups), observed


def main() -> None:
    repo = Path(__file__).resolve().parents[1]
    source = repo / "data" / "raw" / "GSE244832" / "processed_files"
    output = repo / "data" / "interim" / "GSE244832"
    genes = read_single_column(source / "hLIVER_genes.csv")
    cells = read_single_column(source / "hLIVER_cells.csv")
    metadata = pd.read_csv(source / "hLIVER_metadata.csv", index_col=0)
    mapping = pd.read_csv(repo / "metadata" / "gse244832_cluster_mapping.csv")
    if cells != metadata.index.tolist():
        raise RuntimeError("cell order does not match author metadata")
    n_genes, n_cells, n_nonzero, _ = read_header(source / "hLIVER_counts.mtx.gz")
    if (n_genes, n_cells) != (len(genes), len(metadata)):
        raise RuntimeError("matrix dimensions do not match metadata")

    cluster_to_lineage = {
        int(row.seurat_cluster): row.harmonized_lineage
        for row in mapping.itertuples()
        if row.target_included == "yes"
    }
    metadata["harmonized_lineage"] = metadata["seurat_clusters"].map(cluster_to_lineage)
    target = metadata["harmonized_lineage"].notna()
    metadata["source_state"] = "cluster_" + metadata["seurat_clusters"].astype(str)
    keys = (
        metadata.loc[target, "orig.ident"].astype(str)
        + "|||"
        + metadata.loc[target, "harmonized_lineage"].astype(str)
        + "|||"
        + metadata.loc[target, "source_state"].astype(str)
    )
    levels = sorted(
        keys.unique(),
        key=lambda value: (
            int(value.split("|||", 1)[0][2:]),
            value.split("|||", 2)[1],
            int(value.rsplit("_", 1)[1]),
        ),
    )
    all_levels = [*levels, "__EXCLUDED__"]
    lookup = {value: index for index, value in enumerate(all_levels)}
    cell_keys = pd.Series("__EXCLUDED__", index=metadata.index)
    cell_keys.loc[target] = keys
    codes = cell_keys.map(lookup).to_numpy(np.int64)

    aggregated, observed = aggregate_by_cell_group_fast(
        source / "hLIVER_counts.mtx.gz", codes, len(all_levels)
    )
    if observed != n_nonzero:
        raise RuntimeError("Matrix Market nonzero count mismatch")
    aggregated = aggregated[:, : len(levels)]
    group_ids = [f"S{index:05d}" for index in range(1, len(levels) + 1)]
    mmwrite(output / "phase5_donor_state_raw_counts.mtx", csc_matrix(aggregated))
    pd.DataFrame({"gene": genes}).to_csv(output / "phase5_state_genes.csv", index=False)

    rows: list[dict[str, object]] = []
    for index, (group_id, key) in enumerate(zip(group_ids, levels)):
        author_donor, lineage, source_state = key.split("|||", 2)
        selected = (
            target
            & metadata["orig.ident"].eq(author_donor)
            & metadata["harmonized_lineage"].eq(lineage)
            & metadata["source_state"].eq(source_state)
        )
        conditions = metadata.loc[selected, "condition"].astype(str).unique()
        if len(conditions) != 1:
            raise RuntimeError(f"condition is not constant for {key}")
        condition = conditions[0]
        rows.append(
            {
                "state_group_id": group_id,
                "donor_id": canonical_donor_id(author_donor),
                "author_donor_id": author_donor,
                "harmonized_lineage": lineage,
                "source_state": source_state,
                "n_cells": int(selected.sum()),
                "library_size": int(aggregated[:, index].sum()),
                "condition": condition,
                "disease_group": {"NORMAL": "normal", "NAFL": "MASL", "NASH": "MASH"}[condition],
                "etiology": "none" if condition == "NORMAL" else "metabolic",
                "fibrosis_stage": "group-level F2-F4 only" if condition == "NASH" else "not reported",
            }
        )
    manifest = pd.DataFrame(rows)
    manifest.to_csv(output / "phase5_donor_state_manifest.csv", index=False)
    print(
        manifest.groupby(["harmonized_lineage", "source_state", "disease_group"])["n_cells"]
        .agg(["count", "min", "median", "max"])
        .to_string()
    )


if __name__ == "__main__":
    main()
