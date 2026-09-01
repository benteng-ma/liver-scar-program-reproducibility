from __future__ import annotations

import json
from pathlib import Path

import numpy as np
import pandas as pd
from scipy.io import mmread, mmwrite
from scipy.sparse import csc_matrix

from audit_gse290642_lineages import locate_one, read_genes


TARGET_LINEAGES = [
    "endothelial",
    "macrophage_monocyte",
    "mesenchymal_hsc_myofibroblast",
]


def main() -> None:
    repo = Path(__file__).resolve().parents[1]
    source = repo / "data" / "raw" / "GSE290642" / "human_matrices"
    interim = repo / "data" / "interim" / "GSE290642"
    assignments = pd.read_csv(interim / "retained_cell_assignments.csv.gz")
    mapping = pd.read_csv(repo / "metadata" / "gse290642_cluster_mapping.csv")
    cluster_to_lineage = {
        int(row.cluster): row.harmonized_lineage
        for row in mapping.itertuples()
        if row.target_included == "yes"
    }
    if sorted(set(cluster_to_lineage.values())) != sorted(TARGET_LINEAGES):
        raise RuntimeError("frozen mapping does not contain all three target lineages")
    assignments["harmonized_lineage"] = assignments["cluster"].map(cluster_to_lineage)

    donors = pd.read_csv(repo / "metadata" / "donor_manifest.csv", dtype=str).fillna("")
    donors = donors[donors["dataset_id"] == "GSE290642_human"].set_index("sample_id")
    sample_dirs = sorted(path for path in source.iterdir() if path.is_dir())

    sample_inputs: list[tuple[str, Path, Path, Path, list[str]]] = []
    union_symbols: set[str] = set()
    sample_symbol_sets: list[set[str]] = []
    for sample_dir in sample_dirs:
        gsm = sample_dir.name.split("_", 1)[0]
        matrix_path = locate_one(sample_dir, ("matrix.mtx", "matrix.mtx.gz"))
        gene_path = locate_one(sample_dir, ("genes.tsv", "genes.tsv.gz", "features.tsv", "features.tsv.gz"))
        barcode_path = locate_one(sample_dir, ("barcodes.tsv", "barcodes.tsv.gz"))
        _, symbols = read_genes(gene_path)
        upper_symbols = [symbol.upper() for symbol in symbols]
        valid_symbols = {symbol for symbol in upper_symbols if symbol and symbol != "NAN"}
        union_symbols.update(valid_symbols)
        sample_symbol_sets.append(valid_symbols)
        sample_inputs.append((gsm, matrix_path, gene_path, barcode_path, upper_symbols))
    union_genes = sorted(union_symbols)
    shared_genes = sorted(set.intersection(*sample_symbol_sets))
    union_index = {gene: index for index, gene in enumerate(union_genes)}

    columns: list[np.ndarray] = []
    manifest_rows: list[dict[str, object]] = []
    retained_cluster_rows: list[dict[str, object]] = []
    for sample_number, (gsm, matrix_path, gene_path, barcode_path, symbols) in enumerate(sample_inputs, start=1):
        matrix = mmread(matrix_path).tocsr()
        barcodes = pd.read_csv(barcode_path, sep="\t", header=None, dtype=str).iloc[:, 0].astype(str)
        if matrix.shape != (len(symbols), len(barcodes)):
            raise RuntimeError(f"matrix/features/barcodes disagree for {gsm}")
        sample_assignments = assignments[assignments["sample_id"] == gsm].copy()
        if sample_assignments["barcode"].duplicated().any():
            raise RuntimeError(f"duplicate retained barcode in {gsm}")
        barcode_to_position = pd.Series(np.arange(len(barcodes)), index=barcodes)
        if not sample_assignments["barcode"].isin(barcode_to_position.index).all():
            raise RuntimeError(f"retained barcode is absent from raw matrix for {gsm}")
        sample_assignments["matrix_column"] = sample_assignments["barcode"].map(barcode_to_position).astype(int)
        for cluster, values in sample_assignments.groupby("cluster"):
            retained_cluster_rows.append(
                {
                    "sample_id": gsm,
                    "donor_id": donors.loc[gsm, "donor_id"],
                    "cluster": int(cluster),
                    "harmonized_lineage": cluster_to_lineage.get(int(cluster), "excluded"),
                    "n_cells": len(values),
                }
            )
        symbol_positions = np.fromiter((union_index[symbol] for symbol in symbols), dtype=np.int64)
        donor = donors.loc[gsm]
        for lineage in TARGET_LINEAGES:
            cell_columns = sample_assignments.loc[
                sample_assignments["harmonized_lineage"] == lineage, "matrix_column"
            ].to_numpy(dtype=np.int64)
            gene_row_counts = np.asarray(matrix[:, cell_columns].sum(axis=1)).ravel()
            union_counts = np.bincount(
                symbol_positions,
                weights=gene_row_counts,
                minlength=len(union_genes),
            ).astype(np.int64)
            columns.append(union_counts)
            manifest_rows.append(
                {
                    "group_id": f"G{len(manifest_rows) + 1:04d}",
                    "donor_id": donor["donor_id"],
                    "sample_id": gsm,
                    "harmonized_lineage": lineage,
                    "n_cells": len(cell_columns),
                    "disease_group": donor["disease_group"],
                    "etiology": donor["etiology"],
                    "fibrosis_stage": donor["fibrosis_stage"],
                    "age": donor["age"],
                    "sex": donor["sex"],
                }
            )
        print(f"[{sample_number:02d}/24] {gsm} target counts aggregated", flush=True)
        del matrix

    aggregated = np.column_stack(columns)
    manifest = pd.DataFrame(manifest_rows)
    if aggregated.shape != (len(union_genes), 72) or len(manifest) != 72:
        raise RuntimeError("expected union genes x 24 donors x 3 lineages")
    interim.mkdir(parents=True, exist_ok=True)
    mmwrite(interim / "donor_lineage_raw_counts.mtx", csc_matrix(aggregated))
    pd.DataFrame({"gene": union_genes}).to_csv(interim / "genes.csv", index=False)
    pd.DataFrame({"gene": shared_genes}).to_csv(interim / "shared_genes.csv", index=False)
    manifest.to_csv(interim / "donor_lineage_manifest.csv", index=False)
    pd.DataFrame(retained_cluster_rows).to_csv(interim / "donor_cluster_cell_counts.csv", index=False)
    summary = {
        "dataset_id": "GSE290642_human",
        "genes_union": len(union_genes),
        "genes_shared_all_24_donors": len(shared_genes),
        "donors": int(manifest["donor_id"].nunique()),
        "target_cells": int(manifest["n_cells"].sum()),
        "target_donor_lineage_groups": len(manifest),
        "mapping_file": "metadata/gse290642_cluster_mapping.csv",
        "annotation_role": "reconstructed broad-lineage sensitivity only",
        "contrast_roles": {
            "F4_vs_F0": "small directional/sensitivity, three donors per group before lineage gate",
            "all_fibrosis_vs_F0": "sensitivity only",
        },
    }
    (interim / "donor_aggregation_summary.json").write_text(
        json.dumps(summary, indent=2) + "\n", encoding="utf-8"
    )
    print(manifest.groupby(["harmonized_lineage", "fibrosis_stage"])["n_cells"].agg(["count", "min", "median", "max"]).to_string())
    print(json.dumps(summary, indent=2))


if __name__ == "__main__":
    main()
