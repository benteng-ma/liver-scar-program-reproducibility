from __future__ import annotations

import hashlib
import json
from pathlib import Path

import numpy as np
import pandas as pd

from stream_matrix_market import aggregate_by_cell_group, read_header


AUTHOR_MARKERS = {
    "hepatocyte": ["HNF4A", "ALB", "APOA1"],
    "mesenchymal_hsc_myofibroblast": [
        "NGFR",
        "COL1A1",
        "COL1A2",
        "DCN",
        "PDGFRB",
        "RGS5",
        "ACTA2",
        "MYH11",
    ],
    "endothelial": ["PECAM1", "VWF", "EMCN", "KDR", "ENG"],
    "t_cell": ["CD69", "CD3D", "CD3E", "TRBC1"],
    "b_cell": ["CD19", "MS4A1", "CD79A", "CD74"],
    "macrophage_monocyte": ["CD163", "LST1", "TYROBP", "FCER1G", "C1QA"],
    "cholangiocyte": ["KRT19", "KRT7", "KRT8", "KRT18", "EPCAM"],
}


def sha256(path: Path) -> str:
    digest = hashlib.sha256()
    with path.open("rb") as handle:
        for block in iter(lambda: handle.read(8 * 1024 * 1024), b""):
            digest.update(block)
    return digest.hexdigest().upper()


def read_single_column(path: Path) -> list[str]:
    return [line.strip().strip('"') for line in path.read_text().splitlines() if line.strip()]


def main() -> None:
    repo = Path(__file__).resolve().parents[1]
    source = repo / "data" / "raw" / "GSE244832" / "processed_files"
    matrix_path = source / "hLIVER_counts.mtx.gz"
    genes = read_single_column(source / "hLIVER_genes.csv")
    cells = read_single_column(source / "hLIVER_cells.csv")
    metadata = pd.read_csv(source / "hLIVER_metadata.csv", index_col=0)
    if cells != metadata.index.tolist():
        raise RuntimeError("hLIVER cell order does not match metadata row order")

    cluster_labels = sorted(metadata["seurat_clusters"].astype(str).unique(), key=int)
    cluster_to_code = {label: index for index, label in enumerate(cluster_labels)}
    cell_groups = (
        metadata["seurat_clusters"].astype(str).map(cluster_to_code).to_numpy(np.int64)
    )
    n_genes, n_cells, n_nonzero, _ = read_header(matrix_path)
    if n_genes != len(genes) or n_cells != len(cells):
        raise RuntimeError("matrix dimensions do not match gene/cell files")

    counts, observed_nonzero = aggregate_by_cell_group(
        matrix_path, cell_groups, len(cluster_labels)
    )
    library = counts.sum(axis=0)
    cpm = counts / np.maximum(library, 1) * 1_000_000
    log_cpm = np.log1p(cpm)
    gene_to_index = {gene: index for index, gene in enumerate(genes)}

    rows = []
    for code, cluster in enumerate(cluster_labels):
        cell_subset = metadata[metadata["seurat_clusters"].astype(str) == cluster]
        row = {
            "cluster": cluster,
            "n_cells": len(cell_subset),
            "n_donors": cell_subset["orig.ident"].nunique(),
            "conditions": ";".join(sorted(cell_subset["condition"].unique())),
            "library_size": int(library[code]),
        }
        for lineage, markers in AUTHOR_MARKERS.items():
            available = [gene for gene in markers if gene in gene_to_index]
            values = [log_cpm[gene_to_index[gene], code] for gene in available]
            row[f"score_{lineage}"] = float(np.mean(values))
            for gene in available:
                row[f"log1p_cpm_{gene}"] = float(log_cpm[gene_to_index[gene], code])
        score_columns = [key for key in row if key.startswith("score_")]
        winner = max(score_columns, key=lambda key: row[key])
        ordered = sorted((row[key], key) for key in score_columns)
        row["provisional_marker_label"] = winner.removeprefix("score_")
        row["top_minus_second_score"] = float(ordered[-1][0] - ordered[-2][0])
        rows.append(row)

    audit = pd.DataFrame(rows).sort_values("cluster", key=lambda col: col.astype(int))
    qc_dir = repo / "results" / "qc"
    interim_dir = repo / "data" / "interim" / "GSE244832"
    qc_dir.mkdir(parents=True, exist_ok=True)
    interim_dir.mkdir(parents=True, exist_ok=True)
    audit.to_csv(qc_dir / "gse244832_cluster_marker_audit.csv", index=False)
    pd.DataFrame(
        counts,
        index=genes,
        columns=[f"cluster_{cluster}" for cluster in cluster_labels],
    ).to_csv(interim_dir / "gene_by_cluster_raw_counts.csv.gz", compression="gzip")
    summary = {
        "dataset": "GSE244832",
        "cells": len(cells),
        "genes": len(genes),
        "clusters": len(cluster_labels),
        "donors": int(metadata["orig.ident"].nunique()),
        "condition_counts": metadata.groupby("condition")["orig.ident"].nunique().to_dict(),
        "matrix_declared_nonzero": n_nonzero,
        "matrix_observed_nonzero": observed_nonzero,
        "matrix_sha256": sha256(matrix_path),
        "metadata_sha256": sha256(source / "hLIVER_metadata.csv"),
        "mapping_status": "provisional marker reconstruction; manual freeze required",
        "author_marker_source": "Kim et al. PMID 39522884 main text Figure 1B-C",
    }
    (qc_dir / "gse244832_ingest_summary.json").write_text(
        json.dumps(summary, indent=2) + "\n", encoding="utf-8"
    )
    print(audit[["cluster", "n_cells", "n_donors", "provisional_marker_label", "top_minus_second_score"]].to_string(index=False))
    print(json.dumps(summary, indent=2))


if __name__ == "__main__":
    main()
