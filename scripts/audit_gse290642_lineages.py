from __future__ import annotations

import hashlib
import json
from pathlib import Path

import numpy as np
import pandas as pd
import yaml
from scipy.io import mmread
from sklearn.cluster import MiniBatchKMeans
from sklearn.decomposition import PCA
from sklearn.preprocessing import StandardScaler


SEED = 20260830
MIN_UMI = 500
MIN_GENES = 200
MAX_GENES = 6000
MAX_MT_FRACTION = 0.20
N_CLUSTERS = 40
N_PCS = 20
MIN_SCORE_MARGIN = 0.5
MIN_ANCHORS = 2
MIN_ANCHOR_Z = 0.5


def sha256(path: Path) -> str:
    digest = hashlib.sha256()
    with path.open("rb") as handle:
        for block in iter(lambda: handle.read(8 * 1024 * 1024), b""):
            digest.update(block)
    return digest.hexdigest().upper()


def locate_one(sample_dir: Path, names: tuple[str, ...]) -> Path:
    matches = [path for name in names for path in sample_dir.rglob(name)]
    if len(matches) != 1:
        raise RuntimeError(
            f"Expected exactly one of {names} below {sample_dir}, found {len(matches)}"
        )
    return matches[0]


def read_genes(path: Path) -> tuple[list[str], list[str]]:
    table = pd.read_csv(path, sep="\t", header=None, dtype=str)
    if table.shape[1] < 2:
        raise RuntimeError(f"Expected two-column genes.tsv: {path}")
    return table.iloc[:, 0].tolist(), table.iloc[:, 1].tolist()


def marker_expression(
    matrix,
    symbols: list[str],
    marker_order: list[str],
    keep: np.ndarray,
    retained_library: np.ndarray,
) -> np.ndarray:
    symbol_to_rows: dict[str, list[int]] = {}
    marker_set = set(marker_order)
    for row, symbol in enumerate(symbols):
        if symbol in marker_set:
            symbol_to_rows.setdefault(symbol, []).append(row)

    result = np.zeros((int(keep.sum()), len(marker_order)), dtype=np.float32)
    for column, marker in enumerate(marker_order):
        rows = symbol_to_rows.get(marker, [])
        if not rows:
            continue
        values = np.asarray(matrix[rows, :][:, keep].sum(axis=0)).ravel()
        result[:, column] = np.log1p(
            values.astype(np.float64) / np.maximum(retained_library, 1.0) * 10_000.0
        ).astype(np.float32)
    return result


def main() -> None:
    repo = Path(__file__).resolve().parents[1]
    source = repo / "data" / "raw" / "GSE290642" / "human_matrices"
    marker_path = repo / "config" / "gse290642_identity_markers.yaml"
    donor_path = repo / "metadata" / "donor_manifest.csv"
    marker_panels: dict[str, list[str]] = yaml.safe_load(marker_path.read_text())
    marker_order = list(dict.fromkeys(sum(marker_panels.values(), [])))

    donors = pd.read_csv(donor_path, dtype=str).fillna("")
    donors = donors[donors["dataset_id"] == "GSE290642_human"].copy()
    donors_by_gsm = donors.set_index("sample_id")
    if len(donors_by_gsm) != 24:
        raise RuntimeError(f"Expected 24 human donors, found {len(donors_by_gsm)}")

    sample_dirs = sorted(path for path in source.iterdir() if path.is_dir())
    if len(sample_dirs) != 24:
        raise RuntimeError(f"Expected 24 extracted human samples, found {len(sample_dirs)}")

    marker_blocks: list[np.ndarray] = []
    retained_frames: list[pd.DataFrame] = []
    qc_rows: list[dict[str, object]] = []
    input_hashes: dict[str, dict[str, str]] = {}
    gene_counts: dict[str, int] = {}
    available_marker_sets: list[set[str]] = []
    analysis_offset = 0

    for sample_number, sample_dir in enumerate(sample_dirs, start=1):
        gsm = sample_dir.name.split("_", 1)[0]
        if gsm not in donors_by_gsm.index:
            raise RuntimeError(f"No donor-manifest entry for {gsm}")
        matrix_path = locate_one(sample_dir, ("matrix.mtx", "matrix.mtx.gz"))
        gene_path = locate_one(
            sample_dir,
            ("genes.tsv", "genes.tsv.gz", "features.tsv", "features.tsv.gz"),
        )
        barcode_path = locate_one(sample_dir, ("barcodes.tsv", "barcodes.tsv.gz"))
        gene_ids, symbols = read_genes(gene_path)
        barcodes = pd.read_csv(barcode_path, sep="\t", header=None, dtype=str).iloc[:, 0]

        gene_counts[gsm] = len(symbols)
        available_marker_sets.append(set(symbols).intersection(marker_order))

        matrix = mmread(matrix_path).tocsr()
        if matrix.shape != (len(symbols), len(barcodes)):
            raise RuntimeError(
                f"Matrix dimensions {matrix.shape} disagree with features/barcodes for {gsm}"
            )
        library = np.asarray(matrix.sum(axis=0)).ravel().astype(np.float64)
        detected = np.asarray(matrix.getnnz(axis=0)).ravel().astype(np.int32)
        mt_rows = np.fromiter(
            (index for index, symbol in enumerate(symbols) if symbol.upper().startswith("MT-")),
            dtype=np.int64,
        )
        mt_counts = np.asarray(matrix[mt_rows, :].sum(axis=0)).ravel().astype(np.float64)
        mt_fraction = mt_counts / np.maximum(library, 1.0)
        keep = (
            (library >= MIN_UMI)
            & (detected >= MIN_GENES)
            & (detected <= MAX_GENES)
            & (mt_fraction <= MAX_MT_FRACTION)
        )
        retained_library = library[keep]
        block = marker_expression(
            matrix, symbols, marker_order, keep, retained_library
        )
        marker_blocks.append(block)

        donor = donors_by_gsm.loc[gsm]
        retained_n = int(keep.sum())
        retained = pd.DataFrame(
            {
                "analysis_row": np.arange(
                    analysis_offset, analysis_offset + retained_n, dtype=np.int64
                ),
                "sample_id": gsm,
                "donor_id": donor["donor_id"],
                "disease_group": donor["disease_group"],
                "fibrosis_stage": donor["fibrosis_stage"],
                "barcode": barcodes.to_numpy()[keep],
                "umi": library[keep].astype(np.int64),
                "detected_genes": detected[keep],
                "mt_fraction": mt_fraction[keep],
            }
        )
        retained_frames.append(retained)
        analysis_offset += retained_n

        qc_rows.append(
            {
                "sample_id": gsm,
                "donor_id": donor["donor_id"],
                "disease_group": donor["disease_group"],
                "fibrosis_stage": donor["fibrosis_stage"],
                "cells_total": len(barcodes),
                "cells_retained": retained_n,
                "retained_fraction": retained_n / len(barcodes),
                "cells_fail_low_umi": int((library < MIN_UMI).sum()),
                "cells_fail_low_genes": int((detected < MIN_GENES).sum()),
                "cells_fail_high_genes": int((detected > MAX_GENES).sum()),
                "cells_fail_mt": int((mt_fraction > MAX_MT_FRACTION).sum()),
                "median_umi_retained": float(np.median(retained_library)),
                "median_genes_retained": float(np.median(detected[keep])),
                "median_mt_fraction_retained": float(np.median(mt_fraction[keep])),
            }
        )
        input_hashes[gsm] = {
            "matrix_sha256": sha256(matrix_path),
            "genes_sha256": sha256(gene_path),
            "barcodes_sha256": sha256(barcode_path),
        }
        print(
            f"[{sample_number:02d}/24] {gsm}: {retained_n:,}/{len(barcodes):,} cells retained",
            flush=True,
        )
        del matrix, block

    available_markers = sorted(set.intersection(*available_marker_sets))
    missing_markers = sorted(set(marker_order).difference(available_markers))
    if len(available_markers) < 0.9 * len(marker_order):
        raise RuntimeError(
            f"Too many identity markers are absent: {len(missing_markers)}/{len(marker_order)}"
        )

    marker_values = np.vstack(marker_blocks)
    retained_cells = pd.concat(retained_frames, ignore_index=True)
    if len(retained_cells) != marker_values.shape[0]:
        raise RuntimeError("Retained cell metadata and marker matrix are misaligned")

    scaled = StandardScaler().fit_transform(marker_values).astype(np.float32)
    n_pcs = min(N_PCS, scaled.shape[1], scaled.shape[0] - 1)
    pca = PCA(n_components=n_pcs, svd_solver="randomized", random_state=SEED)
    coordinates = pca.fit_transform(scaled).astype(np.float32)
    clustering = MiniBatchKMeans(
        n_clusters=N_CLUSTERS,
        random_state=SEED,
        n_init=20,
        batch_size=8192,
        max_no_improvement=30,
    )
    retained_cells["cluster"] = clustering.fit_predict(coordinates)

    marker_frame = pd.DataFrame(marker_values, columns=marker_order)
    marker_frame["cluster"] = retained_cells["cluster"].to_numpy()
    cluster_means = marker_frame.groupby("cluster", sort=True)[marker_order].mean()
    marker_sd = cluster_means.std(axis=0, ddof=0).replace(0.0, np.nan)
    marker_z = ((cluster_means - cluster_means.mean(axis=0)) / marker_sd).fillna(0.0)

    audit_rows: list[dict[str, object]] = []
    for cluster in cluster_means.index:
        cell_subset = retained_cells[retained_cells["cluster"] == cluster]
        lineage_scores: dict[str, float] = {}
        lineage_anchors: dict[str, int] = {}
        for lineage, markers in marker_panels.items():
            available = [marker for marker in markers if marker in marker_z.columns]
            lineage_scores[lineage] = float(marker_z.loc[cluster, available].mean())
            lineage_anchors[lineage] = int(
                (marker_z.loc[cluster, available] > MIN_ANCHOR_Z).sum()
            )
        ordered = sorted(lineage_scores.items(), key=lambda item: item[1], reverse=True)
        winner, top_score = ordered[0]
        second, second_score = ordered[1]
        margin = top_score - second_score
        passes = margin >= MIN_SCORE_MARGIN and lineage_anchors[winner] >= MIN_ANCHORS
        row: dict[str, object] = {
            "cluster": int(cluster),
            "n_cells": len(cell_subset),
            "n_donors": int(cell_subset["donor_id"].nunique()),
            "n_control_donors": int(
                cell_subset.loc[cell_subset["fibrosis_stage"] == "F0", "donor_id"].nunique()
            ),
            "n_f4_donors": int(
                cell_subset.loc[cell_subset["fibrosis_stage"] == "F4", "donor_id"].nunique()
            ),
            "winner": winner,
            "second": second,
            "top_score": top_score,
            "second_score": second_score,
            "top_minus_second_score": margin,
            "winner_anchors_z_gt_0_5": lineage_anchors[winner],
            "passes_automatic_rule": passes,
            "provisional_label": winner if passes else "ambiguous",
        }
        for lineage in marker_panels:
            row[f"score_{lineage}"] = lineage_scores[lineage]
            row[f"anchors_{lineage}"] = lineage_anchors[lineage]
        for marker in marker_order:
            row[f"mean_log1p_cp10k_{marker}"] = float(cluster_means.loc[cluster, marker])
            row[f"z_{marker}"] = float(marker_z.loc[cluster, marker])
        audit_rows.append(row)

    audit = pd.DataFrame(audit_rows).sort_values("cluster")
    label_map = audit.set_index("cluster")["provisional_label"]
    retained_cells["provisional_label"] = retained_cells["cluster"].map(label_map)

    qc_dir = repo / "results" / "qc"
    log_dir = repo / "results" / "logs"
    interim_dir = repo / "data" / "interim" / "GSE290642"
    for directory in (qc_dir, log_dir, interim_dir):
        directory.mkdir(parents=True, exist_ok=True)
    qc_table = pd.DataFrame(qc_rows)
    qc_table.to_csv(qc_dir / "gse290642_sample_cell_qc.csv", index=False)
    audit.to_csv(qc_dir / "gse290642_cluster_marker_audit.csv", index=False)
    retained_cells.to_csv(
        interim_dir / "retained_cell_assignments.csv.gz", index=False, compression="gzip"
    )
    pd.DataFrame(marker_values, columns=marker_order).to_csv(
        interim_dir / "retained_cell_marker_log1p_cp10k.csv.gz",
        index=False,
        compression="gzip",
    )

    summary = {
        "dataset": "GSE290642_human",
        "annotation_role": "reconstructed broad-lineage sensitivity only",
        "human_donors": len(sample_dirs),
        "cells_total": int(qc_table["cells_total"].sum()),
        "cells_retained": int(qc_table["cells_retained"].sum()),
        "retained_fraction": float(
            qc_table["cells_retained"].sum() / qc_table["cells_total"].sum()
        ),
        "genes_per_sample": gene_counts,
        "gene_count_min": min(gene_counts.values()),
        "gene_count_max": max(gene_counts.values()),
        "marker_panel_total": len(marker_order),
        "marker_panel_available": len(available_markers),
        "missing_markers": missing_markers,
        "clusters": N_CLUSTERS,
        "clusters_passing_automatic_rule": int(audit["passes_automatic_rule"].sum()),
        "provisional_label_counts": audit["provisional_label"].value_counts().to_dict(),
        "qc_thresholds": {
            "min_umi": MIN_UMI,
            "min_genes": MIN_GENES,
            "max_genes": MAX_GENES,
            "max_mt_fraction": MAX_MT_FRACTION,
        },
        "annotation_parameters": {
            "normalization": "log1p(count / cell_library * 10000), identity markers only",
            "scaling": "marker-wise z score across retained cells",
            "pca_components": n_pcs,
            "pca_explained_variance_fraction": float(pca.explained_variance_ratio_.sum()),
            "clustering": "MiniBatchKMeans",
            "n_clusters": N_CLUSTERS,
            "n_init": 20,
            "seed": SEED,
            "minimum_top_minus_second_score": MIN_SCORE_MARGIN,
            "minimum_anchor_markers": MIN_ANCHORS,
            "anchor_z_threshold": MIN_ANCHOR_Z,
        },
        "marker_config_sha256": sha256(marker_path),
        "donor_manifest_sha256": sha256(donor_path),
        "input_hashes": input_hashes,
        "mapping_status": "provisional automatic labels; manual audit and frozen mapping required",
    }
    (qc_dir / "gse290642_ingest_summary.json").write_text(
        json.dumps(summary, indent=2) + "\n", encoding="utf-8"
    )
    (log_dir / "gse290642_annotation_run.json").write_text(
        json.dumps(summary, indent=2) + "\n", encoding="utf-8"
    )
    print(audit[["cluster", "n_cells", "n_donors", "winner", "top_minus_second_score", "winner_anchors_z_gt_0_5", "provisional_label"]].to_string(index=False))
    print(json.dumps({key: summary[key] for key in ["cells_total", "cells_retained", "retained_fraction", "provisional_label_counts"]}, indent=2))


if __name__ == "__main__":
    main()
