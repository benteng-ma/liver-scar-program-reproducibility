from __future__ import annotations

import json
from pathlib import Path

import numpy as np
import pandas as pd
import yaml
from scipy.io import mmread
from sklearn.cluster import MiniBatchKMeans
from sklearn.decomposition import PCA
from sklearn.preprocessing import StandardScaler

from audit_gse290642_lineages import locate_one, marker_expression, read_genes, sha256


SEED = 20260830
MIN_UMI = 500
MIN_GENES = 200
MAX_GENES = 6000
MAX_MT_FRACTION = 0.20
N_CLUSTERS = 25
N_PCS = 20
MIN_SCORE_MARGIN = 0.5
MIN_ANCHORS = 2
MIN_ANCHOR_Z = 0.5


def main() -> None:
    repo = Path(__file__).resolve().parents[1]
    source = repo / "data" / "raw" / "GSE181483" / "human_matrices"
    marker_path = repo / "config" / "gse181483_identity_markers.yaml"
    donor_path = repo / "metadata" / "donor_manifest.csv"
    marker_panels: dict[str, list[str]] = yaml.safe_load(marker_path.read_text())
    marker_order = list(dict.fromkeys(sum(marker_panels.values(), [])))

    donors = pd.read_csv(donor_path, dtype=str).fillna("")
    donors = donors[donors["dataset_id"] == "GSE181483_human"].copy()
    donors_by_gsm = donors.set_index("sample_id")
    if len(donors_by_gsm) != 4:
        raise RuntimeError(f"Expected four human donors, found {len(donors_by_gsm)}")
    sample_dirs = sorted(path for path in source.iterdir() if path.is_dir())
    if len(sample_dirs) != 4:
        raise RuntimeError(f"Expected four extracted human samples, found {len(sample_dirs)}")

    marker_blocks: list[np.ndarray] = []
    retained_frames: list[pd.DataFrame] = []
    qc_rows: list[dict[str, object]] = []
    input_hashes: dict[str, dict[str, str]] = {}
    feature_sets: list[set[str]] = []
    analysis_offset = 0

    for sample_number, sample_dir in enumerate(sample_dirs, start=1):
        gsm = sample_dir.name.split("_", 1)[0]
        matrix_path = locate_one(sample_dir, ("matrix.mtx", "matrix.mtx.gz"))
        gene_path = locate_one(sample_dir, ("genes.tsv", "genes.tsv.gz", "features.tsv", "features.tsv.gz"))
        barcode_path = locate_one(sample_dir, ("barcodes.tsv", "barcodes.tsv.gz"))
        _, symbols_original = read_genes(gene_path)
        symbols = [symbol.upper() for symbol in symbols_original]
        barcodes = pd.read_csv(barcode_path, sep="\t", header=None, dtype=str).iloc[:, 0]
        feature_sets.append({symbol for symbol in symbols if symbol and symbol != "NAN"})
        matrix = mmread(matrix_path).tocsr()
        if matrix.shape != (len(symbols), len(barcodes)):
            raise RuntimeError(f"Matrix dimensions disagree for {gsm}")
        library = np.asarray(matrix.sum(axis=0)).ravel().astype(np.float64)
        detected = np.asarray(matrix.getnnz(axis=0)).ravel().astype(np.int32)
        mt_rows = np.fromiter((i for i, symbol in enumerate(symbols) if symbol.startswith("MT-")), dtype=np.int64)
        mt_counts = np.asarray(matrix[mt_rows, :].sum(axis=0)).ravel().astype(np.float64)
        mt_fraction = mt_counts / np.maximum(library, 1.0)
        keep = (
            (library >= MIN_UMI)
            & (detected >= MIN_GENES)
            & (detected <= MAX_GENES)
            & (mt_fraction <= MAX_MT_FRACTION)
        )
        retained_library = library[keep]
        block = marker_expression(matrix, symbols, marker_order, keep, retained_library)
        marker_blocks.append(block)
        donor = donors_by_gsm.loc[gsm]
        retained_n = int(keep.sum())
        retained_frames.append(
            pd.DataFrame(
                {
                    "analysis_row": np.arange(analysis_offset, analysis_offset + retained_n, dtype=np.int64),
                    "sample_id": gsm,
                    "donor_id": donor["donor_id"],
                    "disease_group": donor["disease_group"],
                    "barcode": barcodes.to_numpy()[keep],
                    "umi": library[keep].astype(np.int64),
                    "detected_genes": detected[keep],
                    "mt_fraction": mt_fraction[keep],
                }
            )
        )
        analysis_offset += retained_n
        qc_rows.append(
            {
                "sample_id": gsm,
                "donor_id": donor["donor_id"],
                "disease_group": donor["disease_group"],
                "cells_total": len(barcodes),
                "cells_retained": retained_n,
                "retained_fraction": retained_n / len(barcodes),
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
        print(f"[{sample_number}/4] {gsm}: {retained_n:,}/{len(barcodes):,} cells retained", flush=True)

    marker_values = np.vstack(marker_blocks)
    retained_cells = pd.concat(retained_frames, ignore_index=True)
    if len(retained_cells) != marker_values.shape[0]:
        raise RuntimeError("Retained metadata and marker matrix are misaligned")
    available_markers = set.intersection(*({m for m in marker_order if m in features} for features in feature_sets))
    if len(available_markers) < 0.9 * len(marker_order):
        raise RuntimeError("Identity marker coverage is below 90%")

    scaled = StandardScaler().fit_transform(marker_values).astype(np.float32)
    n_pcs = min(N_PCS, scaled.shape[1], scaled.shape[0] - 1)
    pca = PCA(n_components=n_pcs, svd_solver="randomized", random_state=SEED)
    coordinates = pca.fit_transform(scaled).astype(np.float32)
    clustering = MiniBatchKMeans(
        n_clusters=N_CLUSTERS,
        random_state=SEED,
        n_init=20,
        batch_size=4096,
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
        subset = retained_cells[retained_cells["cluster"] == cluster]
        scores: dict[str, float] = {}
        anchors: dict[str, int] = {}
        for lineage, markers in marker_panels.items():
            scores[lineage] = float(marker_z.loc[cluster, markers].mean())
            anchors[lineage] = int((marker_z.loc[cluster, markers] > MIN_ANCHOR_Z).sum())
        ordered = sorted(scores.items(), key=lambda item: item[1], reverse=True)
        winner, top_score = ordered[0]
        second, second_score = ordered[1]
        margin = top_score - second_score
        passes = margin >= MIN_SCORE_MARGIN and anchors[winner] >= MIN_ANCHORS
        row: dict[str, object] = {
            "cluster": int(cluster),
            "n_cells": len(subset),
            "n_donors": int(subset["donor_id"].nunique()),
            "n_control_donors": int(subset.loc[subset["disease_group"] == "healthy", "donor_id"].nunique()),
            "n_cirrhosis_donors": int(subset.loc[subset["disease_group"] == "cirrhosis", "donor_id"].nunique()),
            "winner": winner,
            "second": second,
            "top_score": top_score,
            "second_score": second_score,
            "top_minus_second_score": margin,
            "winner_anchors_z_gt_0_5": anchors[winner],
            "passes_automatic_rule": passes,
            "provisional_label": winner if passes else "ambiguous",
        }
        for lineage in marker_panels:
            row[f"score_{lineage}"] = scores[lineage]
            row[f"anchors_{lineage}"] = anchors[lineage]
        for marker in marker_order:
            row[f"mean_log1p_cp10k_{marker}"] = float(cluster_means.loc[cluster, marker])
            row[f"z_{marker}"] = float(marker_z.loc[cluster, marker])
        audit_rows.append(row)

    audit = pd.DataFrame(audit_rows).sort_values("cluster")
    retained_cells["provisional_label"] = retained_cells["cluster"].map(audit.set_index("cluster")["provisional_label"])
    qc_dir = repo / "results" / "qc"
    log_dir = repo / "results" / "logs"
    interim_dir = repo / "data" / "interim" / "GSE181483"
    for directory in (qc_dir, log_dir, interim_dir):
        directory.mkdir(parents=True, exist_ok=True)
    pd.DataFrame(qc_rows).to_csv(qc_dir / "gse181483_sample_cell_qc.csv", index=False)
    audit.to_csv(qc_dir / "gse181483_cluster_marker_audit.csv", index=False)
    retained_cells.to_csv(interim_dir / "retained_cell_assignments.csv.gz", index=False, compression="gzip")
    summary = {
        "dataset": "GSE181483_human",
        "annotation_role": "reconstructed broad-lineage directional support only",
        "human_donors": 4,
        "cells_total": int(sum(row["cells_total"] for row in qc_rows)),
        "cells_retained": len(retained_cells),
        "clusters": N_CLUSTERS,
        "marker_panel_total": len(marker_order),
        "marker_panel_available_all_donors": len(available_markers),
        "clusters_passing_automatic_rule": int(audit["passes_automatic_rule"].sum()),
        "provisional_label_counts": audit["provisional_label"].value_counts().to_dict(),
        "qc_thresholds": {"min_umi": MIN_UMI, "min_genes": MIN_GENES, "max_genes": MAX_GENES, "max_mt_fraction": MAX_MT_FRACTION},
        "annotation_parameters": {
            "normalization": "log1p(count / cell_library * 10000), identity markers only",
            "pca_components": n_pcs,
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
        "mapping_status": "provisional automatic labels; frozen mapping required before aggregation",
    }
    (qc_dir / "gse181483_ingest_summary.json").write_text(json.dumps(summary, indent=2) + "\n", encoding="utf-8")
    (log_dir / "gse181483_annotation_run.json").write_text(json.dumps(summary, indent=2) + "\n", encoding="utf-8")
    print(audit[["cluster", "n_cells", "n_donors", "winner", "top_minus_second_score", "winner_anchors_z_gt_0_5", "provisional_label"]].to_string(index=False))
    print(json.dumps(summary, indent=2))


if __name__ == "__main__":
    main()
