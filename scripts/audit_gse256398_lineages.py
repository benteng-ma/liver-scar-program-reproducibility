from __future__ import annotations

import hashlib
import json
import re
import xml.etree.ElementTree as ET
from pathlib import Path

import h5py
import numpy as np
import pandas as pd
import yaml
from scipy.sparse import csc_matrix
from sklearn.cluster import MiniBatchKMeans
from sklearn.decomposition import PCA
from sklearn.preprocessing import StandardScaler


SEED = 20260831
MIN_UMI = 500
MIN_GENES = 200
MAX_GENES = 6000
MAX_MT_FRACTION = 0.20
N_CLUSTERS = 40
N_PCS = 20
MIN_SCORE_MARGIN = 0.5
MIN_ANCHORS = 2
MIN_ANCHOR_Z = 0.5
NS = {"m": "http://www.ncbi.nlm.nih.gov/geo/info/MINiML"}


def sha256(path: Path) -> str:
    digest = hashlib.sha256()
    with path.open("rb") as handle:
        for block in iter(lambda: handle.read(8 * 1024 * 1024), b""):
            digest.update(block)
    return digest.hexdigest().upper()


def decode(values) -> list[str]:
    return [value.decode("utf-8") if isinstance(value, bytes) else str(value) for value in values]


def read_10x_h5(path: Path):
    with h5py.File(path, "r") as handle:
        group = handle["matrix"]
        shape = tuple(int(value) for value in group["shape"][:])
        matrix = csc_matrix(
            (
                group["data"][:],
                group["indices"][:],
                group["indptr"][:],
            ),
            shape=shape,
        )
        symbols = decode(group["features"]["name"][:])
        barcodes = decode(group["barcodes"][:])
    if matrix.shape != (len(symbols), len(barcodes)):
        raise RuntimeError(f"Matrix dimensions disagree with features/barcodes: {path}")
    return matrix, symbols, barcodes


def classify_group(title: str) -> tuple[str, str, str, int | None]:
    if "Alcohol-associated Cirrhosis" in title:
        return "alcohol_cirrhosis", "alcohol", "cirrhosis", None
    if "Alcohol-associated Hepatitis" in title:
        return "alcohol_hepatitis", "alcohol", "hepatitis", None
    if "Healthy Control" in title:
        # Healthy controls remain available for healthy-versus-disease contrasts,
        # but are not part of the prespecified metabolic-only ordinal trajectory.
        return "healthy", "none", "F0", None
    if "MASH Cirrhosis" in title:
        return "mash_cirrhosis", "metabolic", "F4", 2
    if "MASLD F0" in title:
        return "masld_f0", "metabolic", "F0", 0
    if "MASH Fibrosis F2-3" in title:
        return "mash_fibrosis", "metabolic", "F2-3", 1
    if "MASH Fibrosis F3" in title:
        return "mash_fibrosis", "metabolic", "F3", 1
    raise RuntimeError(f"Unrecognized frozen human sample title: {title}")


def sample_manifest(xml_path: Path, h5_dir: Path) -> pd.DataFrame:
    root = ET.parse(xml_path).getroot()
    h5_by_gsm = {path.name.split("_", 1)[0]: path for path in sorted(h5_dir.glob("*.h5"))}
    rows: list[dict[str, object]] = []
    for sample in root.findall("m:Sample", NS):
        gsm = sample.attrib["iid"]
        title = (sample.findtext("m:Title", default="", namespaces=NS) or "").strip()
        if ", Human," not in title:
            continue
        if gsm not in h5_by_gsm:
            raise RuntimeError(f"Missing H5 for {gsm}")
        characteristics: dict[str, str] = {}
        for node in sample.findall("m:Channel/m:Characteristics", NS):
            characteristics[(node.attrib.get("tag") or "").strip().lower()] = (node.text or "").strip()
        disease_group, etiology, fibrosis_stage, metabolic_order = classify_group(title)
        sample_label = title.split(",", 1)[0].strip()
        age_match = re.search(r"(\d+)", characteristics.get("age", ""))
        rows.append(
            {
                "dataset_id": "GSE256398_human",
                "sample_id": gsm,
                "donor_id": sample_label,
                "title": title,
                "disease_group": disease_group,
                "etiology": etiology,
                "fibrosis_stage": fibrosis_stage,
                "metabolic_order": metabolic_order,
                "age": int(age_match.group(1)) if age_match else np.nan,
                "sex": characteristics.get("sex", "").lower(),
                "h5_file": h5_by_gsm[gsm].name,
            }
        )
    result = pd.DataFrame(rows).sort_values("sample_id").reset_index(drop=True)
    if len(result) != 26 or result["donor_id"].nunique() != 26:
        raise RuntimeError(f"Expected 26 independent human donors, found {len(result)}")
    expected = {
        "healthy": 6,
        "alcohol_cirrhosis": 4,
        "alcohol_hepatitis": 5,
        "masld_f0": 3,
        "mash_fibrosis": 4,
        "mash_cirrhosis": 4,
    }
    observed = result["disease_group"].value_counts().to_dict()
    if observed != expected:
        raise RuntimeError(f"Frozen group counts disagree: {observed}")
    return result


def marker_expression(matrix, symbols, marker_order, keep, retained_library):
    symbol_to_rows: dict[str, list[int]] = {}
    marker_set = set(marker_order)
    for row, symbol in enumerate(symbols):
        symbol = symbol.upper()
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
    raw = repo / "data" / "raw" / "GSE256398"
    h5_dir = raw / "human_h5"
    xml_path = raw / "GSE256398_family.xml"
    marker_path = repo / "config" / "gse290642_identity_markers.yaml"
    marker_panels: dict[str, list[str]] = yaml.safe_load(marker_path.read_text())
    marker_panels = {key: [gene.upper() for gene in genes] for key, genes in marker_panels.items()}
    marker_order = list(dict.fromkeys(sum(marker_panels.values(), [])))
    donors = sample_manifest(xml_path, h5_dir)

    marker_blocks: list[np.ndarray] = []
    retained_frames: list[pd.DataFrame] = []
    qc_rows: list[dict[str, object]] = []
    input_hashes: dict[str, str] = {}
    available_marker_sets: list[set[str]] = []
    analysis_offset = 0

    for number, donor in enumerate(donors.itertuples(index=False), start=1):
        path = h5_dir / donor.h5_file
        matrix, symbols, barcodes = read_10x_h5(path)
        symbols_upper = [symbol.upper() for symbol in symbols]
        available_marker_sets.append(set(symbols_upper).intersection(marker_order))
        library = np.asarray(matrix.sum(axis=0)).ravel().astype(np.float64)
        detected = np.asarray(matrix.getnnz(axis=0)).ravel().astype(np.int32)
        mt_rows = np.fromiter(
            (index for index, symbol in enumerate(symbols_upper) if symbol.startswith("MT-")),
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
        block = marker_expression(matrix, symbols_upper, marker_order, keep, retained_library)
        marker_blocks.append(block)
        retained_n = int(keep.sum())
        retained_frames.append(
            pd.DataFrame(
                {
                    "analysis_row": np.arange(analysis_offset, analysis_offset + retained_n, dtype=np.int64),
                    "sample_id": donor.sample_id,
                    "donor_id": donor.donor_id,
                    "disease_group": donor.disease_group,
                    "etiology": donor.etiology,
                    "fibrosis_stage": donor.fibrosis_stage,
                    "metabolic_order": donor.metabolic_order,
                    "age": donor.age,
                    "sex": donor.sex,
                    "barcode": np.asarray(barcodes, dtype=object)[keep],
                    "umi": library[keep].astype(np.int64),
                    "detected_genes": detected[keep],
                    "mt_fraction": mt_fraction[keep],
                }
            )
        )
        analysis_offset += retained_n
        qc_rows.append(
            {
                "sample_id": donor.sample_id,
                "donor_id": donor.donor_id,
                "disease_group": donor.disease_group,
                "fibrosis_stage": donor.fibrosis_stage,
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
        input_hashes[donor.sample_id] = sha256(path)
        print(f"[{number:02d}/26] {donor.sample_id} {donor.donor_id}: {retained_n:,}/{len(barcodes):,}", flush=True)
        del matrix, block

    available_markers = sorted(set.intersection(*available_marker_sets))
    missing_markers = sorted(set(marker_order).difference(available_markers))
    if len(available_markers) < 0.9 * len(marker_order):
        raise RuntimeError(f"Too many identity markers absent: {len(missing_markers)}/{len(marker_order)}")

    marker_values = np.vstack(marker_blocks)
    retained_cells = pd.concat(retained_frames, ignore_index=True)
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
        subset = retained_cells[retained_cells["cluster"] == cluster]
        scores: dict[str, float] = {}
        anchors: dict[str, int] = {}
        for lineage, markers in marker_panels.items():
            available = [marker for marker in markers if marker in marker_z.columns]
            scores[lineage] = float(marker_z.loc[cluster, available].mean())
            anchors[lineage] = int((marker_z.loc[cluster, available] > MIN_ANCHOR_Z).sum())
        ordered = sorted(scores.items(), key=lambda item: item[1], reverse=True)
        winner, top_score = ordered[0]
        second, second_score = ordered[1]
        margin = top_score - second_score
        passes = margin >= MIN_SCORE_MARGIN and anchors[winner] >= MIN_ANCHORS
        row: dict[str, object] = {
            "cluster": int(cluster),
            "n_cells": len(subset),
            "n_donors": int(subset["donor_id"].nunique()),
            "winner": winner,
            "second": second,
            "top_score": top_score,
            "second_score": second_score,
            "top_minus_second_score": margin,
            "winner_anchors_z_gt_0_5": anchors[winner],
            "passes_automatic_rule": passes,
            "provisional_label": winner if passes else "ambiguous",
        }
        for group, count in subset["disease_group"].value_counts().items():
            row[f"donors_{group}"] = int(subset.loc[subset["disease_group"] == group, "donor_id"].nunique())
        for lineage in marker_panels:
            row[f"score_{lineage}"] = scores[lineage]
            row[f"anchors_{lineage}"] = anchors[lineage]
        for marker in marker_order:
            row[f"mean_log1p_cp10k_{marker}"] = float(cluster_means.loc[cluster, marker])
            row[f"z_{marker}"] = float(marker_z.loc[cluster, marker])
        audit_rows.append(row)

    audit = pd.DataFrame(audit_rows).sort_values("cluster")
    label_map = audit.set_index("cluster")["provisional_label"]
    retained_cells["provisional_label"] = retained_cells["cluster"].map(label_map)

    qc_dir = repo / "results" / "qc"
    log_dir = repo / "results" / "logs"
    interim = repo / "data" / "interim" / "GSE256398"
    metadata_dir = repo / "metadata"
    for directory in (qc_dir, log_dir, interim, metadata_dir):
        directory.mkdir(parents=True, exist_ok=True)
    donors.to_csv(metadata_dir / "gse256398_donor_manifest.csv", index=False)
    pd.DataFrame(qc_rows).to_csv(qc_dir / "gse256398_sample_cell_qc.csv", index=False)
    audit.to_csv(qc_dir / "gse256398_cluster_marker_audit.csv", index=False)
    retained_cells.to_csv(
        interim / "retained_cell_assignments.csv.gz",
        index=False,
        compression={"method": "gzip", "mtime": 0},
    )
    pd.DataFrame(marker_values, columns=marker_order).to_csv(
        interim / "retained_cell_marker_log1p_cp10k.csv.gz",
        index=False,
        compression={"method": "gzip", "mtime": 0},
    )
    summary = {
        "dataset": "GSE256398_human",
        "annotation_role": "post-lock reconstructed broad-lineage external validation",
        "human_donors": 26,
        "group_counts": donors["disease_group"].value_counts().to_dict(),
        "cells_total": int(sum(row["cells_total"] for row in qc_rows)),
        "cells_retained": int(sum(row["cells_retained"] for row in qc_rows)),
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
            "normalization": "log1p(count / cell_library * 10000), frozen identity markers only",
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
        "raw_tar_sha256": sha256(raw / "GSE256398_RAW.tar"),
        "miniml_sha256": sha256(raw / "GSE256398_family.xml.tgz"),
        "marker_config_sha256": sha256(marker_path),
        "input_h5_sha256": input_hashes,
        "mapping_status": "automatic provisional clusters; manual marker audit and mapping freeze required",
    }
    (qc_dir / "gse256398_ingest_summary.json").write_text(json.dumps(summary, indent=2) + "\n", encoding="utf-8")
    (log_dir / "gse256398_annotation_run.json").write_text(json.dumps(summary, indent=2) + "\n", encoding="utf-8")
    print(audit[["cluster", "n_cells", "n_donors", "winner", "top_minus_second_score", "winner_anchors_z_gt_0_5", "provisional_label"]].to_string(index=False))
    print(json.dumps({key: summary[key] for key in ["cells_total", "cells_retained", "group_counts", "provisional_label_counts"]}, indent=2))


if __name__ == "__main__":
    main()
