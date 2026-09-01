from __future__ import annotations

import hashlib
import json
from pathlib import Path

import anndata as ad
import numpy as np
import pandas as pd
from scipy.io import mmwrite
from scipy.sparse import csc_matrix


LABEL_MAP = {
    ("Normal", "LSEC"): "endothelial",
    ("Normal", "HSC_1"): "mesenchymal_hsc_myofibroblast",
    ("Normal", "HSC_2"): "mesenchymal_hsc_myofibroblast",
    ("Normal", "Mac_1"): "macrophage_monocyte",
    ("Normal", "Mac_2"): "macrophage_monocyte",
    ("Disease", "LSEC"): "endothelial",
    ("Disease", "HSC"): "mesenchymal_hsc_myofibroblast",
    ("Disease", "Mac"): "macrophage_monocyte",
}
TARGET_LINEAGES = [
    "endothelial",
    "macrophage_monocyte",
    "mesenchymal_hsc_myofibroblast",
]


def sha256(path: Path) -> str:
    digest = hashlib.sha256()
    with path.open("rb") as handle:
        for block in iter(lambda: handle.read(8 * 1024 * 1024), b""):
            digest.update(block)
    return digest.hexdigest().upper()


def main() -> None:
    repo = Path(__file__).resolve().parents[1]
    source = repo / "data" / "raw" / "GSE210077_Watson6" / "processed" / "adata_healthy_diseased_nucseq.h5ad"
    output = repo / "data" / "interim" / "GSE210077_Watson6"
    output.mkdir(parents=True, exist_ok=True)
    object_data = ad.read_h5ad(source, backed="r")
    obs = object_data.obs.copy()
    if object_data.n_obs != 27_692 or obs["sample_id"].nunique() != 6:
        raise RuntimeError("Watson processed object dimensions/donors changed")
    if object_data.raw is None:
        raise RuntimeError("author object lacks raw expression layer")
    all_genes = pd.Index(object_data.raw.var_names.astype(str).str.upper())
    if all_genes.has_duplicates:
        raise RuntimeError("raw gene symbols are duplicated")

    obs["author_final_label"] = np.where(
        obs["Condition"].astype(str).eq("Normal"),
        obs["cell_type_final_healthy"].astype("string"),
        obs["cell_type_final_injured"].astype("string"),
    )
    obs["harmonized_lineage"] = [
        LABEL_MAP.get((str(condition), str(label)), np.nan)
        for condition, label in zip(obs["Condition"], obs["author_final_label"])
    ]
    target = obs["harmonized_lineage"].notna()

    donor_manifest = pd.read_csv(repo / "metadata" / "donor_manifest.csv", dtype=str).fillna("")
    donor_manifest = donor_manifest[donor_manifest["dataset_id"] == "GSE210077_Watson6"].set_index("donor_id")
    if set(obs["sample_id"].astype(str)) != set(donor_manifest.index):
        raise RuntimeError("author object donor IDs disagree with frozen Watson subset")

    feature_rows: list[dict[str, object]] = []
    availability_masks: list[np.ndarray] = []
    for donor_id in donor_manifest.index:
        donor_indices = np.flatnonzero(obs["sample_id"].astype(str).eq(donor_id).to_numpy())
        probe = np.asarray(object_data.raw.X[donor_indices[: min(5, len(donor_indices))], :])
        nan_pattern = np.isnan(probe)
        if not np.all(nan_pattern == nan_pattern[0, :]):
            raise RuntimeError(f"raw feature availability varies between cells for {donor_id}")
        available = ~nan_pattern[0, :]
        availability_masks.append(available)
        feature_rows.append(
            {
                "donor_id": donor_id,
                "raw_features_available": int(available.sum()),
                "raw_features_nan_filled": int((~available).sum()),
            }
        )
    shared_mask = np.logical_and.reduce(availability_masks)
    genes = all_genes[shared_mask]
    if len(genes) < 20_000:
        raise RuntimeError(f"unexpectedly small shared feature space: {len(genes)}")

    columns: list[np.ndarray] = []
    manifest_rows: list[dict[str, object]] = []
    for donor_id in donor_manifest.index:
        donor = donor_manifest.loc[donor_id]
        for lineage in TARGET_LINEAGES:
            selected = target & obs["sample_id"].astype(str).eq(donor_id) & obs["harmonized_lineage"].eq(lineage)
            indices = np.flatnonzero(selected.to_numpy())
            raw_values = np.asarray(object_data.raw.X[indices, :])[:, shared_mask]
            if np.isnan(raw_values).any():
                raise RuntimeError("shared raw feature space still contains NaN values")
            if raw_values.size and not np.allclose(raw_values, np.round(raw_values)):
                raise RuntimeError("author raw layer contains non-integer expression")
            counts = np.rint(raw_values.sum(axis=0)).astype(np.int64)
            columns.append(counts)
            manifest_rows.append(
                {
                    "group_id": f"G{len(manifest_rows) + 1:04d}",
                    "donor_id": donor_id,
                    "sample_id": donor["sample_id"],
                    "harmonized_lineage": lineage,
                    "n_cells": len(indices),
                    "disease_group": donor["disease_group"],
                    "etiology": donor["etiology"],
                    "fibrosis_stage": donor["fibrosis_stage"],
                    "age": donor["age"],
                    "sex": donor["sex"],
                    "author_condition": str(obs.loc[selected, "Condition"].iloc[0]) if len(indices) else "",
                }
            )
            print(f"{donor_id} {lineage}: {len(indices):,} nuclei", flush=True)
    object_data.file.close()

    aggregated = np.column_stack(columns)
    manifest = pd.DataFrame(manifest_rows)
    if aggregated.shape != (len(genes), 18) or len(manifest) != 18:
        raise RuntimeError("expected 6 donors x 3 lineages")
    if (manifest["n_cells"] == 0).any():
        raise RuntimeError("a Watson donor-lineage group has zero author-labelled nuclei")
    mmwrite(output / "donor_lineage_raw_counts.mtx", csc_matrix(aggregated))
    pd.DataFrame({"gene": genes}).to_csv(output / "genes.csv", index=False)
    pd.DataFrame(feature_rows).to_csv(output / "donor_raw_feature_availability.csv", index=False)
    manifest.to_csv(output / "donor_lineage_manifest.csv", index=False)
    (
        obs.groupby(
            ["sample_id", "Condition", "author_final_label", "harmonized_lineage"],
            dropna=False,
            observed=True,
        )
        .size()
        .reset_index(name="n_cells")
        .to_csv(output / "donor_author_annotation_cell_counts.csv", index=False)
    )
    summary = {
        "dataset_id": "GSE210077_Watson6",
        "author_object_cells": 27_692,
        "author_object_analysis_genes": 24_619,
        "author_object_raw_genes_union": len(all_genes),
        "author_object_raw_genes_shared_all_six_donors": len(genes),
        "donors": 6,
        "target_cells": int(manifest["n_cells"].sum()),
        "target_donor_lineage_groups": len(manifest),
        "source_bytes": source.stat().st_size,
        "source_sha256": sha256(source),
        "source_archive_sha256": "58AD9E85EABAB28E41ABB825C816C06D1977A1A873BC6501AB5B99FC90C519A7",
        "raw_layer_validation": "integer-valued adata.raw.X after excluding author NaN-filled genes absent from a donor feature space",
        "author_label_mapping": {
            "healthy": {"LSEC": "endothelial", "HSC_1+HSC_2": "mesenchymal_hsc_myofibroblast", "Mac_1+Mac_2": "macrophage_monocyte"},
            "fibrotic": {"LSEC": "endothelial", "HSC": "mesenchymal_hsc_myofibroblast", "Mac": "macrophage_monocyte"},
        },
        "endpoint_status": "three healthy versus F2/F3/F4 fibrosis; small directional/sensitivity only",
    }
    (output / "donor_aggregation_summary.json").write_text(json.dumps(summary, indent=2) + "\n", encoding="utf-8")
    print(manifest.groupby(["harmonized_lineage", "disease_group"])["n_cells"].agg(["count", "min", "median", "max"]).to_string())
    print(json.dumps(summary, indent=2))


if __name__ == "__main__":
    main()
