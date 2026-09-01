from __future__ import annotations

import json
from pathlib import Path

import numpy as np
import pandas as pd
from scipy.io import mmwrite
from scipy.sparse import csc_matrix, hstack

from audit_gse256398_lineages import read_10x_h5


TARGETS = ["endothelial", "macrophage_monocyte", "mesenchymal_hsc_myofibroblast"]


def main() -> None:
    repo = Path(__file__).resolve().parents[1]
    h5_dir = repo / "data" / "raw" / "GSE256398" / "human_h5"
    interim = repo / "data" / "interim" / "GSE256398"
    assignments = pd.read_csv(interim / "retained_cell_assignments.csv.gz")
    mapping = pd.read_csv(repo / "metadata" / "gse256398_cluster_mapping.csv")
    donors = pd.read_csv(repo / "metadata" / "gse256398_donor_manifest.csv")
    cluster_to_lineage = {
        int(row.cluster): row.harmonized_lineage
        for row in mapping.itertuples()
        if row.target_included == "yes"
    }
    assignments["harmonized_lineage"] = assignments["cluster"].map(cluster_to_lineage)
    assignments = assignments[assignments["harmonized_lineage"].isin(TARGETS)].copy()

    columns = []
    manifest_rows: list[dict[str, object]] = []
    reference_genes: list[str] | None = None
    group_number = 0
    for donor in donors.itertuples(index=False):
        matrix, genes, barcodes = read_10x_h5(h5_dir / donor.h5_file)
        if reference_genes is None:
            reference_genes = genes
        elif genes != reference_genes:
            raise RuntimeError(f"Feature order differs for {donor.sample_id}")
        barcode_to_column = {barcode: index for index, barcode in enumerate(barcodes)}
        donor_assignments = assignments[assignments["sample_id"].eq(donor.sample_id)]
        for lineage in TARGETS:
            selected = donor_assignments[donor_assignments["harmonized_lineage"].eq(lineage)]
            indices = np.array([barcode_to_column[value] for value in selected["barcode"]], dtype=np.int64)
            if not len(indices):
                raise RuntimeError(f"No retained {lineage} cells for {donor.sample_id}")
            aggregate = csc_matrix(np.asarray(matrix[:, indices].sum(axis=1), dtype=np.int64))
            columns.append(aggregate)
            group_number += 1
            manifest_rows.append(
                {
                    "group_id": f"G{group_number:04d}",
                    "donor_id": donor.donor_id,
                    "sample_id": donor.sample_id,
                    "harmonized_lineage": lineage,
                    "n_cells": len(indices),
                    "disease_group": donor.disease_group,
                    "etiology": donor.etiology,
                    "fibrosis_stage": donor.fibrosis_stage,
                    "metabolic_order": donor.metabolic_order,
                    "age": donor.age,
                    "sex": donor.sex,
                    "library_size": int(aggregate.sum()),
                }
            )
        del matrix

    if reference_genes is None:
        raise RuntimeError("No genes recovered")
    pseudobulk = hstack(columns, format="csc")
    manifest = pd.DataFrame(manifest_rows)
    if pseudobulk.shape != (len(reference_genes), 78) or len(manifest) != 78:
        raise RuntimeError(f"Expected 36,601 genes x 78 donor-lineages, found {pseudobulk.shape}")
    mmwrite(interim / "donor_lineage_raw_counts.mtx", pseudobulk)
    pd.DataFrame({"gene": reference_genes}).to_csv(interim / "genes.csv", index=False)
    manifest.to_csv(interim / "donor_lineage_manifest.csv", index=False)
    eligibility = manifest.copy()
    eligibility["eligible_30"] = eligibility["n_cells"].ge(30)
    eligibility["eligible_20"] = eligibility["n_cells"].ge(20)
    eligibility.to_csv(repo / "metadata" / "gse256398_donor_lineage_eligibility.csv", index=False)

    summary = {
        "dataset_id": "GSE256398_human",
        "genes": len(reference_genes),
        "human_donors": 26,
        "donor_lineage_groups": len(manifest),
        "target_cells": int(manifest["n_cells"].sum()),
        "group_counts": donors["disease_group"].value_counts().to_dict(),
        "eligible_30_by_group_and_lineage": [
            {
                "disease_group": disease_group,
                "lineage": lineage,
                "n_donors": int(n_donors),
            }
            for (disease_group, lineage), n_donors in eligibility.groupby(
                ["disease_group", "harmonized_lineage"]
            )["eligible_30"].sum().items()
        ],
        "mapping_file": "metadata/gse256398_cluster_mapping.csv",
        "annotation_role": "post-lock reconstructed broad-lineage external validation",
    }
    (interim / "donor_aggregation_summary.json").write_text(
        json.dumps(summary, indent=2, default=str) + "\n", encoding="utf-8"
    )
    print(
        manifest.groupby(["harmonized_lineage", "disease_group"])["n_cells"]
        .agg(["count", "min", "median", "max"])
        .to_string()
    )
    print(json.dumps(summary, indent=2, default=str))


if __name__ == "__main__":
    main()
