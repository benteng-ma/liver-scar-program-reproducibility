from __future__ import annotations

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
    assignments["source_state"] = "cluster_" + assignments["cluster"].astype(str)

    columns = []
    rows: list[dict[str, object]] = []
    reference_genes: list[str] | None = None
    group_number = 0
    for donor in donors.itertuples(index=False):
        matrix, genes, barcodes = read_10x_h5(h5_dir / donor.h5_file)
        if reference_genes is None:
            reference_genes = genes
        elif genes != reference_genes:
            raise RuntimeError(f"feature order differs for {donor.sample_id}")
        barcode_to_column = {barcode: index for index, barcode in enumerate(barcodes)}
        donor_rows = assignments[assignments["sample_id"].eq(donor.sample_id)]
        for (lineage, source_state), selected in donor_rows.groupby(
            ["harmonized_lineage", "source_state"], sort=True
        ):
            indices = np.array(
                [barcode_to_column[value] for value in selected["barcode"]], dtype=np.int64
            )
            aggregate = csc_matrix(np.asarray(matrix[:, indices].sum(axis=1), dtype=np.int64))
            columns.append(aggregate)
            group_number += 1
            rows.append(
                {
                    "state_group_id": f"S{group_number:05d}",
                    "donor_id": donor.donor_id,
                    "sample_id": donor.sample_id,
                    "harmonized_lineage": lineage,
                    "source_state": source_state,
                    "n_cells": len(indices),
                    "library_size": int(aggregate.sum()),
                    "disease_group": donor.disease_group,
                    "etiology": donor.etiology,
                    "fibrosis_stage": donor.fibrosis_stage,
                    "metabolic_order": donor.metabolic_order,
                    "age": donor.age,
                    "sex": donor.sex,
                }
            )
        del matrix

    if reference_genes is None:
        raise RuntimeError("no genes recovered")
    aggregated = hstack(columns, format="csc")
    manifest = pd.DataFrame(rows)
    if aggregated.shape[0] != len(reference_genes) or aggregated.shape[1] != len(manifest):
        raise RuntimeError("state matrix and manifest mismatch")
    mmwrite(interim / "phase5_donor_state_raw_counts.mtx", aggregated)
    pd.DataFrame({"gene": reference_genes}).to_csv(
        interim / "phase5_state_genes.csv", index=False
    )
    manifest.to_csv(interim / "phase5_donor_state_manifest.csv", index=False)
    print(
        manifest.groupby(["harmonized_lineage", "source_state", "disease_group"])["n_cells"]
        .agg(["count", "min", "median", "max"])
        .to_string()
    )


if __name__ == "__main__":
    main()
