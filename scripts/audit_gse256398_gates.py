from __future__ import annotations

import json
from pathlib import Path

import pandas as pd

from audit_gse244832_gates import coverage_tier, read_matrix_market_dimensions


PRIMARY_CELLS = 30
SENSITIVITY_CELLS = 20
MIN_DONORS_PER_GROUP = 3


CONTRASTS = {
    "mash_cirrhosis_vs_healthy": ({"healthy"}, {"mash_cirrhosis"}),
    "alcohol_cirrhosis_vs_healthy": ({"healthy"}, {"alcohol_cirrhosis"}),
    "mash_fibrosis_vs_masld_f0": ({"masld_f0"}, {"mash_fibrosis"}),
    "mash_cirrhosis_vs_masld_f0": ({"masld_f0"}, {"mash_cirrhosis"}),
    "mash_vs_alcohol_cirrhosis_etiology": ({"alcohol_cirrhosis"}, {"mash_cirrhosis"}),
}


def main() -> None:
    repo = Path(__file__).resolve().parents[1]
    interim = repo / "data" / "interim" / "GSE256398"
    manifest = pd.read_csv(interim / "donor_lineage_manifest.csv")
    genes = pd.read_csv(interim / "genes.csv")["gene"].astype(str)
    n_rows, n_columns, n_nonzero = read_matrix_market_dimensions(interim / "donor_lineage_raw_counts.mtx")
    if (n_rows, n_columns) != (len(genes), len(manifest)):
        raise RuntimeError("Pseudobulk matrix dimensions disagree with manifests")
    if manifest["donor_id"].nunique() != 26 or len(manifest) != 78:
        raise RuntimeError("Expected 26 donors x three target lineages")
    manifest["eligible_30"] = manifest["n_cells"].ge(PRIMARY_CELLS)
    manifest["eligible_20"] = manifest["n_cells"].ge(SENSITIVITY_CELLS)
    manifest.to_csv(repo / "metadata" / "gse256398_donor_lineage_eligibility.csv", index=False)

    gate_rows: list[dict[str, object]] = []
    for contrast, (controls, cases) in CONTRASTS.items():
        scoped = manifest[manifest["disease_group"].isin(controls | cases)].copy()
        scoped["comparison_group"] = scoped["disease_group"].map(
            {**{value: "control" for value in controls}, **{value: "case" for value in cases}}
        )
        for lineage, lineage_data in scoped.groupby("harmonized_lineage"):
            rows_for_lineage = []
            for group in ("control", "case"):
                values = lineage_data[lineage_data["comparison_group"].eq(group)]
                rows_for_lineage.append(
                    {
                        "dataset_id": "GSE256398_human",
                        "contrast": contrast,
                        "lineage": lineage,
                        "comparison_group": group,
                        "donors_before_cell_gate": len(values),
                        "donors_30": int(values["eligible_30"].sum()),
                        "donors_20": int(values["eligible_20"].sum()),
                        "minimum_cells": int(values["n_cells"].min()),
                        "median_cells": float(values["n_cells"].median()),
                        "maximum_cells": int(values["n_cells"].max()),
                    }
                )
            pass_30 = all(row["donors_30"] >= MIN_DONORS_PER_GROUP for row in rows_for_lineage)
            pass_20 = all(row["donors_20"] >= MIN_DONORS_PER_GROUP for row in rows_for_lineage)
            for row in rows_for_lineage:
                row["formal_30_cell_gate"] = "PASS" if pass_30 else "FAIL"
                row["formal_20_cell_gate"] = "PASS" if pass_20 else "FAIL"
                row["analysis_role"] = "post_lock_external_validation"
                gate_rows.append(row)
    gates = pd.DataFrame(gate_rows)

    programs = pd.read_csv(repo / "literature" / "program_inventory.csv")
    measured = set(genes.str.upper())
    coverage_rows: list[dict[str, object]] = []
    for (program_id, lineage), rows in programs.groupby(["program_id", "cell_lineage"]):
        program_genes = set(rows["gene_symbol"].dropna().astype(str).str.upper())
        detected = sorted(program_genes & measured)
        missing = sorted(program_genes - measured)
        coverage = len(detected) / len(program_genes)
        coverage_rows.append(
            {
                "dataset_id": "GSE256398_human",
                "program_id": program_id,
                "lineage": lineage,
                "n_program_genes": len(program_genes),
                "n_detected": len(detected),
                "coverage": coverage,
                "coverage_tier": coverage_tier(coverage),
                "missing_genes": ";".join(missing),
            }
        )
    coverage = pd.DataFrame(coverage_rows).sort_values("program_id")
    qc = repo / "results" / "qc"
    gates.to_csv(qc / "gse256398_donor_gate_summary.csv", index=False)
    coverage.to_csv(qc / "gse256398_program_coverage.csv", index=False)
    summary = {
        "dataset_id": "GSE256398_human",
        "genes": len(genes),
        "pseudobulk_columns": len(manifest),
        "pseudobulk_nonzero": n_nonzero,
        "donors": int(manifest["donor_id"].nunique()),
        "contrasts": list(CONTRASTS),
        "all_contrast_lineage_30_cell_gates_pass": bool(gates["formal_30_cell_gate"].eq("PASS").all()),
        "programs_primary_coverage": int(coverage["coverage_tier"].eq("primary").sum()),
        "programs_total": len(coverage),
        "analysis_role": "post-lock reconstructed-label external validation",
    }
    (qc / "gse256398_gate_audit_summary.json").write_text(
        json.dumps(summary, indent=2) + "\n", encoding="utf-8"
    )
    print(gates.to_string(index=False))
    print(coverage[["program_id", "coverage", "coverage_tier"]].to_string(index=False))
    print(json.dumps(summary, indent=2))


if __name__ == "__main__":
    main()
