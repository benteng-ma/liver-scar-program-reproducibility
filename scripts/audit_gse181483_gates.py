from __future__ import annotations

import json
from pathlib import Path

import pandas as pd

from audit_gse244832_gates import coverage_tier, read_matrix_market_dimensions


CONTRAST = "cirrhosis_vs_healthy_directional"


def contrast_groups(manifest: pd.DataFrame) -> dict[str, pd.Series]:
    groups = manifest["disease_group"].map({"healthy": "control", "cirrhosis": "case"}).fillna("excluded")
    return {CONTRAST: groups}


def main() -> None:
    repo = Path(__file__).resolve().parents[1]
    interim = repo / "data" / "interim" / "GSE181483"
    manifest = pd.read_csv(interim / "donor_lineage_manifest.csv")
    genes = pd.read_csv(interim / "genes.csv")["gene"].astype(str)
    n_rows, n_columns, n_nonzero = read_matrix_market_dimensions(interim / "donor_lineage_raw_counts.mtx")
    if (n_rows, n_columns) != (len(genes), len(manifest)):
        raise RuntimeError("Pseudobulk matrix dimensions do not match genes/groups")
    if manifest["donor_id"].nunique() != 4 or len(manifest) != 12:
        raise RuntimeError("Expected four donors x three target lineages")
    manifest["eligible_30"] = manifest["n_cells"].ge(30)
    manifest["eligible_20"] = manifest["n_cells"].ge(20)
    manifest.to_csv(repo / "metadata" / "gse181483_donor_lineage_eligibility.csv", index=False)
    groups = contrast_groups(manifest)[CONTRAST]
    gate_rows: list[dict[str, object]] = []
    for lineage, lineage_data in manifest.assign(comparison_group=groups).groupby("harmonized_lineage"):
        counts: dict[tuple[str, int], int] = {}
        for group in ("control", "case"):
            values = lineage_data[lineage_data["comparison_group"] == group]
            counts[(group, 30)] = int(values["eligible_30"].sum())
            counts[(group, 20)] = int(values["eligible_20"].sum())
            gate_rows.append(
                {
                    "dataset_id": "GSE181483_human",
                    "contrast": CONTRAST,
                    "lineage": lineage,
                    "comparison_group": group,
                    "donors_before_cell_gate": len(values),
                    "donors_30": counts[(group, 30)],
                    "donors_20": counts[(group, 20)],
                    "minimum_cells": int(values["n_cells"].min()),
                    "median_cells": float(values["n_cells"].median()),
                    "maximum_cells": int(values["n_cells"].max()),
                }
            )
        directional30 = all(counts[(group, 30)] == 2 for group in ("control", "case"))
        directional20 = all(counts[(group, 20)] == 2 for group in ("control", "case"))
        for row in gate_rows[-2:]:
            row["formal_three_donor_gate"] = "FAIL"
            row["directional_30_cell_gate"] = "PASS" if directional30 else "FAIL"
            row["directional_20_cell_gate"] = "PASS" if directional20 else "FAIL"
            row["analysis_role"] = "two_plus_two_directional_only"
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
                "dataset_id": "GSE181483_human",
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
    qc_dir = repo / "results" / "qc"
    gates.to_csv(qc_dir / "gse181483_donor_gate_summary.csv", index=False)
    coverage.to_csv(qc_dir / "gse181483_program_coverage.csv", index=False)
    summary = {
        "dataset_id": "GSE181483_human",
        "genes_shared_all_four_donors": len(genes),
        "pseudobulk_columns": len(manifest),
        "pseudobulk_nonzero": n_nonzero,
        "donors": int(manifest["donor_id"].nunique()),
        "programs_primary_coverage": int(coverage["coverage_tier"].eq("primary").sum()),
        "programs_total": len(coverage),
        "formal_gate_status": "FAIL for every lineage because n=2 per group",
        "endpoint_status": "cirrhosis versus healthy directional display only",
    }
    (qc_dir / "gse181483_gate_audit_summary.json").write_text(json.dumps(summary, indent=2) + "\n", encoding="utf-8")
    print(gates.to_string(index=False))
    print(coverage[["program_id", "coverage", "coverage_tier"]].to_string(index=False))
    print(json.dumps(summary, indent=2))


if __name__ == "__main__":
    main()
