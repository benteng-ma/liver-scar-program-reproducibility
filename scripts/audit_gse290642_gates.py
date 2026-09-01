from __future__ import annotations

import json
from pathlib import Path

import pandas as pd

from audit_gse244832_gates import coverage_tier, read_matrix_market_dimensions


PRIMARY_CELLS = 30
SENSITIVITY_CELLS = 20
MIN_DONORS_PER_GROUP = 3


def contrast_groups(manifest: pd.DataFrame) -> dict[str, pd.Series]:
    stage = manifest["fibrosis_stage"].astype(str)
    f4 = pd.Series("excluded", index=manifest.index)
    f4.loc[stage.eq("F0")] = "control"
    f4.loc[stage.eq("F4")] = "case"
    all_fibrosis = pd.Series("excluded", index=manifest.index)
    all_fibrosis.loc[stage.eq("F0")] = "control"
    all_fibrosis.loc[stage.ne("F0")] = "case"
    return {
        "f4_vs_f0_reconstructed_label_sensitivity": f4,
        "all_fibrosis_vs_f0_reconstructed_label_sensitivity": all_fibrosis,
    }


def main() -> None:
    repo = Path(__file__).resolve().parents[1]
    interim = repo / "data" / "interim" / "GSE290642"
    manifest = pd.read_csv(interim / "donor_lineage_manifest.csv")
    genes = pd.read_csv(interim / "genes.csv")["gene"].astype(str)
    shared_genes = pd.read_csv(interim / "shared_genes.csv")["gene"].astype(str)
    n_rows, n_columns, n_nonzero = read_matrix_market_dimensions(interim / "donor_lineage_raw_counts.mtx")
    if (n_rows, n_columns) != (len(genes), len(manifest)):
        raise RuntimeError("pseudobulk matrix dimensions do not match genes/groups")
    if manifest["donor_id"].nunique() != 24 or len(manifest) != 72:
        raise RuntimeError("expected 24 donors x 3 target lineages")
    manifest["eligible_30"] = manifest["n_cells"].ge(PRIMARY_CELLS)
    manifest["eligible_20"] = manifest["n_cells"].ge(SENSITIVITY_CELLS)
    manifest.to_csv(repo / "metadata" / "gse290642_donor_lineage_eligibility.csv", index=False)

    gate_rows: list[dict[str, object]] = []
    for contrast, groups in contrast_groups(manifest).items():
        scoped = manifest.assign(comparison_group=groups)
        scoped = scoped[scoped["comparison_group"] != "excluded"]
        for lineage, lineage_data in scoped.groupby("harmonized_lineage"):
            counts: dict[tuple[str, int], int] = {}
            for group in ("control", "case"):
                values = lineage_data[lineage_data["comparison_group"] == group]
                counts[(group, 30)] = int(values["eligible_30"].sum())
                counts[(group, 20)] = int(values["eligible_20"].sum())
                gate_rows.append(
                    {
                        "dataset_id": "GSE290642_human",
                        "contrast": contrast,
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
            pass30 = all(counts[(group, 30)] >= MIN_DONORS_PER_GROUP for group in ("control", "case"))
            pass20 = all(counts[(group, 20)] >= MIN_DONORS_PER_GROUP for group in ("control", "case"))
            for row in gate_rows[-2:]:
                row["formal_30_cell_gate"] = "PASS" if pass30 else "FAIL"
                row["formal_20_cell_gate"] = "PASS" if pass20 else "FAIL"
                row["analysis_role"] = "reconstructed_label_sensitivity_only"
    gates = pd.DataFrame(gate_rows)

    programs = pd.read_csv(repo / "literature" / "program_inventory.csv")
    measured = set(shared_genes.str.upper())
    coverage_rows: list[dict[str, object]] = []
    for (program_id, lineage), rows in programs.groupby(["program_id", "cell_lineage"]):
        program_genes = set(rows["gene_symbol"].dropna().astype(str).str.upper())
        detected = sorted(program_genes & measured)
        missing = sorted(program_genes - measured)
        coverage = len(detected) / len(program_genes)
        coverage_rows.append(
            {
                "dataset_id": "GSE290642_human",
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
    gates.to_csv(qc_dir / "gse290642_donor_gate_summary.csv", index=False)
    coverage.to_csv(qc_dir / "gse290642_program_coverage.csv", index=False)
    summary = {
        "dataset_id": "GSE290642_human",
        "genes": len(genes),
        "genes_shared_all_24_donors": len(shared_genes),
        "pseudobulk_columns": len(manifest),
        "pseudobulk_nonzero": n_nonzero,
        "donors": int(manifest["donor_id"].nunique()),
        "stage_donors": manifest.drop_duplicates("donor_id")["fibrosis_stage"].value_counts().to_dict(),
        "programs_primary_coverage": int(coverage["coverage_tier"].eq("primary").sum()),
        "programs_total": len(coverage),
        "annotation_role": "reconstructed broad-lineage sensitivity only",
        "etiology": "unresolved",
    }
    (qc_dir / "gse290642_gate_audit_summary.json").write_text(
        json.dumps(summary, indent=2) + "\n", encoding="utf-8"
    )
    print(gates.to_string(index=False))
    print(coverage[["program_id", "coverage", "coverage_tier"]].to_string(index=False))
    print(json.dumps(summary, indent=2))


if __name__ == "__main__":
    main()
