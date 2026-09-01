from __future__ import annotations

import json
from pathlib import Path

import pandas as pd


PRIMARY_CELLS = 30
SENSITIVITY_CELLS = 20
MIN_DONORS_PER_GROUP = 3


def coverage_tier(coverage: float) -> str:
    if coverage >= 0.80:
        return "primary"
    if coverage >= 0.60:
        return "evaluable_flagged"
    if coverage >= 0.40:
        return "sensitivity_only"
    return "not_evaluated"


def read_matrix_market_dimensions(path: Path) -> tuple[int, int, int]:
    with path.open("rt", encoding="ascii") as handle:
        banner = handle.readline().strip()
        if not banner.startswith("%%MatrixMarket matrix coordinate"):
            raise ValueError(f"unexpected Matrix Market banner: {banner}")
        for line in handle:
            if not line.startswith("%"):
                return tuple(map(int, line.split()))
    raise ValueError("Matrix Market dimensions are missing")


def main() -> None:
    repo = Path(__file__).resolve().parents[1]
    interim = repo / "data" / "interim" / "GSE244832"
    manifest_path = interim / "donor_lineage_manifest.csv"
    matrix_path = interim / "donor_lineage_raw_counts.mtx"
    genes_path = interim / "genes.csv"
    manifest = pd.read_csv(manifest_path)
    genes = pd.read_csv(genes_path)["gene"].astype(str)
    n_rows, n_columns, n_nonzero = read_matrix_market_dimensions(matrix_path)
    if n_rows != len(genes) or n_columns != len(manifest):
        raise RuntimeError("pseudobulk matrix dimensions do not match genes/groups")
    if manifest["donor_id"].nunique() != 18 or len(manifest) != 54:
        raise RuntimeError("expected 18 donors x 3 target lineages")
    if not manifest.groupby("donor_id")["harmonized_lineage"].nunique().eq(3).all():
        raise RuntimeError("each donor must have all three target lineage rows")

    manifest["eligible_30"] = manifest["n_cells"].ge(PRIMARY_CELLS)
    manifest["eligible_20"] = manifest["n_cells"].ge(SENSITIVITY_CELLS)
    eligibility_path = repo / "metadata" / "gse244832_donor_lineage_eligibility.csv"
    manifest.to_csv(eligibility_path, index=False)

    scoped = manifest[manifest["disease_group"].isin(["normal", "MASH"])].copy()
    scoped["comparison_group"] = scoped["disease_group"].map(
        {"normal": "control", "MASH": "case"}
    )
    gate_rows: list[dict[str, object]] = []
    for lineage, lineage_data in scoped.groupby("harmonized_lineage"):
        group_counts: dict[tuple[str, int], int] = {}
        for group in ("control", "case"):
            values = lineage_data[lineage_data["comparison_group"] == group]
            group_counts[(group, 30)] = int(values["eligible_30"].sum())
            group_counts[(group, 20)] = int(values["eligible_20"].sum())
            gate_rows.append(
                {
                    "dataset_id": "GSE244832",
                    "contrast": "mash_f2f4_group_vs_normal_sensitivity",
                    "lineage": lineage,
                    "comparison_group": group,
                    "donors_before_cell_gate": len(values),
                    "donors_30": group_counts[(group, 30)],
                    "donors_20": group_counts[(group, 20)],
                    "minimum_cells": int(values["n_cells"].min()),
                    "median_cells": float(values["n_cells"].median()),
                    "maximum_cells": int(values["n_cells"].max()),
                }
            )
        pass_30 = all(group_counts[(group, 30)] >= MIN_DONORS_PER_GROUP for group in ("control", "case"))
        pass_20 = all(group_counts[(group, 20)] >= MIN_DONORS_PER_GROUP for group in ("control", "case"))
        for row in gate_rows[-2:]:
            row["formal_30_cell_gate"] = "PASS" if pass_30 else "FAIL"
            row["formal_20_cell_gate"] = "PASS" if pass_20 else "FAIL"
            row["analysis_role"] = "sensitivity_only_no_donor_stage"
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
                "dataset_id": "GSE244832",
                "program_id": program_id,
                "lineage": lineage,
                "n_program_genes": len(program_genes),
                "n_detected": len(detected),
                "coverage": coverage,
                "coverage_tier": coverage_tier(coverage),
                "missing_genes": ";".join(missing),
            }
        )
    coverage_table = pd.DataFrame(coverage_rows).sort_values("program_id")
    qc_dir = repo / "results" / "qc"
    qc_dir.mkdir(parents=True, exist_ok=True)
    gates.to_csv(qc_dir / "gse244832_donor_gate_summary.csv", index=False)
    coverage_table.to_csv(qc_dir / "gse244832_program_coverage.csv", index=False)
    summary = {
        "dataset_id": "GSE244832",
        "genes": len(genes),
        "pseudobulk_columns": len(manifest),
        "pseudobulk_nonzero": n_nonzero,
        "donors": int(manifest["donor_id"].nunique()),
        "condition_donors": manifest.drop_duplicates("donor_id")["disease_group"].value_counts().to_dict(),
        "primary_cell_gate": PRIMARY_CELLS,
        "sensitivity_cell_gate": SENSITIVITY_CELLS,
        "minimum_donors_per_group": MIN_DONORS_PER_GROUP,
        "programs_primary_coverage": int(coverage_table["coverage_tier"].eq("primary").sum()),
        "programs_total": len(coverage_table),
        "endpoint_status": "MASH F2-F4 group versus normal sensitivity only; donor-level fibrosis stage unavailable",
    }
    (qc_dir / "gse244832_gate_audit_summary.json").write_text(
        json.dumps(summary, indent=2) + "\n", encoding="utf-8"
    )
    print(gates.to_string(index=False))
    print(coverage_table[["program_id", "coverage", "coverage_tier"]].to_string(index=False))
    print(json.dumps(summary, indent=2))


if __name__ == "__main__":
    main()
