from __future__ import annotations

import json
from pathlib import Path

import pandas as pd


PRIMARY_CELLS = 30
SENSITIVITY_CELLS = 20
MIN_DONORS_PER_GROUP = 3


def canonical_donor_id(value: object) -> str:
    donor = str(value)
    return donor if donor.startswith("P") else f"P{donor}"


def coverage_tier(coverage: float) -> str:
    if coverage >= 0.80:
        return "primary"
    if coverage >= 0.60:
        return "evaluable_flagged"
    if coverage >= 0.40:
        return "sensitivity_only"
    return "not_evaluated"


def contrast_groups(manifest: pd.DataFrame) -> dict[str, pd.Series]:
    status = manifest["Disease.status"].astype(str)
    stage = manifest["Fibrosis.score..F0.4."].astype(int)
    clinical = status.map(
        lambda value: (
            "control"
            if value == "Healthy control"
            else "case" if value == "NASH with cirrhosis" else "excluded"
        )
    )
    advanced_non_end = pd.Series("excluded", index=manifest.index)
    advanced_non_end.loc[stage.eq(0)] = "control"
    advanced_non_end.loc[stage.isin([3, 4]) & status.ne("end stage")] = "case"
    end_stage = status.map(
        lambda value: (
            "control"
            if value == "Healthy control"
            else "case" if value == "end stage" else "excluded"
        )
    )
    all_fibrosis_non_end = pd.Series("excluded", index=manifest.index)
    all_fibrosis_non_end.loc[stage.eq(0)] = "control"
    all_fibrosis_non_end.loc[stage.isin([1, 2, 3, 4]) & status.ne("end stage")] = "case"
    advanced_pooled = pd.Series("excluded", index=manifest.index)
    advanced_pooled.loc[stage.eq(0)] = "control"
    advanced_pooled.loc[stage.isin([3, 4])] = "case"
    return {
        "clinical_cirrhosis_vs_healthy": clinical,
        "advanced_f3f4_vs_f0_non_end_stage": advanced_non_end,
        "end_stage_vs_healthy_separate_stratum": end_stage,
        "all_f1f4_vs_f0_non_end_stage_sensitivity": all_fibrosis_non_end,
        "advanced_f3f4_vs_f0_pooled_end_stage_sensitivity": advanced_pooled,
    }


def read_matrix_market_dimensions(path: Path) -> tuple[int, int, int]:
    with path.open("rt", encoding="ascii") as handle:
        banner = handle.readline().strip()
        if not banner.startswith("%%MatrixMarket matrix coordinate"):
            raise ValueError(f"unexpected Matrix Market banner: {banner}")
        for line in handle:
            if not line.startswith("%"):
                return tuple(map(int, line.split()))
    raise ValueError("Matrix Market dimensions are missing")


def dataframe_to_markdown(frame: pd.DataFrame) -> str:
    def clean(value: object) -> str:
        return str(value).replace("|", "\\|").replace("\n", " ")

    header = "| " + " | ".join(map(clean, frame.columns)) + " |"
    divider = "| " + " | ".join("---" for _ in frame.columns) + " |"
    rows = [
        "| " + " | ".join(clean(value) for value in row) + " |"
        for row in frame.itertuples(index=False, name=None)
    ]
    return "\n".join([header, divider, *rows])


def main() -> None:
    repo = Path(__file__).resolve().parents[1]
    interim = repo / "data" / "interim" / "GSE202379"
    manifest_path = interim / "donor_lineage_manifest.csv"
    matrix_path = interim / "donor_lineage_raw_counts.mtx"
    genes_path = interim / "genes.csv"
    inventory_path = repo / "literature" / "program_inventory.csv"
    object_inventory_path = interim / "object_inventory.tsv"

    manifest = pd.read_csv(manifest_path)
    genes = pd.read_csv(genes_path)["gene"].astype(str)
    n_rows, n_columns, n_nonzero = read_matrix_market_dimensions(matrix_path)
    if n_rows != len(genes) or n_columns != len(manifest):
        raise RuntimeError("pseudobulk matrix dimensions do not match genes/groups")
    if manifest["donor_id"].nunique() != 47:
        raise RuntimeError("expected 47 biological donors in the author object")
    per_donor_lineages = manifest.groupby("donor_id")["harmonized_lineage"].nunique()
    if not per_donor_lineages.eq(3).all():
        raise RuntimeError("each donor must have all three target lineage rows before gating")

    manifest.insert(
        2, "canonical_donor_id", manifest["donor_id"].map(canonical_donor_id)
    )
    manifest["eligible_primary_30"] = manifest["n_cells"].ge(PRIMARY_CELLS)
    manifest["eligible_sensitivity_20"] = manifest["n_cells"].ge(SENSITIVITY_CELLS)
    eligibility_columns = [
        "group_id",
        "donor_id",
        "canonical_donor_id",
        "harmonized_lineage",
        "n_cells",
        "eligible_primary_30",
        "eligible_sensitivity_20",
        "orig.ident",
        "manuscript.expt",
        "Disease.status",
        "Fibrosis.score..F0.4.",
        "Steatosis",
        "Ballooning",
        "Inflammation",
        "Age",
        "Gender",
    ]
    eligibility = manifest[[column for column in eligibility_columns if column in manifest]]
    eligibility_path = repo / "metadata" / "gse202379_donor_lineage_eligibility.csv"
    eligibility.to_csv(eligibility_path, index=False)

    gate_rows: list[dict[str, object]] = []
    for contrast, group in contrast_groups(manifest).items():
        scoped = manifest.assign(comparison_group=group)
        scoped = scoped[scoped["comparison_group"].ne("excluded")]
        for lineage, lineage_data in scoped.groupby("harmonized_lineage"):
            counts: dict[tuple[str, str], int] = {}
            for label in ("control", "case"):
                group_data = lineage_data[lineage_data["comparison_group"].eq(label)]
                counts[(label, "all")] = len(group_data)
                counts[(label, "primary")] = int(group_data["eligible_primary_30"].sum())
                counts[(label, "sensitivity")] = int(
                    group_data["eligible_sensitivity_20"].sum()
                )
                gate_rows.append(
                    {
                        "dataset_id": "GSE202379",
                        "contrast": contrast,
                        "lineage": lineage,
                        "comparison_group": label,
                        "donors_before_cell_gate": counts[(label, "all")],
                        "donors_primary_30": counts[(label, "primary")],
                        "donors_sensitivity_20": counts[(label, "sensitivity")],
                        "minimum_cells": int(group_data["n_cells"].min()),
                        "median_cells": float(group_data["n_cells"].median()),
                        "maximum_cells": int(group_data["n_cells"].max()),
                    }
                )
            primary_pass = all(
                counts[(label, "primary")] >= MIN_DONORS_PER_GROUP
                for label in ("control", "case")
            )
            sensitivity_pass = all(
                counts[(label, "sensitivity")] >= MIN_DONORS_PER_GROUP
                for label in ("control", "case")
            )
            for row in gate_rows[-2:]:
                row["formal_primary_gate"] = "PASS" if primary_pass else "FAIL"
                row["formal_sensitivity_gate"] = (
                    "PASS" if sensitivity_pass else "FAIL"
                )

    gates = pd.DataFrame(gate_rows)
    qc_dir = repo / "results" / "qc"
    qc_dir.mkdir(parents=True, exist_ok=True)
    gates_path = qc_dir / "gse202379_donor_gate_summary.csv"
    gates.to_csv(gates_path, index=False)

    programs = pd.read_csv(inventory_path)
    measured = set(genes.str.upper())
    coverage_rows = []
    for (program_id, lineage), rows in programs.groupby(["program_id", "cell_lineage"]):
        program_genes = set(rows["gene_symbol"].dropna().astype(str).str.upper())
        detected = sorted(program_genes & measured)
        missing = sorted(program_genes - measured)
        coverage = len(detected) / len(program_genes)
        coverage_rows.append(
            {
                "dataset_id": "GSE202379",
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
    coverage_path = qc_dir / "gse202379_program_coverage.csv"
    coverage.to_csv(coverage_path, index=False)

    donor_summary = (
        manifest.drop_duplicates("donor_id")[
            ["donor_id", "Disease.status", "Fibrosis.score..F0.4."]
        ]
        .groupby(["Disease.status", "Fibrosis.score..F0.4."], dropna=False)
        .size()
        .reset_index(name="donors")
    )
    source_path = repo / "data" / "raw" / "GSE202379" / "GSE202379_SeuratObject_AllCells.rds.gz"
    unwrapped_path = repo / "data" / "raw" / "GSE202379" / "GSE202379_SeuratObject_AllCells.rds"
    download_info = json.loads(
        source_path.with_suffix(source_path.suffix + ".download.json").read_text()
    )
    unwrap_info = json.loads(
        unwrapped_path.with_suffix(unwrapped_path.suffix + ".unwrap.json").read_text()
    )
    object_inventory = dict(
        line.split("\t", 1)
        for line in object_inventory_path.read_text(encoding="utf-8").splitlines()
        if "\t" in line
    )
    ingest_summary = {
        "dataset_id": "GSE202379",
        "author_object_cells": int(object_inventory["cells"]),
        "author_object_genes": len(genes),
        "author_object_donors": int(manifest["donor_id"].nunique()),
        "target_cells": int(manifest["n_cells"].sum()),
        "target_donor_lineage_groups": len(manifest),
        "pseudobulk_nonzero": n_nonzero,
        "source_bytes": source_path.stat().st_size,
        "source_sha256": download_info["sha256"],
        "unwrapped_bytes": unwrapped_path.stat().st_size,
        "unwrapped_sha256": unwrap_info["output_sha256"],
        "r_version": object_inventory["r_version"],
        "SeuratObject_version": object_inventory["SeuratObject_version"],
        "Matrix_version": object_inventory["Matrix_version"],
        "author_annotations_used": {
            "Macrophages": "macrophage_monocyte",
            "Endothelial": "endothelial",
            "Stellate": "mesenchymal_hsc_myofibroblast",
        },
        "primary_cell_gate": PRIMARY_CELLS,
        "sensitivity_cell_gate": SENSITIVITY_CELLS,
        "minimum_donors_per_group": MIN_DONORS_PER_GROUP,
        "programs_primary_coverage": int(coverage["coverage_tier"].eq("primary").sum()),
        "programs_total": len(coverage),
    }
    summary_path = qc_dir / "gse202379_ingest_summary.json"
    summary_path.write_text(
        json.dumps(ingest_summary, indent=2) + "\n", encoding="utf-8"
    )

    primary_gate_table = (
        gates.drop_duplicates(["contrast", "lineage"])[
            ["contrast", "lineage", "formal_primary_gate", "formal_sensitivity_gate"]
        ]
        .sort_values(["contrast", "lineage"])
    )
    report_lines = [
        "# GSE202379 annotation and donor-gate audit",
        "",
        "The author Seurat object was ingested without redefining cell labels. Author `Macrophages`, `Endothelial`, and `Stellate` annotations map exactly to the three frozen broad lineages. Technical regions/lobes were collapsed inside the author object by `Patient.ID` before inference.",
        "",
        f"- Author object: {int(object_inventory['cells']):,} nuclei, {len(genes):,} genes, 47 biological donors.",
        f"- Target lineages: {int(manifest['n_cells'].sum()):,} nuclei in 141 donor-lineage pseudobulks.",
        f"- Frozen cell gates: {PRIMARY_CELLS} nuclei primary; {SENSITIVITY_CELLS} nuclei sensitivity.",
        f"- Program coverage: {int(coverage['coverage_tier'].eq('primary').sum())}/{len(coverage)} programs meet the 80% primary threshold.",
        "- End-stage donors are not pooled into the clinical cirrhosis endpoint; they remain a named separate stratum. A pooled F3-F4 analysis is labelled sensitivity only.",
        "",
        "## Donor strata recovered from the author object",
        "",
        dataframe_to_markdown(donor_summary),
        "",
        "## Contrast-by-lineage gate status",
        "",
        dataframe_to_markdown(primary_gate_table),
        "",
        "`PASS` means both groups retain at least three donors after the applicable donor×lineage cell-count gate. This audit unlocks only the passing GSE202379 contrast-lineage combinations; it does not unlock other cohorts.",
        "",
    ]
    report_path = repo / "reports" / "phase2_gse202379_gate_audit.md"
    report_path.write_text("\n".join(report_lines), encoding="utf-8")

    print(primary_gate_table.to_string(index=False))
    print(coverage[["program_id", "coverage", "coverage_tier"]].to_string(index=False))
    print(json.dumps(ingest_summary, indent=2))


if __name__ == "__main__":
    main()
