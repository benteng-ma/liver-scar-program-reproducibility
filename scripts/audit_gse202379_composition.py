from __future__ import annotations

import json
from pathlib import Path

import pandas as pd
from scipy.stats import spearmanr

from analyze_gse202379_programs import effect_statistics, stable_seed
from audit_gse202379_gates import contrast_groups


def canonical_donor_id(value: object) -> str:
    donor = str(value)
    return donor if donor.startswith("P") else f"P{donor}"


def main() -> None:
    repo = Path(__file__).resolve().parents[1]
    interim = repo / "data" / "interim" / "GSE202379"
    metadata = pd.read_csv(interim / "cell_metadata.csv.gz", low_memory=False)
    metadata["canonical_donor_id"] = metadata["Patient.ID"].map(canonical_donor_id)
    lineage_map = {
        "Macrophages": "macrophage_monocyte",
        "Endothelial": "endothelial",
        "Stellate": "mesenchymal_hsc_myofibroblast",
    }
    metadata["harmonized_lineage"] = metadata["cell.annotation"].map(lineage_map)
    donor_totals = metadata.groupby("canonical_donor_id").size().rename("total_nuclei")
    target_counts = (
        metadata.dropna(subset=["harmonized_lineage"])
        .groupby(["canonical_donor_id", "harmonized_lineage"])
        .size()
        .rename("lineage_nuclei")
        .reset_index()
    )
    proportions = target_counts.merge(
        donor_totals.reset_index(), on="canonical_donor_id", validate="many_to_one"
    )
    proportions["lineage_proportion"] = (
        proportions["lineage_nuclei"] / proportions["total_nuclei"]
    )
    donor_metadata = (
        metadata[
            [
                "canonical_donor_id",
                "Disease.status",
                "Fibrosis.score..F0.4.",
                "Age",
                "Gender",
            ]
        ]
        .drop_duplicates()
    )
    if donor_metadata["canonical_donor_id"].duplicated().any():
        raise RuntimeError("donor clinical metadata are not constant")
    proportions = proportions.merge(
        donor_metadata, on="canonical_donor_id", validate="many_to_one"
    )
    output_dir = repo / "results" / "composition"
    output_dir.mkdir(parents=True, exist_ok=True)
    proportions.to_csv(
        output_dir / "gse202379_donor_lineage_proportions.csv", index=False
    )

    scores = pd.read_csv(
        repo / "results" / "primary" / "gse202379_donor_program_scores.csv.gz"
    )
    scores = scores.merge(
        proportions[
            ["canonical_donor_id", "harmonized_lineage", "lineage_proportion"]
        ],
        on=["canonical_donor_id", "harmonized_lineage"],
        validate="many_to_one",
    )
    correlation_rows = []
    for keys, rows in scores[scores["eligible_sensitivity_20"]].groupby(
        ["program_id", "lineage", "score_method"]
    ):
        rho, p_value = spearmanr(rows["score"], rows["lineage_proportion"])
        correlation_rows.append(
            {
                "dataset_id": "GSE202379",
                "program_id": keys[0],
                "lineage": keys[1],
                "score_method": keys[2],
                "eligible_donors": len(rows),
                "spearman_rho_score_vs_lineage_proportion": float(rho),
                "spearman_p_two_sided": float(p_value),
                "strong_descriptive_tracking": bool(abs(rho) >= 0.5 and p_value < 0.05),
            }
        )
    correlations = pd.DataFrame(correlation_rows).sort_values(
        ["lineage", "program_id", "score_method"]
    )
    correlations.to_csv(
        output_dir / "gse202379_score_composition_correlations.csv", index=False
    )

    manifest = pd.read_csv(interim / "donor_lineage_manifest.csv")
    manifest["canonical_donor_id"] = manifest["donor_id"].map(canonical_donor_id)
    manifest["eligible_primary_30"] = manifest["n_cells"].ge(30)
    manifest["eligible_sensitivity_20"] = manifest["n_cells"].ge(20)
    manifest = manifest.merge(
        proportions[
            ["canonical_donor_id", "harmonized_lineage", "lineage_proportion"]
        ],
        on=["canonical_donor_id", "harmonized_lineage"],
        validate="one_to_one",
    )
    all_effects = pd.concat(
        [
            pd.read_csv(repo / "results" / "primary" / "gse202379_primary_effects.csv"),
            pd.read_csv(
                repo
                / "results"
                / "sensitivity"
                / "gse202379_sensitivity_effects.csv"
            ),
        ],
        ignore_index=True,
    )
    effect_designs = all_effects[
        ["contrast", "lineage", "analysis_tier", "cell_gate"]
    ].drop_duplicates()
    groups = contrast_groups(manifest)
    proportion_effect_rows = []
    for design in effect_designs.itertuples(index=False):
        eligibility = (
            "eligible_primary_30" if design.cell_gate == 30 else "eligible_sensitivity_20"
        )
        group = groups[design.contrast]
        selected = (
            manifest["harmonized_lineage"].eq(design.lineage)
            & manifest[eligibility]
            & group.ne("excluded")
        )
        selected_rows = manifest[selected]
        selected_group = group[selected]
        control = selected_rows.loc[
            selected_group.eq("control"), "lineage_proportion"
        ].to_numpy()
        case = selected_rows.loc[
            selected_group.eq("case"), "lineage_proportion"
        ].to_numpy()
        stats = effect_statistics(
            control,
            case,
            stable_seed("GSE202379", design.contrast, design.lineage, "composition"),
        )
        proportion_effect_rows.append(
            {
                "dataset_id": "GSE202379",
                "contrast": design.contrast,
                "lineage": design.lineage,
                "analysis_tier": design.analysis_tier,
                "cell_gate": design.cell_gate,
                **stats,
            }
        )
    proportion_effects = pd.DataFrame(proportion_effect_rows).sort_values(
        ["analysis_tier", "contrast", "lineage"]
    )
    proportion_effects.to_csv(
        output_dir / "gse202379_lineage_proportion_effects.csv", index=False
    )
    summary = {
        "dataset_id": "GSE202379",
        "donors": proportions["canonical_donor_id"].nunique(),
        "lineages": proportions["harmonized_lineage"].nunique(),
        "strong_score_proportion_tracking_rows": int(
            correlations["strong_descriptive_tracking"].sum()
        ),
        "correlation_rows": len(correlations),
        "interpretation_boundary": (
            "Correlation with lineage abundance is descriptive; all program effects were "
            "estimated inside author-defined lineages and cannot be reduced to cell composition alone."
        ),
    }
    (repo / "results" / "logs" / "gse202379_composition_audit.json").write_text(
        json.dumps(summary, indent=2) + "\n", encoding="utf-8"
    )
    print(proportion_effects.to_string(index=False))
    print(correlations[correlations["strong_descriptive_tracking"]].to_string(index=False))
    print(json.dumps(summary, indent=2))


if __name__ == "__main__":
    main()
