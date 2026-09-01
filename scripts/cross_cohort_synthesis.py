from __future__ import annotations

import hashlib
import json
from pathlib import Path

import numpy as np
import pandas as pd
from scipy.optimize import minimize_scalar


COHORT_META = {
    "GSE202379": {
        "assay": "snRNA-seq",
        "annotation_role": "author labels",
    },
    "GSE290642_human": {
        "assay": "scRNA-seq",
        "annotation_role": "reconstructed broad labels; sensitivity only",
    },
    "GSE244832": {
        "assay": "snRNA-seq",
        "annotation_role": "author cluster mapping; mixed-stage MASH sensitivity",
    },
    "GSE210077_Watson6": {
        "assay": "snRNA-seq",
        "annotation_role": "author labels; mixed F2/F3/F4 sensitivity",
    },
    "GSE181483_human": {
        "assay": "scRNA-seq",
        "annotation_role": "reconstructed broad labels; 2+2 directional only",
    },
}


def endpoint_family(dataset: str, contrast: str) -> str:
    if dataset == "GSE202379" and contrast == "clinical_cirrhosis_vs_healthy":
        return "clinical_cirrhosis_vs_healthy"
    if dataset == "GSE181483_human":
        return "clinical_cirrhosis_vs_healthy_directional"
    if dataset == "GSE202379" and contrast == "advanced_f3f4_vs_f0_non_end_stage":
        return "advanced_f3f4_vs_f0"
    if dataset == "GSE290642_human" and contrast == "f4_vs_f0_reconstructed_label_sensitivity":
        return "advanced_f3f4_vs_f0_reconstructed_sensitivity"
    if dataset == "GSE290642_human":
        return "all_fibrosis_vs_f0_sensitivity"
    if dataset == "GSE244832":
        return "mash_f2f4_group_vs_normal_sensitivity"
    if dataset == "GSE210077_Watson6":
        return "mixed_f2f4_vs_healthy_sensitivity"
    return "unmapped"


def reml_tau2(y: np.ndarray, variance: np.ndarray) -> float:
    upper = max(10.0, float(np.var(y, ddof=1) * 100), float(np.max(variance) * 100))

    def objective(tau2: float) -> float:
        weights = 1.0 / (variance + tau2)
        mean = float(np.sum(weights * y) / np.sum(weights))
        return 0.5 * (
            float(np.log(variance + tau2).sum())
            + float(np.log(weights.sum()))
            + float(np.sum(weights * (y - mean) ** 2))
        )

    result = minimize_scalar(objective, bounds=(0.0, upper), method="bounded")
    at_zero = objective(0.0)
    if not result.success or at_zero <= result.fun:
        return 0.0
    return float(max(result.x, 0.0))


def meta_row(program_id: str, score_method: str, inputs: pd.DataFrame) -> dict[str, object]:
    y = inputs["hedges_g"].to_numpy(dtype=float)
    se = inputs["robust_se_g_hc3"].to_numpy(dtype=float)
    variance = se**2
    fixed_weights = 1.0 / variance
    fixed_mean = float(np.sum(fixed_weights * y) / fixed_weights.sum())
    fixed_se = float(np.sqrt(1.0 / fixed_weights.sum()))
    q = float(np.sum(fixed_weights * (y - fixed_mean) ** 2))
    df = len(y) - 1
    i2 = float(max(0.0, (q - df) / q) * 100) if q > 0 else 0.0
    tau2 = reml_tau2(y, variance)
    random_weights = 1.0 / (variance + tau2)
    random_mean = float(np.sum(random_weights * y) / random_weights.sum())
    random_se = float(np.sqrt(1.0 / random_weights.sum()))
    return {
        "endpoint_family": "advanced_f3f4_vs_f0_sensitivity_meta",
        "lineage": "endothelial",
        "program_id": program_id,
        "score_method": score_method,
        "k": len(y),
        "fixed_hedges_g": fixed_mean,
        "fixed_se": fixed_se,
        "fixed_ci95_low": fixed_mean - 1.96 * fixed_se,
        "fixed_ci95_high": fixed_mean + 1.96 * fixed_se,
        "random_reml_hedges_g": random_mean,
        "random_reml_se": random_se,
        "random_reml_ci95_low": random_mean - 1.96 * random_se,
        "random_reml_ci95_high": random_mean + 1.96 * random_se,
        "q": q,
        "i2_percent": i2,
        "tau2_reml": tau2,
        "maximum_fixed_weight_fraction": float(fixed_weights.max() / fixed_weights.sum()),
        "maximum_random_weight_fraction": float(random_weights.max() / random_weights.sum()),
        "hartung_knapp_status": "not applied: only two cohorts",
        "leave_one_study_out_status": "not applied: only two cohorts",
        "formal_replication_eligible": False,
        "eligibility_reason": "one input uses reconstructed labels and is sensitivity only; k=2",
        "input_datasets": ";".join(inputs["dataset_id"].astype(str)),
        "input_effects": ";".join(f"{value:.6g}" for value in y),
        "input_robust_se": ";".join(f"{value:.6g}" for value in se),
    }


def main() -> None:
    repo = Path(__file__).resolve().parents[1]
    files = [
        repo / "results" / "primary" / "gse202379_primary_effects.csv",
        repo / "results" / "sensitivity" / "gse202379_sensitivity_effects.csv",
        repo / "results" / "sensitivity" / "gse290642_sensitivity_effects.csv",
        repo / "results" / "sensitivity" / "gse244832_sensitivity_effects.csv",
        repo / "results" / "sensitivity" / "gse210077_watson6_sensitivity_effects.csv",
        repo / "results" / "sensitivity" / "gse181483_directional_effects.csv",
    ]
    frames: list[pd.DataFrame] = []
    for path in files:
        frame = pd.read_csv(path)
        frame["source_file"] = path.relative_to(repo).as_posix()
        frames.append(frame)
    effects = pd.concat(frames, ignore_index=True, sort=False)
    effects["assay"] = effects["dataset_id"].map(lambda value: COHORT_META[value]["assay"])
    effects["annotation_role"] = effects["dataset_id"].map(lambda value: COHORT_META[value]["annotation_role"])
    effects["endpoint_family"] = [endpoint_family(d, c) for d, c in zip(effects["dataset_id"], effects["contrast"])]
    effects["expected_direction_observed"] = effects["hedges_g"].gt(0)
    effects["formal_primary_row"] = effects["analysis_tier"].eq("primary")
    effects["comparable_advanced_sensitivity_input"] = (
        ((effects["dataset_id"] == "GSE202379") & (effects["contrast"] == "advanced_f3f4_vs_f0_non_end_stage"))
        | ((effects["dataset_id"] == "GSE290642_human") & (effects["contrast"] == "f4_vs_f0_reconstructed_label_sensitivity"))
    ) & effects["lineage"].eq("endothelial")
    meta_dir = repo / "results" / "meta"
    meta_dir.mkdir(parents=True, exist_ok=True)
    matrix_path = meta_dir / "cross_cohort_effect_matrix.csv"
    effects.to_csv(matrix_path, index=False)

    meta_inputs = effects[effects["comparable_advanced_sensitivity_input"]].copy()
    meta_inputs = meta_inputs[meta_inputs["program_id"] != "RAM2019_ENDO_7_SAENDO2"]
    meta_rows: list[dict[str, object]] = []
    for (program_id, score_method), values in meta_inputs.groupby(["program_id", "score_method"], sort=True):
        if values["dataset_id"].nunique() != 2:
            raise RuntimeError(f"Expected two advanced-fibrosis cohorts for {program_id}/{score_method}")
        meta_rows.append(meta_row(program_id, score_method, values.sort_values("dataset_id")))
    meta = pd.DataFrame(meta_rows).sort_values(["program_id", "score_method"])
    meta_path = meta_dir / "advanced_endothelial_sensitivity_meta.csv"
    meta.to_csv(meta_path, index=False)

    programs = pd.read_csv(repo / "literature" / "program_inventory.csv")
    classification_rows: list[dict[str, object]] = []
    for program_id, program_rows in programs.groupby("program_id", sort=True):
        lineage = str(program_rows["cell_lineage"].iloc[0])
        subset = effects[effects["program_id"] == program_id]
        cohort_method = (
            subset.groupby(["dataset_id", "score_method"])["expected_direction_observed"]
            .any()
            .unstack("score_method", fill_value=False)
        )
        for method in ("singscore", "standardized_mean"):
            if method not in cohort_method:
                cohort_method[method] = False
        both = cohort_method["singscore"] & cohort_method["standardized_mean"]
        support_datasets = sorted(cohort_method.index[both].astype(str))
        formal_primary_datasets = sorted(
            subset.loc[subset["formal_primary_row"] & subset["expected_direction_observed"], "dataset_id"].unique()
        )
        classification = "INDEPENDENT_DIRECTIONAL_SUPPORT" if support_datasets else "NO_EXPECTED_DIRECTION_SUPPORT"
        meta_subset = meta[meta["program_id"] == program_id]
        classification_rows.append(
            {
                "program_id": program_id,
                "lineage": lineage,
                "independent_datasets_both_scores_positive": len(support_datasets),
                "support_datasets": ";".join(support_datasets),
                "formal_primary_expected_direction_datasets": ";".join(formal_primary_datasets),
                "classification": classification,
                "within_cell_state_replicated": False,
                "within_cell_state_reason": "no program has two comparable author-label formal-primary cohorts with an eligible meta-analysis",
                "advanced_sensitivity_meta_available": not meta_subset.empty,
                "advanced_sensitivity_fixed_ci_excludes_zero_both_scores": bool(
                    len(meta_subset) == 2 and (meta_subset["fixed_ci95_low"] > 0).all()
                ),
                "advanced_sensitivity_random_ci_excludes_zero_both_scores": bool(
                    len(meta_subset) == 2 and (meta_subset["random_reml_ci95_low"] > 0).all()
                ),
                "advanced_sensitivity_weight_rule_passes_both_scores": bool(
                    len(meta_subset) == 2 and (meta_subset["maximum_fixed_weight_fraction"] <= 0.70).all()
                ),
                "pan_cirrhotic_transportable": False,
                "pan_cirrhotic_reason": "within-state replication prerequisite not met; etiology and assay requirements unresolved",
                "assay_transfer_class": "UNRESOLVED",
                "assay_transfer_reason": "no comparable eligible scRNA/snRNA endpoint pair supports a formal assay-transfer label",
                "spatial_status": "NOT_EVALUATED_INELIGIBLE_PANEL",
            }
        )
    classifications = pd.DataFrame(classification_rows)
    classification_path = meta_dir / "program_classification_table.csv"
    classifications.to_csv(classification_path, index=False)

    summary = {
        "effect_rows": len(effects),
        "datasets": sorted(effects["dataset_id"].unique()),
        "advanced_endothelial_sensitivity_meta_rows": len(meta),
        "programs": len(classifications),
        "independent_directional_support_programs": int(classifications["classification"].eq("INDEPENDENT_DIRECTIONAL_SUPPORT").sum()),
        "within_cell_state_replicated_programs": int(classifications["within_cell_state_replicated"].sum()),
        "pan_cirrhotic_transportable_programs": int(classifications["pan_cirrhotic_transportable"].sum()),
        "assay_robust_programs": int(classifications["assay_transfer_class"].eq("ASSAY_ROBUST").sum()),
        "primary_meta_status": "not run: no endpoint has two comparable author-label formal-primary cohorts",
        "sensitivity_meta_status": "advanced endothelial only: GSE202379 F3-F4 vs F0 plus reconstructed-label GSE290642 F4 vs F0",
        "matrix_sha256": hashlib.sha256(matrix_path.read_bytes()).hexdigest().upper(),
        "meta_sha256": hashlib.sha256(meta_path.read_bytes()).hexdigest().upper(),
        "classification_sha256": hashlib.sha256(classification_path.read_bytes()).hexdigest().upper(),
    }
    (repo / "results" / "logs" / "cross_cohort_synthesis_run.json").write_text(json.dumps(summary, indent=2) + "\n", encoding="utf-8")
    print(meta.to_string(index=False))
    print(classifications[["program_id", "classification", "support_datasets", "within_cell_state_replicated", "assay_transfer_class"]].to_string(index=False))
    print(json.dumps(summary, indent=2))


if __name__ == "__main__":
    main()
