from __future__ import annotations

import hashlib
import json
from pathlib import Path

import numpy as np
import pandas as pd
from scipy.stats import spearmanr

from cross_cohort_synthesis import reml_tau2


PAIRS = [
    {
        "pair_id": "mash_cirrhosis_endpoint_aligned_cross_cohort",
        "pair_class": "stage_and_etiology_aligned; author versus reconstructed labels",
        "left": ("GSE202379", "clinical_cirrhosis_vs_healthy"),
        "right": ("GSE256398_human", "mash_cirrhosis_vs_healthy"),
    },
    {
        "pair_id": "mash_cross_cohort_stage_mismatched",
        "pair_class": "etiology aligned; mixed F2-F4 versus cirrhosis endpoint",
        "left": ("GSE244832", "mash_f2f4_group_vs_normal_sensitivity"),
        "right": ("GSE256398_human", "mash_cirrhosis_vs_healthy"),
    },
    {
        "pair_id": "same_assay_cirrhosis_etiology_contrast",
        "pair_class": "same cohort, assay and endpoint; MASH versus alcohol etiology",
        "left": ("GSE256398_human", "alcohol_cirrhosis_vs_healthy"),
        "right": ("GSE256398_human", "mash_cirrhosis_vs_healthy"),
    },
    {
        "pair_id": "same_assay_mash_stage_contrast",
        "pair_class": "same cohort, assay and etiology; fibrosis versus cirrhosis case groups",
        "left": ("GSE256398_human", "mash_fibrosis_vs_masld_f0"),
        "right": ("GSE256398_human", "mash_cirrhosis_vs_masld_f0"),
    },
    {
        "pair_id": "advanced_fibrosis_cross_assay_reference",
        "pair_class": "existing endpoint-aligned scRNA versus snRNA sensitivity reference",
        "left": ("GSE202379", "advanced_f3f4_vs_f0_non_end_stage"),
        "right": ("GSE290642_human", "f4_vs_f0_reconstructed_label_sensitivity"),
    },
]


def combined(repo: Path) -> pd.DataFrame:
    phase2 = pd.read_csv(repo / "results" / "meta" / "cross_cohort_effect_matrix.csv")
    phase3 = pd.read_csv(repo / "results" / "phase3" / "gse256398_program_effects.csv")
    phase3["formal_primary_row"] = False
    phase3["source_file"] = "results/phase3/gse256398_program_effects.csv"
    phase3["assay"] = "snRNA-seq"
    phase3["annotation_role"] = "post-lock reconstructed broad labels"
    return pd.concat([phase2, phase3], ignore_index=True, sort=False)


def pairwise(effects: pd.DataFrame, output: Path) -> pd.DataFrame:
    rows: list[dict[str, object]] = []
    for specification in PAIRS:
        left_dataset, left_contrast = specification["left"]
        right_dataset, right_contrast = specification["right"]
        subset = effects[
            (effects["dataset_id"].eq(left_dataset) & effects["contrast"].eq(left_contrast))
            | (effects["dataset_id"].eq(right_dataset) & effects["contrast"].eq(right_contrast))
        ]
        for (lineage, method), values in subset.groupby(["lineage", "score_method"], sort=True):
            left = values[values["dataset_id"].eq(left_dataset) & values["contrast"].eq(left_contrast)].set_index("program_id")["hedges_g"]
            right = values[values["dataset_id"].eq(right_dataset) & values["contrast"].eq(right_contrast)].set_index("program_id")["hedges_g"]
            paired = pd.concat([left.rename("left"), right.rename("right")], axis=1).dropna()
            if len(paired) < 3:
                continue
            rows.append(
                {
                    "pair_id": specification["pair_id"],
                    "pair_class": specification["pair_class"],
                    "lineage": lineage,
                    "score_method": method,
                    "dataset_left": left_dataset,
                    "contrast_left": left_contrast,
                    "dataset_right": right_dataset,
                    "contrast_right": right_contrast,
                    "n_programs": len(paired),
                    "spearman_rho": float(spearmanr(paired["left"], paired["right"]).statistic),
                    "sign_agreement": float((np.sign(paired["left"]) == np.sign(paired["right"])).mean()),
                    "median_absolute_g_difference": float(np.median(np.abs(paired["left"] - paired["right"]))),
                }
            )
    result = pd.DataFrame(rows)
    result.to_csv(output / "phase3_selected_pair_concordance.csv", index=False)
    return result


def meta_row(program_id: str, lineage: str, method: str, values: pd.DataFrame) -> dict[str, object]:
    y = values["hedges_g"].to_numpy(dtype=float)
    variance = values["robust_se_g_hc3"].to_numpy(dtype=float) ** 2
    fixed_weights = 1 / variance
    fixed = float(np.sum(fixed_weights * y) / fixed_weights.sum())
    fixed_se = float(np.sqrt(1 / fixed_weights.sum()))
    q = float(np.sum(fixed_weights * (y - fixed) ** 2))
    i2 = float(max(0, (q - 1) / q) * 100) if q > 0 else 0.0
    tau2 = reml_tau2(y, variance)
    random_weights = 1 / (variance + tau2)
    random = float(np.sum(random_weights * y) / random_weights.sum())
    random_se = float(np.sqrt(1 / random_weights.sum()))
    return {
        "endpoint_family": "mash_cirrhosis_vs_healthy_post_lock_sensitivity_meta",
        "lineage": lineage,
        "program_id": program_id,
        "score_method": method,
        "k": 2,
        "fixed_hedges_g": fixed,
        "fixed_ci95_low": fixed - 1.96 * fixed_se,
        "fixed_ci95_high": fixed + 1.96 * fixed_se,
        "random_reml_hedges_g": random,
        "random_reml_ci95_low": random - 1.96 * random_se,
        "random_reml_ci95_high": random + 1.96 * random_se,
        "q": q,
        "i2_percent": i2,
        "tau2_reml": tau2,
        "maximum_fixed_weight_fraction": float(fixed_weights.max() / fixed_weights.sum()),
        "maximum_random_weight_fraction": float(random_weights.max() / random_weights.sum()),
        "input_datasets": ";".join(values["dataset_id"].astype(str)),
        "input_effects": ";".join(f"{value:.8g}" for value in y),
        "formal_replication_eligible": False,
        "eligibility_reason": "post-lock k=2 synthesis includes one reconstructed-label cohort",
    }


def cirrhosis_meta(effects: pd.DataFrame, output: Path) -> tuple[pd.DataFrame, pd.DataFrame]:
    selected = effects[
        (effects["dataset_id"].eq("GSE202379") & effects["contrast"].eq("clinical_cirrhosis_vs_healthy"))
        | (effects["dataset_id"].eq("GSE256398_human") & effects["contrast"].eq("mash_cirrhosis_vs_healthy"))
    ].copy()
    rows: list[dict[str, object]] = []
    for (lineage, program_id, method), values in selected.groupby(["lineage", "program_id", "score_method"], sort=True):
        if values["dataset_id"].nunique() != 2:
            continue
        rows.append(meta_row(program_id, lineage, method, values.sort_values("dataset_id")))
    meta = pd.DataFrame(rows)
    meta.to_csv(output / "phase3_mash_cirrhosis_sensitivity_meta.csv", index=False)
    pivot = meta.pivot_table(
        index=["lineage", "program_id"],
        columns="score_method",
        values=["fixed_ci95_low", "random_reml_ci95_low", "maximum_fixed_weight_fraction"],
        aggfunc="first",
    )
    summary_rows: list[dict[str, object]] = []
    for (lineage, program_id), row in pivot.iterrows():
        both_methods = all((metric, method) in row.index for metric in ("fixed_ci95_low", "random_reml_ci95_low", "maximum_fixed_weight_fraction") for method in ("singscore", "standardized_mean"))
        summary_rows.append(
            {
                "lineage": lineage,
                "program_id": program_id,
                "both_score_methods_available": both_methods,
                "fixed_ci_positive_both_scores": bool(both_methods and row[("fixed_ci95_low", "singscore")] > 0 and row[("fixed_ci95_low", "standardized_mean")] > 0),
                "random_reml_ci_positive_both_scores": bool(both_methods and row[("random_reml_ci95_low", "singscore")] > 0 and row[("random_reml_ci95_low", "standardized_mean")] > 0),
                "fixed_weight_max_le_070_both_scores": bool(both_methods and row[("maximum_fixed_weight_fraction", "singscore")] <= 0.70 and row[("maximum_fixed_weight_fraction", "standardized_mean")] <= 0.70),
                "formal_replication_eligible": False,
                "classification": "POST_LOCK_SENSITIVITY_SUPPORT" if both_methods and row[("fixed_ci95_low", "singscore")] > 0 and row[("fixed_ci95_low", "standardized_mean")] > 0 else "NO_DUAL_SCORE_POST_LOCK_META_SUPPORT",
            }
        )
    summary = pd.DataFrame(summary_rows)
    summary.to_csv(output / "phase3_mash_cirrhosis_meta_program_summary.csv", index=False)
    return meta, summary


def main() -> None:
    repo = Path(__file__).resolve().parents[1]
    output = repo / "results" / "phase3"
    output.mkdir(parents=True, exist_ok=True)
    effects = combined(repo)
    matrix_path = output / "phase3_cross_cohort_effect_matrix.csv"
    effects.to_csv(matrix_path, index=False)
    concordance = pairwise(effects, output)
    meta, program_summary = cirrhosis_meta(effects, output)
    summary = {
        "effect_rows": len(effects),
        "datasets": int(effects["dataset_id"].nunique()),
        "selected_pair_rows": len(concordance),
        "mash_cirrhosis_meta_rows": len(meta),
        "post_lock_fixed_dual_score_support_programs": int(program_summary["fixed_ci_positive_both_scores"].sum()),
        "post_lock_random_dual_score_support_programs": int(program_summary["random_reml_ci_positive_both_scores"].sum()),
        "formal_replication_eligible_programs": 0,
        "matrix_sha256": hashlib.sha256(matrix_path.read_bytes()).hexdigest().upper(),
        "interpretation": "post-lock sensitivity synthesis; frozen Phase 2 labels unchanged",
    }
    (repo / "results" / "logs" / "phase3_cross_cohort_synthesis_run.json").write_text(json.dumps(summary, indent=2) + "\n", encoding="utf-8")
    print(concordance.to_string(index=False))
    print(program_summary.to_string(index=False))
    print(json.dumps(summary, indent=2))


if __name__ == "__main__":
    main()
