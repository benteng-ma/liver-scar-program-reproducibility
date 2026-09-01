from __future__ import annotations

import hashlib
import json
from itertools import combinations
from pathlib import Path

import numpy as np
import pandas as pd


SEED = 20260830
N_RESAMPLES = 10_000

REPRESENTATIVE_ENDPOINTS = {
    "GSE202379": "advanced_f3f4_vs_f0_non_end_stage",
    "GSE290642_human": "f4_vs_f0_reconstructed_label_sensitivity",
    "GSE244832": "mash_f2f4_group_vs_normal_sensitivity",
    "GSE210077_Watson6": "mixed_f2f4_fibrosis_vs_healthy_sensitivity",
    "GSE181483_human": "cirrhosis_vs_healthy_directional",
}


def stable_seed(*parts: object) -> int:
    payload = "|".join(str(part) for part in parts).encode("utf-8")
    return int.from_bytes(hashlib.sha256(payload).digest()[:8], "little") % (2**32 - 1)


def bh_adjust(values: pd.Series) -> pd.Series:
    result = pd.Series(np.nan, index=values.index, dtype=float)
    valid = values.dropna().astype(float)
    if valid.empty:
        return result
    order = np.argsort(valid.to_numpy())
    ranked = valid.to_numpy()[order]
    adjusted = ranked * len(ranked) / np.arange(1, len(ranked) + 1)
    adjusted = np.minimum.accumulate(adjusted[::-1])[::-1]
    adjusted = np.clip(adjusted, 0, 1)
    result.loc[valid.index[order]] = adjusted
    return result


def average_ranks(values: np.ndarray) -> np.ndarray:
    values = np.asarray(values, dtype=float)
    order = np.argsort(values, kind="mergesort")
    sorted_values = values[order]
    ranks = np.empty(len(values), dtype=float)
    start = 0
    while start < len(values):
        stop = start + 1
        while stop < len(values) and sorted_values[stop] == sorted_values[start]:
            stop += 1
        ranks[order[start:stop]] = (start + 1 + stop) / 2.0
        start = stop
    return ranks


def spearman_statistic(left: np.ndarray | pd.Series, right: np.ndarray | pd.Series) -> float:
    left_ranks = average_ranks(np.asarray(left, dtype=float))
    right_ranks = average_ranks(np.asarray(right, dtype=float))
    left_centered = left_ranks - left_ranks.mean()
    right_centered = right_ranks - right_ranks.mean()
    denominator = np.sqrt(np.sum(left_centered**2) * np.sum(right_centered**2))
    if denominator == 0:
        return float("nan")
    return float(np.sum(left_centered * right_centered) / denominator)


def hc3_standardized_slope(stage: np.ndarray, score: np.ndarray) -> tuple[float, float, float, float]:
    standardized = (score - score.mean()) / score.std(ddof=1)
    x = np.column_stack([np.ones(len(stage)), stage])
    xtx_inv = np.linalg.inv(x.T @ x)
    beta = xtx_inv @ x.T @ standardized
    residual = standardized - x @ beta
    leverage = np.einsum("ij,jk,ik->i", x, xtx_inv, x)
    scaled = residual / np.maximum(1.0 - leverage, 1e-12)
    meat = x.T @ ((scaled**2)[:, None] * x)
    covariance = xtx_inv @ meat @ xtx_inv
    se = float(np.sqrt(covariance[1, 1]))
    slope = float(beta[1])
    return slope, se, slope - 1.96 * se, slope + 1.96 * se


def stage_trend_row(values: pd.DataFrame, analysis_set: str) -> dict[str, object]:
    stage = values["fibrosis_stage_numeric"].to_numpy(dtype=float)
    score = values["score"].to_numpy(dtype=float)
    observed = spearman_statistic(stage, score)
    rng = np.random.default_rng(
        stable_seed(analysis_set, values["lineage"].iloc[0], values["program_id"].iloc[0], values["score_method"].iloc[0])
    )
    permuted = np.empty(N_RESAMPLES, dtype=float)
    bootstrapped = np.empty(N_RESAMPLES, dtype=float)
    for index in range(N_RESAMPLES):
        permuted[index] = spearman_statistic(rng.permutation(stage), score)
        sample = rng.integers(0, len(stage), size=len(stage))
        bootstrapped[index] = spearman_statistic(stage[sample], score[sample])
    permutation_p = float((1 + np.sum(np.abs(permuted) >= abs(observed))) / (N_RESAMPLES + 1))
    bootstrap_valid = bootstrapped[np.isfinite(bootstrapped)]
    slope, slope_se, slope_low, slope_high = hc3_standardized_slope(stage, score)
    return {
        "analysis_set": analysis_set,
        "lineage": values["lineage"].iloc[0],
        "program_id": values["program_id"].iloc[0],
        "score_method": values["score_method"].iloc[0],
        "n_donors": len(values),
        "n_stages": int(values["fibrosis_stage_numeric"].nunique()),
        "stage_min": int(stage.min()),
        "stage_max": int(stage.max()),
        "spearman_rho": observed,
        "spearman_bootstrap_ci95_low": float(np.quantile(bootstrap_valid, 0.025)),
        "spearman_bootstrap_ci95_high": float(np.quantile(bootstrap_valid, 0.975)),
        "permutation_p_two_sided": permutation_p,
        "standardized_score_slope_per_stage": slope,
        "hc3_se": slope_se,
        "hc3_ci95_low": slope_low,
        "hc3_ci95_high": slope_high,
        "resamples": N_RESAMPLES,
    }


def stage_trends(repo: Path, output: Path) -> tuple[pd.DataFrame, pd.DataFrame]:
    scores = pd.read_csv(repo / "results" / "primary" / "gse202379_donor_program_scores.csv.gz")
    scores["fibrosis_stage_numeric"] = pd.to_numeric(scores["Fibrosis.score..F0.4."], errors="coerce")
    scores = scores[scores["eligible_primary_30"].astype(bool) & scores["fibrosis_stage_numeric"].notna()].copy()
    sets = {
        "non_end_stage_primary_exploratory": scores[~scores["Disease.status"].eq("end stage")].copy(),
        "including_end_stage_sensitivity": scores.copy(),
    }
    rows: list[dict[str, object]] = []
    medians: list[pd.DataFrame] = []
    keys = ["lineage", "program_id", "score_method"]
    for name, frame in sets.items():
        for _, values in frame.groupby(keys, sort=True):
            values = values.drop_duplicates("canonical_donor_id").sort_values("canonical_donor_id")
            if len(values) < 10 or values["fibrosis_stage_numeric"].nunique() < 4:
                continue
            rows.append(stage_trend_row(values, name))
        summary = (
            frame.groupby(keys + ["fibrosis_stage_numeric"], as_index=False)
            .agg(
                n_donors=("canonical_donor_id", "nunique"),
                median_score=("score", "median"),
                q25_score=("score", lambda x: x.quantile(0.25)),
                q75_score=("score", lambda x: x.quantile(0.75)),
            )
        )
        summary.insert(0, "analysis_set", name)
        medians.append(summary)
    trends = pd.DataFrame(rows)
    trends["fdr_bh"] = trends.groupby("analysis_set", group_keys=False)["permutation_p_two_sided"].apply(bh_adjust)
    trends["positive_fdr_005"] = trends["spearman_rho"].gt(0) & trends["fdr_bh"].lt(0.05)
    dual = (
        trends.pivot_table(
            index=["analysis_set", "lineage", "program_id"],
            columns="score_method",
            values="positive_fdr_005",
            aggfunc="first",
            fill_value=False,
        )
        .reindex(columns=["singscore", "standardized_mean"], fill_value=False)
        .all(axis=1)
        .rename("dual_score_positive_fdr_005")
    )
    trends = trends.merge(dual.reset_index(), on=["analysis_set", "lineage", "program_id"], how="left")
    medians_frame = pd.concat(medians, ignore_index=True)
    trends.to_csv(output / "gse202379_stage_trends.csv", index=False)
    medians_frame.to_csv(output / "gse202379_stage_medians.csv", index=False)
    return trends, medians_frame


def representative_effects(repo: Path) -> pd.DataFrame:
    effects = pd.read_csv(repo / "results" / "meta" / "cross_cohort_effect_matrix.csv")
    keep = pd.Series(False, index=effects.index)
    for dataset, contrast in REPRESENTATIVE_ENDPOINTS.items():
        keep |= effects["dataset_id"].eq(dataset) & effects["contrast"].eq(contrast)
    return effects[keep].copy()


def concordance(repo: Path, output: Path) -> tuple[pd.DataFrame, pd.DataFrame]:
    effects = representative_effects(repo)
    within_rows: list[dict[str, object]] = []
    for (dataset, contrast, lineage), values in effects.groupby(["dataset_id", "contrast", "lineage"], sort=True):
        pivot = values.pivot_table(index="program_id", columns="score_method", values="hedges_g", aggfunc="first").dropna()
        if len(pivot) < 3 or not {"singscore", "standardized_mean"}.issubset(pivot.columns):
            continue
        within_rows.append(
            {
                "dataset_id": dataset,
                "contrast": contrast,
                "lineage": lineage,
                "n_programs": len(pivot),
                "spearman_rho": spearman_statistic(pivot["singscore"], pivot["standardized_mean"]),
                "sign_agreement": float((np.sign(pivot["singscore"]) == np.sign(pivot["standardized_mean"])).mean()),
            }
        )
    within = pd.DataFrame(within_rows)
    across_rows: list[dict[str, object]] = []
    for lineage, lineage_values in effects.groupby("lineage", sort=True):
        for score_method, method_values in lineage_values.groupby("score_method", sort=True):
            cohort_values = {
                dataset: values.set_index("program_id")["hedges_g"]
                for dataset, values in method_values.groupby("dataset_id", sort=True)
            }
            for left, right in combinations(sorted(cohort_values), 2):
                paired = pd.concat([cohort_values[left], cohort_values[right]], axis=1, join="inner").dropna()
                paired.columns = ["left", "right"]
                if len(paired) < 3:
                    continue
                across_rows.append(
                    {
                        "lineage": lineage,
                        "score_method": score_method,
                        "dataset_left": left,
                        "dataset_right": right,
                        "n_programs": len(paired),
                        "spearman_rho": spearman_statistic(paired["left"], paired["right"]),
                        "sign_agreement": float((np.sign(paired["left"]) == np.sign(paired["right"])).mean()),
                        "includes_watson": "GSE210077_Watson6" in {left, right},
                    }
                )
    across = pd.DataFrame(across_rows)
    within.to_csv(output / "within_cohort_score_method_concordance.csv", index=False)
    across.to_csv(output / "cross_cohort_program_concordance.csv", index=False)
    return within, across


def evidence_matrix(repo: Path, output: Path) -> pd.DataFrame:
    effects = pd.read_csv(repo / "results" / "meta" / "cross_cohort_effect_matrix.csv")
    classifications = pd.read_csv(repo / "results" / "meta" / "program_classification_table.csv")
    random_frames = [pd.read_csv(path) for path in sorted((repo / "results" / "random_controls").glob("*_random_module_benchmark.csv"))]
    random = pd.concat(random_frames, ignore_index=True, sort=False)
    programs = classifications[["program_id", "lineage"]].copy()

    evaluable = effects.groupby("program_id").size().gt(0).rename("evaluable_independent")
    positive_both = (
        effects.groupby(["program_id", "dataset_id", "score_method"])["hedges_g"].max().gt(0)
        .unstack("score_method", fill_value=False)
        .reindex(columns=["singscore", "standardized_mean"], fill_value=False)
        .all(axis=1)
        .groupby("program_id").any()
        .rename("positive_both_scores_any_cohort")
    )
    random_both = (
        random.groupby(["program_id", "dataset_id", "score_method"])["above_random_95th_percentile"].max().astype(bool)
        .unstack("score_method", fill_value=False)
        .reindex(columns=["singscore", "standardized_mean"], fill_value=False)
        .all(axis=1)
        .groupby("program_id").any()
        .rename("random_specific_both_scores_any_cohort")
    )
    primary = effects[effects["formal_primary_row"].astype(bool)].copy()
    primary_both = (
        primary.groupby(["program_id", "dataset_id", "score_method"])["robust_ci95_low"].max().gt(0)
        .unstack("score_method", fill_value=False)
        .reindex(columns=["singscore", "standardized_mean"], fill_value=False)
        .all(axis=1)
        .groupby("program_id").any()
        .rename("positive_primary_interval_both_scores")
    )
    matrix = programs.set_index("program_id").join([evaluable, positive_both, random_both, primary_both]).fillna(False)
    matrix = matrix.join(
        classifications.set_index("program_id")[[
            "advanced_sensitivity_meta_available",
            "within_cell_state_replicated",
            "pan_cirrhotic_transportable",
            "assay_transfer_class",
        ]]
    )
    matrix["assay_robust"] = matrix["assay_transfer_class"].eq("ASSAY_ROBUST")
    matrix = matrix.drop(columns="assay_transfer_class").reset_index()
    matrix.to_csv(output / "program_evidence_attrition_matrix.csv", index=False)
    return matrix


def main() -> None:
    repo = Path(__file__).resolve().parents[1]
    output = repo / "results" / "exploratory"
    output.mkdir(parents=True, exist_ok=True)
    trends, _ = stage_trends(repo, output)
    within, across = concordance(repo, output)
    evidence = evidence_matrix(repo, output)
    summary = {
        "frozen_plan": "reports/exploratory_value_rescue_plan.md",
        "stage_trend_rows": len(trends),
        "non_end_stage_dual_score_fdr_programs": sorted(
            trends.loc[
                trends["analysis_set"].eq("non_end_stage_primary_exploratory")
                & trends["dual_score_positive_fdr_005"].astype(bool),
                "program_id",
            ].unique()
        ),
        "including_end_stage_dual_score_fdr_programs": sorted(
            trends.loc[
                trends["analysis_set"].eq("including_end_stage_sensitivity")
                & trends["dual_score_positive_fdr_005"].astype(bool),
                "program_id",
            ].unique()
        ),
        "within_cohort_concordance_rows": len(within),
        "within_cohort_median_spearman": float(within["spearman_rho"].median()),
        "cross_cohort_concordance_rows": len(across),
        "cross_cohort_median_spearman": float(across["spearman_rho"].median()),
        "cross_cohort_median_spearman_without_watson": float(across.loc[~across["includes_watson"], "spearman_rho"].median()),
        "evidence_matrix_programs": len(evidence),
        "random_specific_both_scores_programs": int(evidence["random_specific_both_scores_any_cohort"].sum()),
        "positive_primary_interval_both_scores_programs": int(evidence["positive_primary_interval_both_scores"].sum()),
        "interpretation_boundary": "exploratory only; no frozen higher-order classification changed",
    }
    log_dir = repo / "results" / "logs"
    log_dir.mkdir(parents=True, exist_ok=True)
    (log_dir / "exploratory_value_analysis_run.json").write_text(json.dumps(summary, indent=2) + "\n", encoding="utf-8")
    print(json.dumps(summary, indent=2))


if __name__ == "__main__":
    main()
