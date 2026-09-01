from __future__ import annotations

import hashlib
import json
from pathlib import Path

import numpy as np
import pandas as pd
import statsmodels.api as sm
from scipy.stats import rankdata, spearmanr, theilslopes
from statsmodels.stats.outliers_influence import variance_inflation_factor

from analyze_gse202379_programs import stable_seed
from analyze_gse256398_programs import benjamini_hochberg


N_PERMUTATIONS = 10_000
HISTOLOGY = {
    "fibrosis": "Fibrosis.score..F0.4.",
    "steatosis": "Steatosis",
    "ballooning": "Ballooning",
    "inflammation": "Inflammation",
}


def permutation_spearman(x: np.ndarray, y: np.ndarray, seed: int) -> tuple[float, float]:
    x_rank = rankdata(x, method="average")
    y_rank = rankdata(y, method="average")
    x_centered = x_rank - x_rank.mean()
    y_centered = y_rank - y_rank.mean()
    denominator = np.sqrt(np.sum(x_centered**2) * np.sum(y_centered**2))
    observed = float(np.sum(x_centered * y_centered) / denominator)
    rng = np.random.default_rng(seed)
    extreme = 0
    for _ in range(N_PERMUTATIONS):
        permuted = rng.permutation(x_centered)
        statistic = float(np.sum(permuted * y_centered) / denominator)
        extreme += abs(statistic) >= abs(observed) - 1e-15
    return observed, (extreme + 1) / (N_PERMUTATIONS + 1)


def zscore(values: pd.Series) -> pd.Series:
    values = values.astype(float)
    sd = values.std(ddof=1)
    if not np.isfinite(sd) or sd == 0:
        raise ValueError("cannot standardize invariant covariate")
    return (values - values.mean()) / sd


def incremental_r_squared(y: pd.Series, full_x: pd.DataFrame, predictor: str) -> float:
    full = sm.OLS(y, sm.add_constant(full_x, has_constant="add")).fit()
    reduced = sm.OLS(y, sm.add_constant(full_x.drop(columns=predictor), has_constant="add")).fit()
    return float(full.rsquared - reduced.rsquared)


def model_rows(values: pd.DataFrame, analysis_set: str, model_name: str) -> list[dict[str, object]]:
    base_columns = list(HISTOLOGY.values()) + ["Age", "Gender", "score"]
    frame = values[base_columns].copy()
    for column in list(HISTOLOGY.values()) + ["Age", "score"]:
        frame[column] = pd.to_numeric(frame[column], errors="coerce")
    frame["sex_female"] = frame["Gender"].astype(str).str.upper().map({"M": 0.0, "F": 1.0})
    frame = frame.dropna()
    if len(frame) < 12 or frame["sex_female"].nunique() < 2:
        return []
    y = zscore(frame["score"])
    standardized = pd.DataFrame(index=frame.index)
    for name, column in HISTOLOGY.items():
        standardized[name] = zscore(frame[column])
    standardized["age"] = zscore(frame["Age"])
    standardized["sex_female"] = frame["sex_female"].astype(float)
    predictor_sets = (
        {name: [name, "age", "sex_female"] for name in HISTOLOGY}
        if model_name == "minimal_age_sex"
        else {name: list(HISTOLOGY) + ["age", "sex_female"] for name in HISTOLOGY}
    )
    rows: list[dict[str, object]] = []
    for histology_name, predictors in predictor_sets.items():
        x = standardized[predictors]
        fit = sm.OLS(y, sm.add_constant(x, has_constant="add")).fit(cov_type="HC3")
        predictor_index = predictors.index(histology_name)
        vif = float(variance_inflation_factor(x.to_numpy(), predictor_index))
        rows.append(
            {
                "dataset_id": "GSE202379",
                "analysis_set": analysis_set,
                "model": model_name,
                "lineage": values["lineage"].iloc[0],
                "program_id": values["program_id"].iloc[0],
                "score_method": values["score_method"].iloc[0],
                "histology_axis": histology_name,
                "n_donors": len(frame),
                "standardized_beta": float(fit.params[histology_name]),
                "hc3_se": float(fit.bse[histology_name]),
                "hc3_ci95_low": float(fit.conf_int().loc[histology_name, 0]),
                "hc3_ci95_high": float(fit.conf_int().loc[histology_name, 1]),
                "hc3_p_two_sided": float(fit.pvalues[histology_name]),
                "model_r_squared": float(fit.rsquared),
                "incremental_r_squared": incremental_r_squared(y, x, histology_name),
                "vif": vif,
            }
        )
    return rows


def gse202379(repo: Path, output: Path) -> tuple[pd.DataFrame, pd.DataFrame]:
    scores = pd.read_csv(repo / "results" / "primary" / "gse202379_donor_program_scores.csv.gz")
    scores = scores[scores["eligible_primary_30"].astype(bool)].copy()
    analysis_sets = {
        "non_end_stage_primary_exploratory": scores[~scores["Disease.status"].eq("end stage")].copy(),
        "including_end_stage_sensitivity": scores.copy(),
    }
    univariate_rows: list[dict[str, object]] = []
    model_output: list[dict[str, object]] = []
    keys = ["lineage", "program_id", "score_method"]
    for analysis_set, frame in analysis_sets.items():
        for _, values in frame.groupby(keys, sort=True):
            values = values.drop_duplicates("canonical_donor_id")
            for histology_name, column in HISTOLOGY.items():
                x = pd.to_numeric(values[column], errors="coerce")
                y = pd.to_numeric(values["score"], errors="coerce")
                valid = x.notna() & y.notna()
                if valid.sum() < 10 or x[valid].nunique() < 3:
                    continue
                rho, p_value = permutation_spearman(
                    x[valid].to_numpy(dtype=float),
                    y[valid].to_numpy(dtype=float),
                    stable_seed("phase3_histology", analysis_set, values["lineage"].iloc[0], values["program_id"].iloc[0], values["score_method"].iloc[0], histology_name),
                )
                univariate_rows.append(
                    {
                        "dataset_id": "GSE202379",
                        "analysis_set": analysis_set,
                        "lineage": values["lineage"].iloc[0],
                        "program_id": values["program_id"].iloc[0],
                        "score_method": values["score_method"].iloc[0],
                        "histology_axis": histology_name,
                        "n_donors": int(valid.sum()),
                        "n_levels": int(x[valid].nunique()),
                        "spearman_rho": rho,
                        "permutation_p_two_sided": p_value,
                        "permutations": N_PERMUTATIONS,
                    }
                )
            model_output.extend(model_rows(values, analysis_set, "minimal_age_sex"))
            model_output.extend(model_rows(values, analysis_set, "joint_histology_age_sex"))
    univariate = pd.DataFrame(univariate_rows)
    univariate["fdr_within_set_lineage_method_axis"] = univariate.groupby(
        ["analysis_set", "lineage", "score_method", "histology_axis"], group_keys=False
    )["permutation_p_two_sided"].apply(benjamini_hochberg)
    models = pd.DataFrame(model_output)
    models["fdr_within_set_model_lineage_method_axis"] = models.groupby(
        ["analysis_set", "model", "lineage", "score_method", "histology_axis"], group_keys=False
    )["hc3_p_two_sided"].apply(benjamini_hochberg)
    univariate.to_csv(output / "gse202379_histology_univariate.csv", index=False)
    models.to_csv(output / "gse202379_histology_models.csv", index=False)
    return univariate, models


def gse244832(repo: Path, output: Path) -> tuple[pd.DataFrame, pd.DataFrame]:
    scores = pd.read_csv(repo / "results" / "sensitivity" / "gse244832_donor_program_scores.csv.gz")
    scores = scores[scores["eligible_30"].astype(bool)].copy()
    scores["metabolic_order"] = scores["disease_group"].map({"normal": 0, "MASL": 1, "MASH": 2})
    rows: list[dict[str, object]] = []
    keys = ["lineage", "program_id", "score_method"]
    for _, values in scores.groupby(keys, sort=True):
        values = values.drop_duplicates("donor_id").dropna(subset=["metabolic_order", "score"])
        order = values["metabolic_order"].to_numpy(dtype=float)
        score = values["score"].to_numpy(dtype=float)
        if len(values) < 12 or len(np.unique(order)) != 3:
            continue
        rho, p_value = permutation_spearman(
            order,
            score,
            stable_seed("phase3_gse244832_progression", values["lineage"].iloc[0], values["program_id"].iloc[0], values["score_method"].iloc[0]),
        )
        slope, intercept, low, high = theilslopes(score, order, alpha=0.95)
        rows.append(
            {
                "dataset_id": "GSE244832",
                "trend": "NORMAL_to_MASL_to_MASH",
                "lineage": values["lineage"].iloc[0],
                "program_id": values["program_id"].iloc[0],
                "score_method": values["score_method"].iloc[0],
                "n_donors": len(values),
                "n_normal": int((order == 0).sum()),
                "n_masl": int((order == 1).sum()),
                "n_mash": int((order == 2).sum()),
                "spearman_rho": rho,
                "permutation_p_two_sided": p_value,
                "permutations": N_PERMUTATIONS,
                "theil_sen_slope_per_level": float(slope),
                "theil_sen_intercept": float(intercept),
                "theil_sen_ci95_low": float(low),
                "theil_sen_ci95_high": float(high),
            }
        )
    trends = pd.DataFrame(rows)
    trends["fdr_within_lineage_method"] = trends.groupby(
        ["lineage", "score_method"], group_keys=False
    )["permutation_p_two_sided"].apply(benjamini_hochberg)
    medians = (
        scores.groupby(["lineage", "program_id", "score_method", "disease_group", "metabolic_order"], as_index=False)
        .agg(n_donors=("donor_id", "nunique"), median_score=("score", "median"), q25=("score", lambda x: x.quantile(0.25)), q75=("score", lambda x: x.quantile(0.75)))
    )
    trends.to_csv(output / "gse244832_metabolic_progression.csv", index=False)
    medians.to_csv(output / "gse244832_metabolic_progression_medians.csv", index=False)
    return trends, medians


def main() -> None:
    repo = Path(__file__).resolve().parents[1]
    output = repo / "results" / "phase3"
    output.mkdir(parents=True, exist_ok=True)
    univariate, models = gse202379(repo, output)
    trends, _ = gse244832(repo, output)
    summary = {
        "gse202379_univariate_rows": len(univariate),
        "gse202379_univariate_fdr_lt_005": int((univariate["fdr_within_set_lineage_method_axis"] < 0.05).sum()),
        "gse202379_model_rows": len(models),
        "gse202379_model_fdr_lt_005": int((models["fdr_within_set_model_lineage_method_axis"] < 0.05).sum()),
        "gse244832_progression_rows": len(trends),
        "gse244832_progression_fdr_lt_005": int((trends["fdr_within_lineage_method"] < 0.05).sum()),
        "permutations_per_row": N_PERMUTATIONS,
        "interpretation": "post-lock specificity and progression analyses; descriptive, non-causal",
    }
    (repo / "results" / "logs" / "phase3_clinical_specificity_run.json").write_text(
        json.dumps(summary, indent=2) + "\n", encoding="utf-8"
    )
    print(json.dumps(summary, indent=2))


if __name__ == "__main__":
    main()
