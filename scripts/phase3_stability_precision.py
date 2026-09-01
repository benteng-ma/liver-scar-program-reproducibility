from __future__ import annotations

import hashlib
import json
from itertools import combinations
from pathlib import Path

import numpy as np
import pandas as pd
from scipy.cluster.hierarchy import leaves_list, linkage
from scipy.spatial.distance import squareform
from scipy.stats import norm, rankdata
from statsmodels.stats.power import TTestIndPower

from analyze_gse202379_programs import stable_seed
from audit_gse202379_gates import contrast_groups as gse202379_contrast_groups
from audit_gse256398_gates import CONTRASTS as GSE256398_CONTRASTS


N_BOOTSTRAPS = 10_000
N_WEIGHT_DRAWS = 100_000


def overlap_and_effective_tests(repo: Path, output: Path) -> tuple[pd.DataFrame, pd.DataFrame, pd.DataFrame]:
    inventory = pd.read_csv(repo / "literature" / "program_inventory.csv")
    sets = {
        program_id: set(rows["gene_symbol"].dropna().astype(str).str.upper())
        for program_id, rows in inventory.groupby("program_id")
    }
    lineages = inventory.groupby("program_id")["cell_lineage"].first().to_dict()
    rows: list[dict[str, object]] = []
    order_rows: list[dict[str, object]] = []
    for lineage in sorted(set(lineages.values())):
        programs = sorted(program for program, value in lineages.items() if value == lineage)
        distance = np.zeros((len(programs), len(programs)), dtype=float)
        for i, left in enumerate(programs):
            for j, right in enumerate(programs):
                shared = len(sets[left] & sets[right])
                union = len(sets[left] | sets[right])
                smaller = min(len(sets[left]), len(sets[right]))
                jaccard = shared / union if union else np.nan
                overlap_coefficient = shared / smaller if smaller else np.nan
                distance[i, j] = 1 - jaccard
                if i <= j:
                    rows.append(
                        {
                            "lineage": lineage,
                            "program_left": left,
                            "program_right": right,
                            "genes_left": len(sets[left]),
                            "genes_right": len(sets[right]),
                            "shared_genes": shared,
                            "jaccard": jaccard,
                            "overlap_coefficient": overlap_coefficient,
                        }
                    )
        if len(programs) > 1:
            tree = linkage(squareform(distance, checks=False), method="average")
            ordered = [programs[index] for index in leaves_list(tree)]
        else:
            ordered = programs
        order_rows.extend(
            {"lineage": lineage, "cluster_order": index + 1, "program_id": program}
            for index, program in enumerate(ordered)
        )
    overlaps = pd.DataFrame(rows)
    orders = pd.DataFrame(order_rows)
    overlaps.to_csv(output / "program_gene_overlap.csv", index=False)
    orders.to_csv(output / "program_overlap_cluster_order.csv", index=False)

    effective_rows: list[dict[str, object]] = []
    for path in sorted((repo / "results").glob("**/*donor_program_scores.csv.gz")):
        scores = pd.read_csv(path)
        dataset = str(scores["dataset_id"].iloc[0])
        eligibility = "eligible_sensitivity_20" if "eligible_sensitivity_20" in scores else "eligible_20"
        donor_column = "canonical_donor_id" if "canonical_donor_id" in scores else "donor_id"
        scores = scores[scores[eligibility].astype(bool)].copy()
        for (lineage, method), values in scores.groupby(["lineage", "score_method"], sort=True):
            pivot = values.pivot_table(index=donor_column, columns="program_id", values="score", aggfunc="first").dropna(axis=0)
            if pivot.shape[0] < 3 or pivot.shape[1] < 2:
                continue
            correlation = pivot.corr().fillna(0).to_numpy()
            np.fill_diagonal(correlation, 1)
            eigenvalues = np.clip(np.linalg.eigvalsh(correlation), 0, None)
            effective = float(eigenvalues.sum() ** 2 / np.sum(eigenvalues**2))
            upper = np.abs(correlation[np.triu_indices_from(correlation, k=1)])
            effective_rows.append(
                {
                    "dataset_id": dataset,
                    "score_method": method,
                    "lineage": lineage,
                    "n_donors": pivot.shape[0],
                    "n_programs": pivot.shape[1],
                    "effective_program_tests_participation_ratio": effective,
                    "effective_fraction": effective / pivot.shape[1],
                    "median_absolute_score_correlation": float(np.median(upper)),
                    "maximum_absolute_score_correlation": float(np.max(upper)),
                    "eigenvalues_descending": ";".join(f"{value:.8g}" for value in eigenvalues[::-1]),
                }
            )
    effective = pd.DataFrame(effective_rows)
    effective.to_csv(output / "program_effective_test_counts.csv", index=False)
    return overlaps, orders, effective


def bootstrap_g(score_matrix: np.ndarray, n_control: int, n_case: int, seed: int) -> np.ndarray:
    rng = np.random.default_rng(seed)
    control = score_matrix[:, :n_control]
    case = score_matrix[:, n_control:]
    control_indices = rng.integers(0, n_control, size=(N_BOOTSTRAPS, n_control))
    case_indices = rng.integers(0, n_case, size=(N_BOOTSTRAPS, n_case))
    sampled_control = control[:, control_indices]
    sampled_case = case[:, case_indices]
    degrees_freedom = n_control + n_case - 2
    pooled = (
        (n_control - 1) * sampled_control.var(axis=2, ddof=1)
        + (n_case - 1) * sampled_case.var(axis=2, ddof=1)
    ) / degrees_freedom
    correction = 1 - 3 / (4 * degrees_freedom - 1)
    with np.errstate(divide="ignore", invalid="ignore"):
        effects = correction * (sampled_case.mean(axis=2) - sampled_control.mean(axis=2)) / np.sqrt(pooled)
    effects[~np.isfinite(effects)] = np.nan
    return effects


def bootstrap_contexts(repo: Path) -> list[dict[str, object]]:
    contexts: list[dict[str, object]] = []
    scores = pd.read_csv(repo / "results" / "primary" / "gse202379_donor_program_scores.csv.gz")
    manifest = pd.read_csv(repo / "data" / "interim" / "GSE202379" / "donor_lineage_manifest.csv")
    groups = gse202379_contrast_groups(manifest)
    for contrast in ("clinical_cirrhosis_vs_healthy", "advanced_f3f4_vs_f0_non_end_stage"):
        mapping = pd.Series(groups[contrast].to_numpy(), index=manifest["group_id"]).to_dict()
        contexts.append({"dataset_id": "GSE202379", "contrast": contrast, "scores": scores, "group": scores["group_id"].map(mapping).fillna("excluded"), "eligibility": "eligible_primary_30", "donor": "canonical_donor_id"})
    scores = pd.read_csv(repo / "results" / "sensitivity" / "gse244832_donor_program_scores.csv.gz")
    contexts.append({"dataset_id": "GSE244832", "contrast": "mash_f2f4_group_vs_normal_sensitivity", "scores": scores, "group": scores["disease_group"].map({"normal": "control", "MASH": "case"}).fillna("excluded"), "eligibility": "eligible_30", "donor": "donor_id"})
    scores = pd.read_csv(repo / "results" / "phase3" / "gse256398_donor_program_scores.csv.gz")
    for contrast, (controls, cases) in GSE256398_CONTRASTS.items():
        groups = scores["disease_group"].map({**{value: "control" for value in controls}, **{value: "case" for value in cases}}).fillna("excluded")
        contexts.append({"dataset_id": "GSE256398_human", "contrast": contrast, "scores": scores, "group": groups, "eligibility": "eligible_30", "donor": "donor_id"})
    return contexts


def donor_bootstrap(repo: Path, output: Path) -> tuple[pd.DataFrame, pd.DataFrame]:
    summary_rows: list[dict[str, object]] = []
    pair_rows: list[dict[str, object]] = []
    for context in bootstrap_contexts(repo):
        scores = context["scores"].copy()
        scores["comparison_group"] = context["group"].to_numpy()
        scores = scores[scores[context["eligibility"]].astype(bool) & scores["comparison_group"].ne("excluded")]
        for (lineage, method), values in scores.groupby(["lineage", "score_method"], sort=True):
            donor_groups = values[[context["donor"], "comparison_group"]].drop_duplicates().sort_values(["comparison_group", context["donor"]])
            controls = donor_groups.loc[donor_groups["comparison_group"].eq("control"), context["donor"]].tolist()
            cases = donor_groups.loc[donor_groups["comparison_group"].eq("case"), context["donor"]].tolist()
            if len(controls) < 3 or len(cases) < 3:
                continue
            donor_order = controls + cases
            pivot = values.pivot_table(index="program_id", columns=context["donor"], values="score", aggfunc="first").reindex(columns=donor_order).dropna(axis=0)
            effects = bootstrap_g(pivot.to_numpy(dtype=float), len(controls), len(cases), stable_seed("phase3_bootstrap", context["dataset_id"], context["contrast"], lineage, method))
            ranks = rankdata(-effects, axis=0, method="average", nan_policy="omit")
            programs = pivot.index.tolist()
            for index, program in enumerate(programs):
                valid_effect = effects[index, np.isfinite(effects[index])]
                valid_rank = ranks[index, np.isfinite(ranks[index])]
                summary_rows.append(
                    {
                        "dataset_id": context["dataset_id"],
                        "contrast": context["contrast"],
                        "lineage": lineage,
                        "score_method": method,
                        "program_id": program,
                        "n_control": len(controls),
                        "n_case": len(cases),
                        "valid_bootstraps": len(valid_effect),
                        "positive_effect_probability": float((valid_effect > 0).mean()),
                        "effect_median": float(np.median(valid_effect)),
                        "effect_ci95_low": float(np.quantile(valid_effect, 0.025)),
                        "effect_ci95_high": float(np.quantile(valid_effect, 0.975)),
                        "rank_median": float(np.median(valid_rank)),
                        "rank_ci95_low": float(np.quantile(valid_rank, 0.025)),
                        "rank_ci95_high": float(np.quantile(valid_rank, 0.975)),
                        "top_five_probability": float((valid_rank <= 5).mean()),
                    }
                )
            for left_index, right_index in combinations(range(len(programs)), 2):
                valid = np.isfinite(effects[left_index]) & np.isfinite(effects[right_index])
                pair_rows.append(
                    {
                        "dataset_id": context["dataset_id"],
                        "contrast": context["contrast"],
                        "lineage": lineage,
                        "score_method": method,
                        "program_left": programs[left_index],
                        "program_right": programs[right_index],
                        "valid_bootstraps": int(valid.sum()),
                        "probability_left_effect_greater": float((effects[left_index, valid] > effects[right_index, valid]).mean()),
                    }
                )
    summary = pd.DataFrame(summary_rows)
    pairs = pd.DataFrame(pair_rows)
    summary.to_csv(output / "donor_bootstrap_program_rank_stability.csv", index=False)
    pairs.to_csv(output / "donor_bootstrap_pairwise_ordering.csv.gz", index=False, compression="gzip")
    return summary, pairs


def weight_sensitivity(repo: Path, output: Path) -> tuple[pd.DataFrame, pd.DataFrame]:
    card = pd.read_csv(repo / "results" / "deep_benchmark" / "program_transportability_report_card.csv")
    domains = [
        "measurement_domain_0_20",
        "score_method_domain_0_20",
        "directional_transfer_domain_0_20",
        "matched_random_specificity_domain_0_20",
        "endpoint_evidence_domain_0_20",
    ]
    rng = np.random.default_rng(stable_seed("phase3_report_card_weights"))
    weights = rng.dirichlet(np.ones(len(domains)), size=N_WEIGHT_DRAWS)
    totals = card[domains].to_numpy(dtype=float) / 20 @ weights.T * 100
    ranks = rankdata(-totals, axis=0, method="average")
    rows: list[dict[str, object]] = []
    for index, program in enumerate(card["program_id"]):
        rows.append(
            {
                "program_id": program,
                "lineage": card.loc[index, "lineage"],
                "frozen_equal_weight_total": card.loc[index, "transportability_readiness_total_0_100"],
                "rank_median": float(np.median(ranks[index])),
                "rank_ci95_low": float(np.quantile(ranks[index], 0.025)),
                "rank_ci95_high": float(np.quantile(ranks[index], 0.975)),
                "top_one_probability": float((ranks[index] == 1).mean()),
                "top_three_probability": float((ranks[index] <= 3).mean()),
                "top_five_probability": float((ranks[index] <= 5).mean()),
                "score_median": float(np.median(totals[index])),
                "score_ci95_low": float(np.quantile(totals[index], 0.025)),
                "score_ci95_high": float(np.quantile(totals[index], 0.975)),
                "dirichlet_alpha_each_domain": 1.0,
                "weight_draws": N_WEIGHT_DRAWS,
            }
        )
    result = pd.DataFrame(rows).sort_values(["top_five_probability", "rank_median"], ascending=[False, True])
    result.to_csv(output / "report_card_weight_sensitivity.csv", index=False)
    weight_summary = pd.DataFrame(
        {
            "domain": domains,
            "mean_weight": weights.mean(axis=0),
            "q025_weight": np.quantile(weights, 0.025, axis=0),
            "median_weight": np.quantile(weights, 0.5, axis=0),
            "q975_weight": np.quantile(weights, 0.975, axis=0),
        }
    )
    weight_summary.to_csv(output / "report_card_weight_draw_summary.csv", index=False)
    return result, weight_summary


def combined_effects(repo: Path) -> pd.DataFrame:
    effects = pd.read_csv(repo / "results" / "meta" / "cross_cohort_effect_matrix.csv")
    phase3 = pd.read_csv(repo / "results" / "phase3" / "gse256398_program_effects.csv")
    phase3["formal_primary_row"] = False
    return pd.concat([effects, phase3], ignore_index=True, sort=False)


def threshold_grid(repo: Path, output: Path) -> tuple[pd.DataFrame, pd.DataFrame]:
    effects = combined_effects(repo)
    random = pd.concat(
        [pd.read_csv(path) for path in sorted((repo / "results" / "random_controls").glob("*_random_module_benchmark.csv"))],
        ignore_index=True,
        sort=False,
    )
    keys = ["dataset_id", "contrast", "lineage", "program_id", "score_method"]
    merged = effects.merge(
        random[keys + ["real_effect_percentile"]].drop_duplicates(keys),
        on=keys,
        how="inner",
    )
    summary_rows: list[dict[str, object]] = []
    detail_rows: list[dict[str, object]] = []
    for cell_gate in (20, 30):
        for coverage in (0.60, 0.80, 0.90):
            for random_percentile in (0.90, 0.95, 0.99):
                subset = merged[merged["cell_gate"].ge(cell_gate) & merged["program_coverage"].ge(coverage)].copy()
                subset["positive"] = subset["hedges_g"].gt(0)
                subset["random_specific"] = subset["real_effect_percentile"].ge(random_percentile)
                subset["positive_ci"] = subset["robust_ci95_low"].gt(0)
                contexts = []
                for identifiers, rows in subset.groupby(["dataset_id", "contrast", "lineage", "program_id"], sort=False):
                    if set(rows["score_method"]) != {"singscore", "standardized_mean"}:
                        continue
                    dual_positive_random = bool((rows["positive"] & rows["random_specific"]).all())
                    dual_positive_ci_random = bool((rows["positive"] & rows["positive_ci"] & rows["random_specific"]).all())
                    contexts.append((identifiers, dual_positive_random, dual_positive_ci_random))
                    detail_rows.append(
                        {
                            "cell_gate": cell_gate,
                            "coverage_threshold": coverage,
                            "random_percentile_threshold": random_percentile,
                            "dataset_id": identifiers[0],
                            "contrast": identifiers[1],
                            "lineage": identifiers[2],
                            "program_id": identifiers[3],
                            "dual_score_positive_and_random_specific": dual_positive_random,
                            "dual_score_positive_ci_and_random_specific": dual_positive_ci_random,
                        }
                    )
                summary_rows.append(
                    {
                        "cell_gate": cell_gate,
                        "coverage_threshold": coverage,
                        "random_percentile_threshold": random_percentile,
                        "evaluable_program_contexts": len(contexts),
                        "dual_score_positive_random_contexts": sum(value[1] for value in contexts),
                        "dual_score_positive_ci_random_contexts": sum(value[2] for value in contexts),
                        "programs_with_any_dual_positive_random_context": len({value[0][3] for value in contexts if value[1]}),
                    }
                )
    summary = pd.DataFrame(summary_rows)
    details = pd.DataFrame(detail_rows)
    summary.to_csv(output / "threshold_sensitivity_grid.csv", index=False)
    details.to_csv(output / "threshold_sensitivity_contexts.csv.gz", index=False, compression="gzip")
    return summary, details


def precision_and_planning(repo: Path, output: Path, bootstrap: pd.DataFrame) -> tuple[pd.DataFrame, pd.DataFrame, pd.DataFrame]:
    effects = combined_effects(repo)
    power_model = TTestIndPower()
    normal_multiplier = float(norm.ppf(0.975) + norm.ppf(0.80))

    def approximate_balanced_n(effect_size: float) -> float:
        if not np.isfinite(effect_size) or effect_size < 0.10:
            return np.nan
        return float(max(3, np.ceil(2 * normal_multiplier**2 / effect_size**2)))

    rows: list[dict[str, object]] = []
    for effect in effects.itertuples(index=False):
        n_control = int(effect.n_control)
        n_case = int(effect.n_case)
        if n_control < 3 or n_case < 3:
            continue
        mde = normal_multiplier * np.sqrt(1 / n_control + 1 / n_case)
        observed = abs(float(effect.hedges_g))
        planned = approximate_balanced_n(observed)
        rows.append(
            {
                "dataset_id": effect.dataset_id,
                "contrast": effect.contrast,
                "lineage": effect.lineage,
                "program_id": effect.program_id,
                "score_method": effect.score_method,
                "n_control": n_control,
                "n_case": n_case,
                "hedges_g": effect.hedges_g,
                "robust_ci95_width": float(effect.robust_ci95_high - effect.robust_ci95_low) if pd.notna(effect.robust_ci95_high) else np.nan,
                "minimum_detectable_cohens_d_80_power": mde,
                "balanced_n_per_group_for_observed_absolute_g_80_power": planned,
            }
        )
    precision = pd.DataFrame(rows)
    precision.to_csv(output / "contrast_precision_and_mde.csv", index=False)

    curve_rows: list[dict[str, object]] = []
    for effect_size in (0.3, 0.5, 0.8, 1.0, 1.5):
        for n_per_group in range(3, 61):
            curve_rows.append(
                {
                    "standardized_effect": effect_size,
                    "balanced_n_per_group": n_per_group,
                    "two_sided_alpha": 0.05,
                    "power": float(power_model.power(effect_size=effect_size, nobs1=n_per_group, alpha=0.05, ratio=1.0, alternative="two-sided")),
                }
            )
    curves = pd.DataFrame(curve_rows)
    curves.to_csv(output / "balanced_sample_size_power_curves.csv", index=False)

    selected = bootstrap[
        bootstrap["program_id"].isin(["RAM2019_MAC_SIG_A_SAM", "RAM2019_MAC_SIG_B_SAM", "RAM2019_MAC_SIG_E_TMO"])
        | bootstrap["lineage"].eq("endothelial")
    ].copy()
    selected["balanced_n_per_group_for_bootstrap_median_effect_80_power"] = selected["effect_median"].abs().map(
        approximate_balanced_n
    )
    selected["planning_interpretation"] = "bootstrap estimate for prospective planning; not validation evidence"
    selected.to_csv(output / "selected_program_sample_size_planning.csv", index=False)
    return precision, curves, selected


def main() -> None:
    repo = Path(__file__).resolve().parents[1]
    output = repo / "results" / "phase3"
    output.mkdir(parents=True, exist_ok=True)
    overlaps, orders, effective = overlap_and_effective_tests(repo, output)
    bootstrap, pairs = donor_bootstrap(repo, output)
    weights, _ = weight_sensitivity(repo, output)
    grid, details = threshold_grid(repo, output)
    precision, curves, planning = precision_and_planning(repo, output, bootstrap)
    summary = {
        "overlap_rows": len(overlaps),
        "effective_test_rows": len(effective),
        "bootstrap_program_rows": len(bootstrap),
        "bootstrap_pairwise_rows": len(pairs),
        "bootstrap_resamples": N_BOOTSTRAPS,
        "weight_sensitivity_programs": len(weights),
        "dirichlet_weight_draws": N_WEIGHT_DRAWS,
        "threshold_grid_cells": len(grid),
        "precision_rows": len(precision),
        "planning_rows": len(planning),
        "interpretation": "post-lock robustness and planning analyses; frozen Phase 2 labels unchanged",
    }
    (repo / "results" / "logs" / "phase3_stability_precision_run.json").write_text(
        json.dumps(summary, indent=2) + "\n", encoding="utf-8"
    )
    print(json.dumps(summary, indent=2))


if __name__ == "__main__":
    main()
