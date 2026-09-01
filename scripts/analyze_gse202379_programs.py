from __future__ import annotations

import hashlib
import itertools
import json
import math
from pathlib import Path

import numpy as np
import pandas as pd
import statsmodels.api as sm
from scipy.io import mmread
from scipy.stats import rankdata

from audit_gse202379_gates import contrast_groups


SEED = 20260830
MAX_EXACT_PERMUTATIONS = 100_000
MONTE_CARLO_PERMUTATIONS = 10_000


def stable_seed(*parts: str) -> int:
    payload = "|||".join(parts).encode("utf-8")
    offset = int.from_bytes(hashlib.sha256(payload).digest()[:8], "little")
    return (SEED + offset) % (2**32)


def singscore_up(ranks: np.ndarray, gene_indices: np.ndarray) -> np.ndarray:
    n_genes = ranks.shape[0]
    module_size = len(gene_indices)
    observed_mean = ranks[gene_indices, :].mean(axis=0)
    theoretical_min = (module_size + 1) / 2
    theoretical_max = (2 * n_genes - module_size + 1) / 2
    return 2 * (observed_mean - theoretical_min) / (
        theoretical_max - theoretical_min
    ) - 1


def effect_statistics(
    control: np.ndarray,
    case: np.ndarray,
    permutation_seed: int,
) -> dict[str, object]:
    control = np.asarray(control, dtype=float)
    case = np.asarray(case, dtype=float)
    n_control = len(control)
    n_case = len(case)
    degrees_freedom = n_control + n_case - 2
    variance_control = control.var(ddof=1)
    variance_case = case.var(ddof=1)
    pooled_variance = (
        (n_control - 1) * variance_control + (n_case - 1) * variance_case
    ) / degrees_freedom
    if pooled_variance <= 0 or not np.isfinite(pooled_variance):
        raise ValueError("pooled donor-level score variance is not positive")
    pooled_sd = math.sqrt(pooled_variance)
    mean_difference = float(case.mean() - control.mean())
    cohens_d = mean_difference / pooled_sd
    correction = 1 - 3 / (4 * degrees_freedom - 1)
    hedges_g = correction * cohens_d

    values = np.concatenate([control, case])
    labels = np.concatenate(
        [np.zeros(n_control, dtype=float), np.ones(n_case, dtype=float)]
    )
    design = sm.add_constant(labels)
    robust_fit = sm.OLS(values, design).fit(cov_type="HC3")
    robust_se_g = correction * float(robust_fit.bse[1]) / pooled_sd
    robust_low = hedges_g - 1.96 * robust_se_g
    robust_high = hedges_g + 1.96 * robust_se_g
    conventional_se = math.sqrt(
        (n_control + n_case) / (n_control * n_case)
        + hedges_g**2 / (2 * degrees_freedom)
    )
    conventional_low = hedges_g - 1.96 * conventional_se
    conventional_high = hedges_g + 1.96 * conventional_se

    observed_difference = mean_difference
    n_total = n_control + n_case
    allocations = math.comb(n_total, n_case)
    if allocations <= MAX_EXACT_PERMUTATIONS:
        extreme = 0
        total_sum = values.sum()
        for case_indices in itertools.combinations(range(n_total), n_case):
            case_sum = values[list(case_indices)].sum()
            difference = case_sum / n_case - (total_sum - case_sum) / n_control
            extreme += abs(difference) >= abs(observed_difference) - 1e-15
        permutation_p = extreme / allocations
        permutation_mode = f"exact_{allocations}_allocations"
    else:
        rng = np.random.default_rng(permutation_seed)
        extreme = 0
        for _ in range(MONTE_CARLO_PERMUTATIONS):
            permuted = rng.permutation(labels)
            difference = values[permuted == 1].mean() - values[permuted == 0].mean()
            extreme += abs(difference) >= abs(observed_difference) - 1e-15
        permutation_p = (extreme + 1) / (MONTE_CARLO_PERMUTATIONS + 1)
        permutation_mode = f"monte_carlo_{MONTE_CARLO_PERMUTATIONS}"

    return {
        "n_control": n_control,
        "n_case": n_case,
        "mean_control": float(control.mean()),
        "mean_case": float(case.mean()),
        "mean_difference": mean_difference,
        "pooled_sd": pooled_sd,
        "cohens_d": cohens_d,
        "hedges_g": hedges_g,
        "robust_se_g_hc3": robust_se_g,
        "robust_ci95_low": robust_low,
        "robust_ci95_high": robust_high,
        "conventional_se_g": conventional_se,
        "conventional_ci95_low": conventional_low,
        "conventional_ci95_high": conventional_high,
        "permutation_p_two_sided": permutation_p,
        "permutation_mode": permutation_mode,
    }


def main() -> None:
    repo = Path(__file__).resolve().parents[1]
    interim = repo / "data" / "interim" / "GSE202379"
    genes = pd.read_csv(interim / "genes.csv")["gene"].astype(str)
    manifest = pd.read_csv(interim / "donor_lineage_manifest.csv")
    manifest["canonical_donor_id"] = manifest["donor_id"].map(
        lambda value: str(value) if str(value).startswith("P") else f"P{value}"
    )
    manifest["eligible_primary_30"] = manifest["n_cells"].ge(30)
    manifest["eligible_sensitivity_20"] = manifest["n_cells"].ge(20)
    counts = mmread(interim / "donor_lineage_raw_counts.mtx").tocsc()
    if counts.shape != (len(genes), len(manifest)):
        raise RuntimeError("count matrix does not match gene/group manifests")
    library_sizes = np.asarray(counts.sum(axis=0)).ravel()
    if (library_sizes <= 0).any():
        raise RuntimeError("donor-lineage pseudobulk with zero library size")
    log_cpm = np.log2(counts.toarray() / library_sizes * 1_000_000 + 1)
    ranks = rankdata(log_cpm, axis=0, method="average")
    gene_to_index = {gene.upper(): index for index, gene in enumerate(genes)}

    programs = pd.read_csv(repo / "literature" / "program_inventory.csv")
    coverage = pd.read_csv(repo / "results" / "qc" / "gse202379_program_coverage.csv")
    coverage_lookup = coverage.set_index("program_id").to_dict("index")
    score_rows: list[dict[str, object]] = []
    score_arrays: dict[tuple[str, str], np.ndarray] = {}
    reference_stats: dict[str, dict[str, np.ndarray]] = {}

    for lineage, lineage_manifest in manifest.groupby("harmonized_lineage", sort=False):
        columns = lineage_manifest.index.to_numpy()
        reference_columns = lineage_manifest.index[
            lineage_manifest["eligible_sensitivity_20"]
        ].to_numpy()
        means = log_cpm[:, reference_columns].mean(axis=1)
        standard_deviations = log_cpm[:, reference_columns].std(axis=1, ddof=1)
        invariant = standard_deviations == 0
        safe_standard_deviations = standard_deviations.copy()
        safe_standard_deviations[invariant] = 1
        standardized = (
            log_cpm[:, columns] - means[:, None]
        ) / safe_standard_deviations[:, None]
        standardized[invariant, :] = 0
        reference_stats[lineage] = {
            "mean_log2_cpm": means,
            "sd_log2_cpm": standard_deviations,
            "detection_rate": np.asarray(
                (counts[:, reference_columns] > 0).mean(axis=1)
            ).ravel(),
        }
        for program_id, rows in programs[programs["cell_lineage"].eq(lineage)].groupby(
            "program_id"
        ):
            coverage_info = coverage_lookup[program_id]
            if coverage_info["coverage_tier"] == "not_evaluated":
                continue
            measured_genes = sorted(
                set(rows["gene_symbol"].astype(str).str.upper()) & set(gene_to_index)
            )
            gene_indices = np.array(
                [gene_to_index[gene] for gene in measured_genes], dtype=int
            )
            sing = singscore_up(ranks[:, columns], gene_indices)
            zmean = standardized[gene_indices, :].mean(axis=0)
            score_arrays[(program_id, "singscore")] = np.full(len(manifest), np.nan)
            score_arrays[(program_id, "standardized_mean")] = np.full(
                len(manifest), np.nan
            )
            score_arrays[(program_id, "singscore")][columns] = sing
            score_arrays[(program_id, "standardized_mean")][columns] = zmean
            for local_index, manifest_index in enumerate(columns):
                base = manifest.loc[manifest_index].to_dict()
                for score_name, score_value in (
                    ("singscore", sing[local_index]),
                    ("standardized_mean", zmean[local_index]),
                ):
                    score_rows.append(
                        {
                            "dataset_id": "GSE202379",
                            "program_id": program_id,
                            "lineage": lineage,
                            "score_method": score_name,
                            "score": float(score_value),
                            "measured_program_genes": len(measured_genes),
                            "program_coverage": coverage_info["coverage"],
                            "coverage_tier": coverage_info["coverage_tier"],
                            **base,
                            "library_size": int(library_sizes[manifest_index]),
                        }
                    )

    scores = pd.DataFrame(score_rows)
    primary_dir = repo / "results" / "primary"
    sensitivity_dir = repo / "results" / "sensitivity"
    primary_dir.mkdir(parents=True, exist_ok=True)
    sensitivity_dir.mkdir(parents=True, exist_ok=True)
    scores.to_csv(
        primary_dir / "gse202379_donor_program_scores.csv.gz",
        index=False,
        compression="gzip",
    )

    gates = pd.read_csv(repo / "results" / "qc" / "gse202379_donor_gate_summary.csv")
    gate_status = gates.drop_duplicates(["contrast", "lineage"]).set_index(
        ["contrast", "lineage"]
    )
    primary_contrasts = {
        "clinical_cirrhosis_vs_healthy",
        "advanced_f3f4_vs_f0_non_end_stage",
    }
    primary_effect_rows: list[dict[str, object]] = []
    sensitivity_effect_rows: list[dict[str, object]] = []

    for contrast, group in contrast_groups(manifest).items():
        for program_id, program_rows in programs.groupby("program_id"):
            lineage = str(program_rows["cell_lineage"].iloc[0])
            coverage_info = coverage_lookup[program_id]
            if coverage_info["coverage_tier"] == "not_evaluated":
                continue
            status = gate_status.loc[(contrast, lineage)]
            primary_allowed = (
                contrast in primary_contrasts
                and coverage_info["coverage_tier"] == "primary"
                and status["formal_primary_gate"] == "PASS"
            )
            sensitivity_allowed = status["formal_sensitivity_gate"] == "PASS"
            if primary_allowed:
                eligibility_column = "eligible_primary_30"
                analysis_tier = "primary"
            elif sensitivity_allowed:
                eligibility_column = "eligible_sensitivity_20"
                analysis_tier = "sensitivity"
            else:
                continue
            selected = (
                manifest["harmonized_lineage"].eq(lineage)
                & manifest[eligibility_column]
                & group.ne("excluded")
            )
            selected_indices = manifest.index[selected]
            selected_groups = group.loc[selected_indices]
            for score_method in ("singscore", "standardized_mean"):
                values = score_arrays[(program_id, score_method)][selected_indices]
                control = values[selected_groups.eq("control").to_numpy()]
                case = values[selected_groups.eq("case").to_numpy()]
                if len(control) < 3 or len(case) < 3:
                    raise RuntimeError("effect was reached without three donors per group")
                stats = effect_statistics(
                    control,
                    case,
                    stable_seed("GSE202379", contrast, program_id, score_method),
                )
                row = {
                    "dataset_id": "GSE202379",
                    "contrast": contrast,
                    "program_id": program_id,
                    "lineage": lineage,
                    "score_method": score_method,
                    "analysis_tier": analysis_tier,
                    "cell_gate": 30 if primary_allowed else 20,
                    "coverage_tier": coverage_info["coverage_tier"],
                    "program_coverage": coverage_info["coverage"],
                    "expected_direction": "positive_case_higher",
                    **stats,
                }
                if primary_allowed:
                    primary_effect_rows.append(row)
                else:
                    sensitivity_effect_rows.append(row)

    primary_effects = pd.DataFrame(primary_effect_rows).sort_values(
        ["contrast", "lineage", "program_id", "score_method"]
    )
    sensitivity_effects = pd.DataFrame(sensitivity_effect_rows).sort_values(
        ["contrast", "lineage", "program_id", "score_method"]
    )
    primary_effects.to_csv(
        primary_dir / "gse202379_primary_effects.csv", index=False
    )
    sensitivity_effects.to_csv(
        sensitivity_dir / "gse202379_sensitivity_effects.csv", index=False
    )
    run_summary = {
        "dataset_id": "GSE202379",
        "normalization": "log2(CPM + 1) within donor-lineage pseudobulk",
        "score_methods": ["direction-aware singscore", "standardized signed mean"],
        "standardization_reference": "all lineage donors passing 20-cell gate; outcome-blind",
        "primary_effect_rows": len(primary_effects),
        "sensitivity_effect_rows": len(sensitivity_effects),
        "seed": SEED,
        "primary_effects_sha256": hashlib.sha256(
            (primary_dir / "gse202379_primary_effects.csv").read_bytes()
        ).hexdigest().upper(),
        "sensitivity_effects_sha256": hashlib.sha256(
            (sensitivity_dir / "gse202379_sensitivity_effects.csv").read_bytes()
        ).hexdigest().upper(),
    }
    (repo / "results" / "logs" / "gse202379_analysis_run.json").write_text(
        json.dumps(run_summary, indent=2) + "\n", encoding="utf-8"
    )
    print(primary_effects[["contrast", "program_id", "score_method", "n_control", "n_case", "hedges_g", "robust_ci95_low", "robust_ci95_high", "permutation_p_two_sided"]].to_string(index=False))
    print(json.dumps(run_summary, indent=2))


if __name__ == "__main__":
    main()
