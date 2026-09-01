from __future__ import annotations

import csv
import gzip
import hashlib
import json
import math
from collections import Counter
from pathlib import Path

import numpy as np
import pandas as pd
from scipy.io import mmread
from scipy.stats import rankdata

from analyze_gse202379_programs import SEED, singscore_up, stable_seed
from audit_gse202379_gates import contrast_groups


N_RANDOM_MODULES = 1_000


def deciles(values: np.ndarray) -> np.ndarray:
    ranks = rankdata(values, method="average")
    return np.minimum(((ranks - 1) * 10 // len(values)).astype(int), 9)


def matched_modules(
    target_indices: np.ndarray,
    expression_decile: np.ndarray,
    detection_decile: np.ndarray,
    rng: np.random.Generator,
) -> tuple[np.ndarray, np.ndarray, np.ndarray]:
    target_set = set(map(int, target_indices))
    all_bins = [
        (expression, detection) for expression in range(10) for detection in range(10)
    ]
    base_pools: dict[tuple[int, int], np.ndarray] = {}
    for key in all_bins:
        candidates = np.where(
            (expression_decile == key[0]) & (detection_decile == key[1])
        )[0]
        base_pools[key] = np.array(
            [index for index in candidates if int(index) not in target_set], dtype=int
        )
    target_bins = [
        (int(expression_decile[index]), int(detection_decile[index]))
        for index in target_indices
    ]
    target_bin_counts = Counter(target_bins)
    ordered_targets = sorted(
        range(len(target_indices)),
        key=lambda position: (
            len(base_pools[target_bins[position]])
            / max(target_bin_counts[target_bins[position]], 1),
            target_bins[position],
            int(target_indices[position]),
        ),
    )
    modules = np.empty((N_RANDOM_MODULES, len(target_indices)), dtype=np.int32)
    distances = np.empty_like(modules, dtype=np.int8)
    target_positions = np.empty_like(modules, dtype=np.int32)
    for module_index in range(N_RANDOM_MODULES):
        pools = {
            key: list(rng.permutation(values)) for key, values in base_pools.items()
        }
        for output_position, target_position in enumerate(ordered_targets):
            target_bin = target_bins[target_position]
            bin_order = sorted(
                all_bins,
                key=lambda key: (
                    abs(key[0] - target_bin[0]) + abs(key[1] - target_bin[1]),
                    rng.random(),
                ),
            )
            chosen_bin = next((key for key in bin_order if pools[key]), None)
            if chosen_bin is None:
                raise RuntimeError("matched background exhausted before module completion")
            modules[module_index, output_position] = int(pools[chosen_bin].pop())
            distances[module_index, output_position] = (
                abs(chosen_bin[0] - target_bin[0])
                + abs(chosen_bin[1] - target_bin[1])
            )
            target_positions[module_index, output_position] = target_position
        if len(np.unique(modules[module_index])) != len(target_indices):
            raise RuntimeError("random module contains a duplicated gene")
    return modules, distances, target_positions


def score_modules(
    modules: np.ndarray,
    ranks: np.ndarray,
    standardized: np.ndarray,
) -> tuple[np.ndarray, np.ndarray]:
    n_modules = modules.shape[0]
    n_donors = ranks.shape[1]
    sing = np.empty((n_modules, n_donors), dtype=float)
    zmean = np.empty((n_modules, n_donors), dtype=float)
    for start in range(0, n_modules, 50):
        stop = min(start + 50, n_modules)
        batch = modules[start:stop]
        mean_ranks = ranks[batch, :].mean(axis=1)
        module_size = modules.shape[1]
        theoretical_min = (module_size + 1) / 2
        theoretical_max = (2 * ranks.shape[0] - module_size + 1) / 2
        sing[start:stop] = (
            2
            * (mean_ranks - theoretical_min)
            / (theoretical_max - theoretical_min)
            - 1
        )
        zmean[start:stop] = standardized[batch, :].mean(axis=1)
    return sing, zmean


def vectorized_hedges_g(
    module_scores: np.ndarray,
    control_columns: np.ndarray,
    case_columns: np.ndarray,
) -> np.ndarray:
    control = module_scores[:, control_columns]
    case = module_scores[:, case_columns]
    n_control = control.shape[1]
    n_case = case.shape[1]
    degrees_freedom = n_control + n_case - 2
    pooled_variance = (
        (n_control - 1) * control.var(axis=1, ddof=1)
        + (n_case - 1) * case.var(axis=1, ddof=1)
    ) / degrees_freedom
    correction = 1 - 3 / (4 * degrees_freedom - 1)
    difference = case.mean(axis=1) - control.mean(axis=1)
    with np.errstate(divide="ignore", invalid="ignore"):
        effects = correction * difference / np.sqrt(pooled_variance)
    effects[~np.isfinite(effects)] = np.nan
    return effects


def main() -> None:
    repo = Path(__file__).resolve().parents[1]
    interim = repo / "data" / "interim" / "GSE202379"
    random_interim = interim / "random_modules"
    random_interim.mkdir(parents=True, exist_ok=True)
    output_dir = repo / "results" / "random_controls"
    output_dir.mkdir(parents=True, exist_ok=True)

    genes = pd.read_csv(interim / "genes.csv")["gene"].astype(str)
    gene_to_index = {gene.upper(): index for index, gene in enumerate(genes)}
    manifest = pd.read_csv(interim / "donor_lineage_manifest.csv")
    manifest["eligible_primary_30"] = manifest["n_cells"].ge(30)
    manifest["eligible_sensitivity_20"] = manifest["n_cells"].ge(20)
    counts = mmread(interim / "donor_lineage_raw_counts.mtx").tocsc()
    library_sizes = np.asarray(counts.sum(axis=0)).ravel()
    log_cpm = np.log2(counts.toarray() / library_sizes * 1_000_000 + 1)
    ranks = rankdata(log_cpm, axis=0, method="average")
    programs = pd.read_csv(repo / "literature" / "program_inventory.csv")
    coverage = pd.read_csv(repo / "results" / "qc" / "gse202379_program_coverage.csv")
    coverage_lookup = coverage.set_index("program_id").to_dict("index")
    primary_effects = pd.read_csv(
        repo / "results" / "primary" / "gse202379_primary_effects.csv"
    )
    sensitivity_effects = pd.read_csv(
        repo / "results" / "sensitivity" / "gse202379_sensitivity_effects.csv"
    )
    real_effects = pd.concat([primary_effects, sensitivity_effects], ignore_index=True)
    groups_by_contrast = contrast_groups(manifest)

    membership_path = random_interim / "matched_random_module_membership.csv.gz"
    random_effect_path = output_dir / "gse202379_random_module_effects.csv.gz"
    benchmark_rows: list[dict[str, object]] = []
    matching_rows: list[dict[str, object]] = []
    with gzip.open(membership_path, "wt", encoding="utf-8", newline="") as member_file, gzip.open(
        random_effect_path, "wt", encoding="utf-8", newline=""
    ) as effect_file:
        member_writer = csv.DictWriter(
            member_file,
            fieldnames=[
                "dataset_id",
                "lineage",
                "program_id",
                "module_id",
                "target_gene",
                "random_gene",
                "expression_detection_bin_distance",
            ],
        )
        effect_writer = csv.DictWriter(
            effect_file,
            fieldnames=[
                "dataset_id",
                "contrast",
                "lineage",
                "program_id",
                "score_method",
                "analysis_tier",
                "module_id",
                "random_hedges_g",
            ],
        )
        member_writer.writeheader()
        effect_writer.writeheader()

        for lineage, lineage_manifest in manifest.groupby("harmonized_lineage", sort=False):
            lineage_columns = lineage_manifest.index.to_numpy()
            global_to_local = {
                int(global_index): local_index
                for local_index, global_index in enumerate(lineage_columns)
            }
            reference_columns = lineage_manifest.index[
                lineage_manifest["eligible_sensitivity_20"]
            ].to_numpy()
            mean_expression = log_cpm[:, reference_columns].mean(axis=1)
            detection_rate = np.asarray(
                (counts[:, reference_columns] > 0).mean(axis=1)
            ).ravel()
            expression_decile = deciles(mean_expression)
            detection_decile = deciles(detection_rate)
            standard_deviations = log_cpm[:, reference_columns].std(axis=1, ddof=1)
            invariant = standard_deviations == 0
            safe_sd = standard_deviations.copy()
            safe_sd[invariant] = 1
            standardized = (
                log_cpm[:, lineage_columns] - mean_expression[:, None]
            ) / safe_sd[:, None]
            standardized[invariant, :] = 0
            lineage_ranks = ranks[:, lineage_columns]

            lineage_programs = programs[programs["cell_lineage"].eq(lineage)]
            for program_id, rows in lineage_programs.groupby("program_id"):
                if coverage_lookup[program_id]["coverage_tier"] == "not_evaluated":
                    continue
                measured_genes = sorted(
                    set(rows["gene_symbol"].astype(str).str.upper())
                    & set(gene_to_index)
                )
                target_indices = np.array(
                    [gene_to_index[gene] for gene in measured_genes], dtype=int
                )
                rng = np.random.default_rng(
                    stable_seed("GSE202379", lineage, program_id, "random_modules")
                )
                modules, distances, target_positions = matched_modules(
                    target_indices,
                    expression_decile,
                    detection_decile,
                    rng,
                )
                random_sing, random_zmean = score_modules(
                    modules, lineage_ranks, standardized
                )
                for module_index in range(N_RANDOM_MODULES):
                    module_id = f"R{module_index + 1:04d}"
                    for output_position in range(modules.shape[1]):
                        target_position = target_positions[module_index, output_position]
                        member_writer.writerow(
                            {
                                "dataset_id": "GSE202379",
                                "lineage": lineage,
                                "program_id": program_id,
                                "module_id": module_id,
                                "target_gene": genes.iloc[target_indices[target_position]],
                                "random_gene": genes.iloc[
                                    modules[module_index, output_position]
                                ],
                                "expression_detection_bin_distance": int(
                                    distances[module_index, output_position]
                                ),
                            }
                        )
                matching_rows.append(
                    {
                        "dataset_id": "GSE202379",
                        "lineage": lineage,
                        "program_id": program_id,
                        "measured_module_size": len(target_indices),
                        "random_modules": N_RANDOM_MODULES,
                        "exact_bin_match_fraction": float((distances == 0).mean()),
                        "expanded_bin_match_fraction": float((distances > 0).mean()),
                        "maximum_bin_distance": int(distances.max()),
                        "within_module_replacement": False,
                    }
                )

                program_effects = real_effects[real_effects["program_id"].eq(program_id)]
                for effect in program_effects.itertuples(index=False):
                    group = groups_by_contrast[effect.contrast]
                    eligibility_column = (
                        "eligible_primary_30" if effect.cell_gate == 30 else "eligible_sensitivity_20"
                    )
                    selected = (
                        manifest["harmonized_lineage"].eq(lineage)
                        & manifest[eligibility_column]
                        & group.ne("excluded")
                    )
                    selected_global = manifest.index[selected].to_numpy()
                    selected_local = np.array(
                        [global_to_local[int(index)] for index in selected_global], dtype=int
                    )
                    selected_groups = group.loc[selected_global]
                    control_columns = selected_local[
                        selected_groups.eq("control").to_numpy()
                    ]
                    case_columns = selected_local[selected_groups.eq("case").to_numpy()]
                    module_scores = (
                        random_sing
                        if effect.score_method == "singscore"
                        else random_zmean
                    )
                    random_g = vectorized_hedges_g(
                        module_scores, control_columns, case_columns
                    )
                    valid = random_g[np.isfinite(random_g)]
                    if len(valid) != N_RANDOM_MODULES:
                        raise RuntimeError("a random module produced undefined Hedges g")
                    q95 = float(np.quantile(valid, 0.95))
                    real_g = float(effect.hedges_g)
                    empirical_p = (1 + int((valid >= real_g).sum())) / (
                        N_RANDOM_MODULES + 1
                    )
                    percentile = float((valid <= real_g).mean())
                    benchmark_rows.append(
                        {
                            "dataset_id": "GSE202379",
                            "contrast": effect.contrast,
                            "lineage": lineage,
                            "program_id": program_id,
                            "score_method": effect.score_method,
                            "analysis_tier": effect.analysis_tier,
                            "cell_gate": effect.cell_gate,
                            "coverage_tier": effect.coverage_tier,
                            "real_hedges_g": real_g,
                            "random_g_median": float(np.median(valid)),
                            "random_g_95th_percentile": q95,
                            "real_effect_percentile": percentile,
                            "empirical_p_one_sided": empirical_p,
                            "above_random_95th_percentile": real_g > q95,
                            "random_modules": N_RANDOM_MODULES,
                        }
                    )
                    for module_index, effect_value in enumerate(random_g, start=1):
                        effect_writer.writerow(
                            {
                                "dataset_id": "GSE202379",
                                "contrast": effect.contrast,
                                "lineage": lineage,
                                "program_id": program_id,
                                "score_method": effect.score_method,
                                "analysis_tier": effect.analysis_tier,
                                "module_id": f"R{module_index:04d}",
                                "random_hedges_g": float(effect_value),
                            }
                        )

    benchmark = pd.DataFrame(benchmark_rows).sort_values(
        ["analysis_tier", "contrast", "lineage", "program_id", "score_method"]
    )
    matching = pd.DataFrame(matching_rows).sort_values(["lineage", "program_id"])
    benchmark_path = output_dir / "gse202379_random_module_benchmark.csv"
    matching_path = output_dir / "gse202379_random_module_matching_qc.csv"
    benchmark.to_csv(benchmark_path, index=False)
    matching.to_csv(matching_path, index=False)
    summary = {
        "dataset_id": "GSE202379",
        "seed": SEED,
        "random_modules_per_program": N_RANDOM_MODULES,
        "programs_benchmarked": matching["program_id"].nunique(),
        "effect_rows_benchmarked": len(benchmark),
        "primary_effect_rows_above_random_95th": int(
            (
                benchmark["analysis_tier"].eq("primary")
                & benchmark["above_random_95th_percentile"]
            ).sum()
        ),
        "membership_sha256": hashlib.sha256(membership_path.read_bytes()).hexdigest().upper(),
        "random_effects_sha256": hashlib.sha256(random_effect_path.read_bytes()).hexdigest().upper(),
        "benchmark_sha256": hashlib.sha256(benchmark_path.read_bytes()).hexdigest().upper(),
    }
    (repo / "results" / "logs" / "gse202379_random_module_run.json").write_text(
        json.dumps(summary, indent=2) + "\n", encoding="utf-8"
    )
    print(benchmark[benchmark["analysis_tier"].eq("primary")].to_string(index=False))
    print(matching.to_string(index=False))
    print(json.dumps(summary, indent=2))


if __name__ == "__main__":
    main()
