from __future__ import annotations

import argparse
import csv
import gzip
import hashlib
import json
from pathlib import Path

import numpy as np
import pandas as pd
from scipy.io import mmread
from scipy.stats import rankdata

from analyze_gse202379_programs import stable_seed
from audit_gse181483_gates import contrast_groups as gse181483_contrast_groups
from audit_gse210077_watson6_gates import contrast_groups as watson_contrast_groups
from audit_gse290642_gates import contrast_groups as gse290642_contrast_groups
from audit_gse256398_gates import CONTRASTS as gse256398_contrasts
from random_module_benchmark_gse202379 import (
    N_RANDOM_MODULES,
    deciles,
    matched_modules,
    score_modules,
    vectorized_hedges_g,
)


def configuration(dataset: str) -> dict[str, object]:
    if dataset == "GSE181483_human":
        return {
            "interim": "GSE181483",
            "prefix": "gse181483",
            "effects": "gse181483_directional_effects.csv",
        }
    if dataset == "GSE244832":
        return {
            "interim": "GSE244832",
            "prefix": "gse244832",
            "effects": "gse244832_sensitivity_effects.csv",
        }
    if dataset == "GSE290642_human":
        return {
            "interim": "GSE290642",
            "prefix": "gse290642",
            "effects": "gse290642_sensitivity_effects.csv",
        }
    if dataset == "GSE210077_Watson6":
        return {
            "interim": "GSE210077_Watson6",
            "prefix": "gse210077_watson6",
            "effects": "gse210077_watson6_sensitivity_effects.csv",
        }
    if dataset == "GSE256398_human":
        return {
            "interim": "GSE256398",
            "prefix": "gse256398",
            "effects": "gse256398_program_effects.csv",
            "effects_directory": "phase3",
        }
    raise ValueError(f"unsupported validation dataset: {dataset}")


def groups_for(dataset: str, manifest: pd.DataFrame) -> dict[str, pd.Series]:
    if dataset == "GSE181483_human":
        return gse181483_contrast_groups(manifest)
    if dataset == "GSE244832":
        groups = manifest["disease_group"].map({"normal": "control", "MASH": "case"}).fillna("excluded")
        return {"mash_f2f4_group_vs_normal_sensitivity": groups}
    if dataset == "GSE210077_Watson6":
        return watson_contrast_groups(manifest)
    if dataset == "GSE256398_human":
        result: dict[str, pd.Series] = {}
        for contrast, (controls, cases) in gse256398_contrasts.items():
            result[contrast] = manifest["disease_group"].map(
                {
                    **{value: "control" for value in controls},
                    **{value: "case" for value in cases},
                }
            ).fillna("excluded")
        return result
    return gse290642_contrast_groups(manifest)


def main() -> None:
    parser = argparse.ArgumentParser()
    parser.add_argument(
        "--dataset",
        choices=["GSE181483_human", "GSE244832", "GSE290642_human", "GSE210077_Watson6", "GSE256398_human"],
        required=True,
    )
    args = parser.parse_args()
    dataset = args.dataset
    cfg = configuration(dataset)
    repo = Path(__file__).resolve().parents[1]
    interim = repo / "data" / "interim" / str(cfg["interim"])
    prefix = str(cfg["prefix"])
    random_interim = interim / "random_modules"
    output_dir = repo / "results" / "random_controls"
    random_interim.mkdir(parents=True, exist_ok=True)
    output_dir.mkdir(parents=True, exist_ok=True)

    union_genes = pd.read_csv(interim / "genes.csv")["gene"].astype(str)
    counts = mmread(interim / "donor_lineage_raw_counts.mtx").tocsc()
    if dataset == "GSE290642_human":
        shared = set(pd.read_csv(interim / "shared_genes.csv")["gene"].astype(str))
        shared_rows = np.array([index for index, gene in enumerate(union_genes) if gene in shared], dtype=int)
        genes = union_genes.iloc[shared_rows].reset_index(drop=True)
        counts = counts[shared_rows, :]
    else:
        genes = union_genes
    gene_to_index = {gene.upper(): index for index, gene in enumerate(genes)}
    manifest = pd.read_csv(interim / "donor_lineage_manifest.csv")
    manifest["eligible_30"] = manifest["n_cells"].ge(30)
    manifest["eligible_20"] = manifest["n_cells"].ge(20)
    library_sizes = np.asarray(counts.sum(axis=0)).ravel()
    log_cpm = np.log2(counts.toarray() / library_sizes * 1_000_000 + 1)
    ranks = rankdata(log_cpm, axis=0, method="average")
    programs = pd.read_csv(repo / "literature" / "program_inventory.csv")
    coverage = pd.read_csv(repo / "results" / "qc" / f"{prefix}_program_coverage.csv")
    coverage_lookup = coverage.set_index("program_id").to_dict("index")
    effects_directory = str(cfg.get("effects_directory", "sensitivity"))
    effects = pd.read_csv(repo / "results" / effects_directory / str(cfg["effects"]))
    groups_by_contrast = groups_for(dataset, manifest)

    membership_path = random_interim / "matched_random_module_membership.csv.gz"
    random_effect_path = output_dir / f"{prefix}_random_module_effects.csv.gz"
    benchmark_rows: list[dict[str, object]] = []
    matching_rows: list[dict[str, object]] = []
    with gzip.open(membership_path, "wt", encoding="utf-8", newline="") as member_file, gzip.open(random_effect_path, "wt", encoding="utf-8", newline="") as effect_file:
        member_writer = csv.DictWriter(member_file, fieldnames=["dataset_id", "lineage", "program_id", "module_id", "target_gene", "random_gene", "expression_detection_bin_distance"])
        effect_writer = csv.DictWriter(effect_file, fieldnames=["dataset_id", "contrast", "lineage", "program_id", "score_method", "analysis_tier", "module_id", "random_hedges_g"])
        member_writer.writeheader()
        effect_writer.writeheader()

        for lineage, lineage_manifest in manifest.groupby("harmonized_lineage", sort=False):
            lineage_columns = lineage_manifest.index.to_numpy()
            global_to_local = {int(global_index): local_index for local_index, global_index in enumerate(lineage_columns)}
            reference_columns = lineage_manifest.index[lineage_manifest["eligible_20"]].to_numpy()
            mean_expression = log_cpm[:, reference_columns].mean(axis=1)
            detection_rate = np.asarray((counts[:, reference_columns] > 0).mean(axis=1)).ravel()
            expression_decile = deciles(mean_expression)
            detection_decile = deciles(detection_rate)
            standard_deviations = log_cpm[:, reference_columns].std(axis=1, ddof=1)
            invariant = standard_deviations == 0
            safe_sd = standard_deviations.copy()
            safe_sd[invariant] = 1
            standardized = (log_cpm[:, lineage_columns] - mean_expression[:, None]) / safe_sd[:, None]
            standardized[invariant, :] = 0
            lineage_ranks = ranks[:, lineage_columns]

            lineage_programs = programs[programs["cell_lineage"].eq(lineage)]
            for program_id, rows in lineage_programs.groupby("program_id"):
                program_effects = effects[effects["program_id"].eq(program_id)]
                if program_effects.empty or coverage_lookup[program_id]["coverage_tier"] == "not_evaluated":
                    continue
                measured_genes = sorted(set(rows["gene_symbol"].astype(str).str.upper()) & set(gene_to_index))
                target_indices = np.array([gene_to_index[gene] for gene in measured_genes], dtype=int)
                rng = np.random.default_rng(stable_seed(dataset, lineage, program_id, "random_modules"))
                modules, distances, target_positions = matched_modules(target_indices, expression_decile, detection_decile, rng)
                random_sing, random_zmean = score_modules(modules, lineage_ranks, standardized)
                for module_index in range(N_RANDOM_MODULES):
                    module_id = f"R{module_index + 1:04d}"
                    for output_position in range(modules.shape[1]):
                        target_position = target_positions[module_index, output_position]
                        member_writer.writerow(
                            {
                                "dataset_id": dataset,
                                "lineage": lineage,
                                "program_id": program_id,
                                "module_id": module_id,
                                "target_gene": genes.iloc[target_indices[target_position]],
                                "random_gene": genes.iloc[modules[module_index, output_position]],
                                "expression_detection_bin_distance": int(distances[module_index, output_position]),
                            }
                        )
                matching_rows.append(
                    {
                        "dataset_id": dataset,
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
                for effect in program_effects.itertuples(index=False):
                    group = groups_by_contrast[effect.contrast]
                    eligibility_column = "eligible_30" if effect.cell_gate == 30 else "eligible_20"
                    selected = manifest["harmonized_lineage"].eq(lineage) & manifest[eligibility_column] & group.ne("excluded")
                    selected_global = manifest.index[selected].to_numpy()
                    selected_local = np.array([global_to_local[int(index)] for index in selected_global], dtype=int)
                    selected_groups = group.loc[selected_global]
                    control_columns = selected_local[selected_groups.eq("control").to_numpy()]
                    case_columns = selected_local[selected_groups.eq("case").to_numpy()]
                    module_scores = random_sing if effect.score_method == "singscore" else random_zmean
                    random_g = vectorized_hedges_g(module_scores, control_columns, case_columns)
                    valid = random_g[np.isfinite(random_g)]
                    if len(valid) != N_RANDOM_MODULES:
                        raise RuntimeError("a random module produced undefined Hedges g")
                    real_g = float(effect.hedges_g)
                    q95 = float(np.quantile(valid, 0.95))
                    benchmark_rows.append(
                        {
                            "dataset_id": dataset,
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
                            "real_effect_percentile": float((valid <= real_g).mean()),
                            "empirical_p_one_sided": (1 + int((valid >= real_g).sum())) / (N_RANDOM_MODULES + 1),
                            "above_random_95th_percentile": real_g > q95,
                            "random_modules": N_RANDOM_MODULES,
                        }
                    )
                    for module_index, effect_value in enumerate(random_g, start=1):
                        effect_writer.writerow(
                            {
                                "dataset_id": dataset,
                                "contrast": effect.contrast,
                                "lineage": lineage,
                                "program_id": program_id,
                                "score_method": effect.score_method,
                                "analysis_tier": effect.analysis_tier,
                                "module_id": f"R{module_index:04d}",
                                "random_hedges_g": float(effect_value),
                            }
                        )

    benchmark = pd.DataFrame(benchmark_rows).sort_values(["contrast", "lineage", "program_id", "score_method"])
    matching = pd.DataFrame(matching_rows).sort_values(["lineage", "program_id"])
    benchmark_path = output_dir / f"{prefix}_random_module_benchmark.csv"
    matching_path = output_dir / f"{prefix}_random_module_matching_qc.csv"
    benchmark.to_csv(benchmark_path, index=False)
    matching.to_csv(matching_path, index=False)
    summary = {
        "dataset_id": dataset,
        "seed": 20260830,
        "random_modules_per_program": N_RANDOM_MODULES,
        "programs_benchmarked": int(matching["program_id"].nunique()),
        "effect_rows_benchmarked": len(benchmark),
        "effect_rows_above_random_95th": int(benchmark["above_random_95th_percentile"].sum()),
        "membership_sha256": hashlib.sha256(membership_path.read_bytes()).hexdigest().upper(),
        "random_effects_sha256": hashlib.sha256(random_effect_path.read_bytes()).hexdigest().upper(),
        "benchmark_sha256": hashlib.sha256(benchmark_path.read_bytes()).hexdigest().upper(),
    }
    (repo / "results" / "logs" / f"{prefix}_random_module_run.json").write_text(json.dumps(summary, indent=2) + "\n", encoding="utf-8")
    print(benchmark.to_string(index=False))
    print(matching.to_string(index=False))
    print(json.dumps(summary, indent=2))


if __name__ == "__main__":
    main()
