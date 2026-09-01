from __future__ import annotations

import csv
import gzip
import hashlib
import json
from pathlib import Path

import numpy as np
import pandas as pd
from scipy.io import mmread
from scipy.stats import rankdata

from analyze_gse202379_programs import SEED, stable_seed
from audit_gse202379_gates import contrast_groups
from random_module_benchmark_gse202379 import vectorized_hedges_g


N_RANDOMIZATIONS = 1_000


def balanced_random_signs(
    modules: int, genes: int, rng: np.random.Generator
) -> np.ndarray:
    signs = np.ones((modules, genes), dtype=np.float32)
    negatives = genes // 2
    for row in signs:
        row[rng.choice(genes, size=negatives, replace=False)] = -1
    return signs


def main() -> None:
    repo = Path(__file__).resolve().parents[1]
    interim = repo / "data" / "interim" / "GSE202379"
    genes = pd.read_csv(interim / "genes.csv")["gene"].astype(str)
    gene_to_index = {gene.upper(): index for index, gene in enumerate(genes)}
    manifest = pd.read_csv(interim / "donor_lineage_manifest.csv")
    manifest["eligible_primary_30"] = manifest["n_cells"].ge(30)
    manifest["eligible_sensitivity_20"] = manifest["n_cells"].ge(20)
    counts = mmread(interim / "donor_lineage_raw_counts.mtx").tocsc()
    library_sizes = np.asarray(counts.sum(axis=0)).ravel()
    log_cpm = np.log2(counts.toarray() / library_sizes * 1_000_000 + 1)
    ranks = rankdata(log_cpm, axis=0, method="average")
    centered_ranks = 2 * (ranks - (len(genes) + 1) / 2) / (len(genes) - 1)
    programs = pd.read_csv(repo / "literature" / "program_inventory.csv")
    effects = pd.concat(
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
    groups = contrast_groups(manifest)
    output_dir = repo / "results" / "random_controls"
    output_dir.mkdir(parents=True, exist_ok=True)
    detail_path = output_dir / "gse202379_direction_randomization_effects.csv.gz"
    benchmark_rows: list[dict[str, object]] = []

    with gzip.open(detail_path, "wt", encoding="utf-8", newline="") as detail_file:
        writer = csv.DictWriter(
            detail_file,
            fieldnames=[
                "dataset_id",
                "contrast",
                "lineage",
                "program_id",
                "score_method",
                "analysis_tier",
                "randomization_id",
                "randomized_direction_hedges_g",
            ],
        )
        writer.writeheader()
        for lineage, lineage_manifest in manifest.groupby("harmonized_lineage", sort=False):
            lineage_columns = lineage_manifest.index.to_numpy()
            global_to_local = {
                int(global_index): local_index
                for local_index, global_index in enumerate(lineage_columns)
            }
            reference_columns = lineage_manifest.index[
                lineage_manifest["eligible_sensitivity_20"]
            ].to_numpy()
            means = log_cpm[:, reference_columns].mean(axis=1)
            standard_deviations = log_cpm[:, reference_columns].std(axis=1, ddof=1)
            invariant = standard_deviations == 0
            safe_sd = standard_deviations.copy()
            safe_sd[invariant] = 1
            standardized = (
                log_cpm[:, lineage_columns] - means[:, None]
            ) / safe_sd[:, None]
            standardized[invariant, :] = 0
            lineage_centered_ranks = centered_ranks[:, lineage_columns]
            for program_id, program_effects in effects[
                effects["lineage"].eq(lineage)
            ].groupby("program_id"):
                rows = programs[programs["program_id"].eq(program_id)]
                measured_genes = sorted(
                    set(rows["gene_symbol"].astype(str).str.upper())
                    & set(gene_to_index)
                )
                indices = np.array(
                    [gene_to_index[gene] for gene in measured_genes], dtype=int
                )
                rng = np.random.default_rng(
                    stable_seed("GSE202379", lineage, program_id, "direction_randomization")
                )
                signs = balanced_random_signs(
                    N_RANDOMIZATIONS, len(indices), rng
                ).astype(float)
                randomized_rank_scores = (
                    signs @ lineage_centered_ranks[indices, :] / len(indices)
                )
                randomized_z_scores = signs @ standardized[indices, :] / len(indices)
                for effect in program_effects.itertuples(index=False):
                    eligibility = (
                        "eligible_primary_30"
                        if effect.cell_gate == 30
                        else "eligible_sensitivity_20"
                    )
                    group = groups[effect.contrast]
                    selected = (
                        manifest["harmonized_lineage"].eq(lineage)
                        & manifest[eligibility]
                        & group.ne("excluded")
                    )
                    selected_global = manifest.index[selected].to_numpy()
                    selected_local = np.array(
                        [global_to_local[int(index)] for index in selected_global],
                        dtype=int,
                    )
                    selected_group = group.loc[selected_global]
                    control_columns = selected_local[
                        selected_group.eq("control").to_numpy()
                    ]
                    case_columns = selected_local[selected_group.eq("case").to_numpy()]
                    randomized_scores = (
                        randomized_rank_scores
                        if effect.score_method == "singscore"
                        else randomized_z_scores
                    )
                    randomized_g = vectorized_hedges_g(
                        randomized_scores, control_columns, case_columns
                    )
                    valid = randomized_g[np.isfinite(randomized_g)]
                    if len(valid) != N_RANDOMIZATIONS:
                        raise RuntimeError("undefined direction-randomized Hedges g")
                    real_g = float(effect.hedges_g)
                    q95 = float(np.quantile(valid, 0.95))
                    benchmark_rows.append(
                        {
                            "dataset_id": "GSE202379",
                            "contrast": effect.contrast,
                            "lineage": lineage,
                            "program_id": program_id,
                            "score_method": effect.score_method,
                            "analysis_tier": effect.analysis_tier,
                            "cell_gate": effect.cell_gate,
                            "real_hedges_g": real_g,
                            "randomized_g_median": float(np.median(valid)),
                            "randomized_g_95th_percentile": q95,
                            "real_effect_percentile": float((valid <= real_g).mean()),
                            "empirical_p_one_sided": (
                                1 + int((valid >= real_g).sum())
                            )
                            / (N_RANDOMIZATIONS + 1),
                            "above_direction_randomized_95th": real_g > q95,
                            "randomizations": N_RANDOMIZATIONS,
                            "negative_directions_per_randomization": len(indices) // 2,
                        }
                    )
                    for randomization_id, effect_value in enumerate(
                        randomized_g, start=1
                    ):
                        writer.writerow(
                            {
                                "dataset_id": "GSE202379",
                                "contrast": effect.contrast,
                                "lineage": lineage,
                                "program_id": program_id,
                                "score_method": effect.score_method,
                                "analysis_tier": effect.analysis_tier,
                                "randomization_id": f"D{randomization_id:04d}",
                                "randomized_direction_hedges_g": float(effect_value),
                            }
                        )

    benchmark = pd.DataFrame(benchmark_rows).sort_values(
        ["analysis_tier", "contrast", "lineage", "program_id", "score_method"]
    )
    benchmark_path = output_dir / "gse202379_direction_randomization_benchmark.csv"
    benchmark.to_csv(benchmark_path, index=False)
    summary = {
        "dataset_id": "GSE202379",
        "seed": SEED,
        "randomizations_per_program": N_RANDOMIZATIONS,
        "balanced_sign_assignment": True,
        "effect_rows": len(benchmark),
        "primary_effect_rows_above_randomized_95th": int(
            (
                benchmark["analysis_tier"].eq("primary")
                & benchmark["above_direction_randomized_95th"]
            ).sum()
        ),
        "detail_sha256": hashlib.sha256(detail_path.read_bytes()).hexdigest().upper(),
        "benchmark_sha256": hashlib.sha256(benchmark_path.read_bytes()).hexdigest().upper(),
    }
    (repo / "results" / "logs" / "gse202379_direction_randomization_run.json").write_text(
        json.dumps(summary, indent=2) + "\n", encoding="utf-8"
    )
    print(benchmark[benchmark["analysis_tier"].eq("primary")].to_string(index=False))
    print(json.dumps(summary, indent=2))


if __name__ == "__main__":
    main()
