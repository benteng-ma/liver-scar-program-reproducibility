from __future__ import annotations

import hashlib
import json
from pathlib import Path

import numpy as np
import pandas as pd
from scipy.io import mmread
from scipy.stats import rankdata

from audit_gse202379_gates import contrast_groups
from random_module_benchmark_gse202379 import vectorized_hedges_g


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
    detail_rows: list[dict[str, object]] = []

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
        lineage_ranks = ranks[:, lineage_columns]

        lineage_effects = effects[effects["lineage"].eq(lineage)]
        for program_id, program_effects in lineage_effects.groupby("program_id"):
            program_rows = programs[programs["program_id"].eq(program_id)]
            measured_genes = sorted(
                set(program_rows["gene_symbol"].astype(str).str.upper())
                & set(gene_to_index)
            )
            if len(measured_genes) < 2:
                continue
            indices = np.array([gene_to_index[gene] for gene in measured_genes], dtype=int)
            module_size = len(indices)
            rank_values = lineage_ranks[indices, :]
            z_values = standardized[indices, :]
            rank_sum = rank_values.sum(axis=0)
            z_sum = z_values.sum(axis=0)
            loo_rank_means = (rank_sum[None, :] - rank_values) / (module_size - 1)
            theoretical_min = module_size / 2
            theoretical_max = (2 * ranks.shape[0] - module_size + 2) / 2
            loo_sing = (
                2
                * (loo_rank_means - theoretical_min)
                / (theoretical_max - theoretical_min)
                - 1
            )
            loo_zmean = (z_sum[None, :] - z_values) / (module_size - 1)

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
                    [global_to_local[int(index)] for index in selected_global], dtype=int
                )
                selected_group = group.loc[selected_global]
                control_columns = selected_local[
                    selected_group.eq("control").to_numpy()
                ]
                case_columns = selected_local[selected_group.eq("case").to_numpy()]
                loo_scores = (
                    loo_sing if effect.score_method == "singscore" else loo_zmean
                )
                loo_g = vectorized_hedges_g(
                    loo_scores, control_columns, case_columns
                )
                for gene, effect_value in zip(measured_genes, loo_g):
                    detail_rows.append(
                        {
                            "dataset_id": "GSE202379",
                            "contrast": effect.contrast,
                            "lineage": lineage,
                            "program_id": program_id,
                            "score_method": effect.score_method,
                            "analysis_tier": effect.analysis_tier,
                            "cell_gate": effect.cell_gate,
                            "real_hedges_g": effect.hedges_g,
                            "dropped_gene": gene,
                            "leave_one_gene_out_hedges_g": float(effect_value),
                        }
                    )

    details = pd.DataFrame(detail_rows)
    output_dir = repo / "results" / "robustness"
    output_dir.mkdir(parents=True, exist_ok=True)
    detail_path = output_dir / "gse202379_leave_one_gene_out.csv.gz"
    details.to_csv(detail_path, index=False, compression="gzip")
    summary = (
        details.groupby(
            [
                "dataset_id",
                "contrast",
                "lineage",
                "program_id",
                "score_method",
                "analysis_tier",
                "cell_gate",
                "real_hedges_g",
            ],
            as_index=False,
        )
        .agg(
            genes_dropped=("dropped_gene", "count"),
            loo_g_min=("leave_one_gene_out_hedges_g", "min"),
            loo_g_median=("leave_one_gene_out_hedges_g", "median"),
            loo_g_max=("leave_one_gene_out_hedges_g", "max"),
            loo_positive_fraction=(
                "leave_one_gene_out_hedges_g",
                lambda values: float((values > 0).mean()),
            ),
        )
    )
    summary["real_expected_direction"] = summary["real_hedges_g"].gt(0)
    summary["all_loo_expected_direction"] = summary["loo_g_min"].gt(0)
    summary["sign_stable_vs_real"] = np.where(
        summary["real_hedges_g"].gt(0),
        summary["loo_g_min"].gt(0),
        summary["loo_g_max"].lt(0),
    )
    summary_path = output_dir / "gse202379_leave_one_gene_out_summary.csv"
    summary.to_csv(summary_path, index=False)
    run_summary = {
        "dataset_id": "GSE202379",
        "effect_rows": len(summary),
        "leave_one_gene_out_rows": len(details),
        "primary_effect_rows_all_loo_expected_direction": int(
            (
                summary["analysis_tier"].eq("primary")
                & summary["all_loo_expected_direction"]
            ).sum()
        ),
        "detail_sha256": hashlib.sha256(detail_path.read_bytes()).hexdigest().upper(),
        "summary_sha256": hashlib.sha256(summary_path.read_bytes()).hexdigest().upper(),
    }
    (repo / "results" / "logs" / "gse202379_leave_one_gene_out_run.json").write_text(
        json.dumps(run_summary, indent=2) + "\n", encoding="utf-8"
    )
    print(summary[summary["analysis_tier"].eq("primary")].to_string(index=False))
    print(json.dumps(run_summary, indent=2))


if __name__ == "__main__":
    main()
