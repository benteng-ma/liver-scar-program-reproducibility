from __future__ import annotations

import json
from pathlib import Path

import numpy as np
import pandas as pd
from scipy.io import mmread
from scipy.stats import rankdata

from analyze_gse202379_programs import effect_statistics, stable_seed
from random_module_benchmark_gse202379 import score_modules, vectorized_hedges_g


CONTRASTS = {
    "alcohol_hepatitis_vs_healthy": ({"healthy"}, {"alcohol_hepatitis"}),
    "alcohol_cirrhosis_vs_alcohol_hepatitis": (
        {"alcohol_hepatitis"},
        {"alcohol_cirrhosis"},
    ),
}
GATES = ((30, "primary"), (20, "sensitivity"))
PRIMARY_PROGRAMS = {"RAM2019_ENDO_2", "RAM2019_ENDO_6_SAENDO1"}


def main() -> None:
    repo = Path(__file__).resolve().parents[1]
    interim = repo / "data" / "interim" / "GSE256398"
    output = repo / "results" / "phase5"
    output.mkdir(parents=True, exist_ok=True)

    genes = pd.read_csv(interim / "genes.csv")["gene"].astype(str)
    gene_to_index = {gene.upper(): index for index, gene in enumerate(genes)}
    manifest = pd.read_csv(interim / "donor_lineage_manifest.csv")
    counts = mmread(interim / "donor_lineage_raw_counts.mtx").tocsc()
    library = np.asarray(counts.sum(axis=0)).ravel()
    log_cpm = np.log2(counts.toarray() / library * 1_000_000 + 1)
    ranks = rankdata(log_cpm, axis=0, method="average")
    scores = pd.read_csv(repo / "results" / "phase3" / "gse256398_donor_program_scores.csv.gz")
    programs = pd.read_csv(repo / "literature" / "program_inventory.csv")
    coverage = pd.read_csv(repo / "results" / "qc" / "gse256398_program_coverage.csv")
    coverage_lookup = coverage.set_index("program_id").to_dict("index")
    membership = pd.read_csv(interim / "random_modules" / "matched_random_module_membership.csv.gz")

    effect_rows: list[dict[str, object]] = []
    benchmark_rows: list[dict[str, object]] = []

    for lineage, lineage_manifest in manifest.groupby("harmonized_lineage", sort=False):
        lineage_columns = lineage_manifest.index.to_numpy()
        global_to_local = {int(value): index for index, value in enumerate(lineage_columns)}
        reference_columns = lineage_manifest.index[lineage_manifest["n_cells"].ge(20)].to_numpy()
        mean_expression = log_cpm[:, reference_columns].mean(axis=1)
        standard_deviation = log_cpm[:, reference_columns].std(axis=1, ddof=1)
        invariant = standard_deviation == 0
        safe_sd = standard_deviation.copy()
        safe_sd[invariant] = 1
        standardized = (log_cpm[:, lineage_columns] - mean_expression[:, None]) / safe_sd[:, None]
        standardized[invariant, :] = 0
        lineage_ranks = ranks[:, lineage_columns]

        for program_id, _ in programs[programs["cell_lineage"].eq(lineage)].groupby("program_id"):
            if coverage_lookup[program_id]["coverage_tier"] == "not_evaluated":
                continue
            module_members = membership[
                membership["lineage"].eq(lineage) & membership["program_id"].eq(program_id)
            ]
            module_ids = sorted(module_members["module_id"].unique())
            modules = []
            for module_id in module_ids:
                random_genes = module_members.loc[
                    module_members["module_id"].eq(module_id), "random_gene"
                ].astype(str)
                indices = [gene_to_index[value.upper()] for value in random_genes]
                modules.append(indices)
            lengths = {len(value) for value in modules}
            if len(lengths) != 1 or len(modules) != 1000:
                raise RuntimeError(f"invalid frozen random membership for {program_id}")
            module_array = np.asarray(modules, dtype=int)
            random_sing, random_zmean = score_modules(module_array, lineage_ranks, standardized)

            program_scores = scores[
                scores["program_id"].eq(program_id)
                & scores["lineage"].eq(lineage)
            ].copy()
            for contrast, (controls, cases) in CONTRASTS.items():
                label = manifest["disease_group"].map(
                    {
                        **{value: "control" for value in controls},
                        **{value: "case" for value in cases},
                    }
                ).fillna("excluded")
                for gate, tier in GATES:
                    selected = (
                        manifest["harmonized_lineage"].eq(lineage)
                        & manifest["n_cells"].ge(gate)
                        & label.ne("excluded")
                    )
                    selected_global = manifest.index[selected].to_numpy()
                    selected_label = label.loc[selected_global]
                    n_control = int(selected_label.eq("control").sum())
                    n_case = int(selected_label.eq("case").sum())
                    if n_control < 3 or n_case < 3:
                        continue
                    selected_local = np.array(
                        [global_to_local[int(value)] for value in selected_global], dtype=int
                    )
                    control_local = selected_local[selected_label.eq("control").to_numpy()]
                    case_local = selected_local[selected_label.eq("case").to_numpy()]
                    selected_groups = manifest.loc[selected_global, "group_id"].astype(str)
                    group_label = dict(zip(selected_groups, selected_label))
                    for method in ("singscore", "standardized_mean"):
                        values = program_scores[program_scores["score_method"].eq(method)].copy()
                        values = values[values["group_id"].isin(selected_groups)]
                        values["comparison_group"] = values["group_id"].map(group_label)
                        control_values = values.loc[
                            values["comparison_group"].eq("control"), "score"
                        ].to_numpy(float)
                        case_values = values.loc[
                            values["comparison_group"].eq("case"), "score"
                        ].to_numpy(float)
                        stats = effect_statistics(
                            control_values,
                            case_values,
                            stable_seed("phase5", contrast, lineage, program_id, method, str(gate)),
                        )
                        effect_rows.append(
                            {
                                "dataset_id": "GSE256398_human",
                                "contrast": contrast,
                                "lineage": lineage,
                                "program_id": program_id,
                                "score_method": method,
                                "cell_gate": gate,
                                "analysis_tier": tier,
                                "coverage_tier": coverage_lookup[program_id]["coverage_tier"],
                                **stats,
                            }
                        )
                        module_scores = random_sing if method == "singscore" else random_zmean
                        random_g = vectorized_hedges_g(
                            module_scores, control_local, case_local
                        )
                        valid = random_g[np.isfinite(random_g)]
                        if len(valid) != 1000:
                            raise RuntimeError("undefined frozen random-module effect")
                        real_g = float(stats["hedges_g"])
                        q95 = float(np.quantile(valid, 0.95))
                        benchmark_rows.append(
                            {
                                "dataset_id": "GSE256398_human",
                                "contrast": contrast,
                                "lineage": lineage,
                                "program_id": program_id,
                                "score_method": method,
                                "cell_gate": gate,
                                "analysis_tier": tier,
                                "real_hedges_g": real_g,
                                "random_g_median": float(np.median(valid)),
                                "random_g_95th_percentile": q95,
                                "real_effect_percentile": float((valid <= real_g).mean()),
                                "empirical_p_one_sided": (1 + int((valid >= real_g).sum())) / 1001,
                                "above_random_95th_percentile": real_g > q95,
                            }
                        )

    effects = pd.DataFrame(effect_rows).sort_values(
        ["analysis_tier", "contrast", "lineage", "program_id", "score_method"]
    )
    benchmark = pd.DataFrame(benchmark_rows).sort_values(
        ["analysis_tier", "contrast", "lineage", "program_id", "score_method"]
    )
    effects.to_csv(output / "alcohol_context_program_effects.csv", index=False)
    benchmark.to_csv(output / "alcohol_context_random_benchmark.csv", index=False)

    joined = effects.merge(
        benchmark[
            [
                "contrast",
                "lineage",
                "program_id",
                "score_method",
                "cell_gate",
                "above_random_95th_percentile",
                "real_effect_percentile",
            ]
        ],
        on=["contrast", "lineage", "program_id", "score_method", "cell_gate"],
        how="left",
    )
    primary = joined[
        joined["analysis_tier"].eq("primary")
        & joined["program_id"].isin(PRIMARY_PROGRAMS)
    ]
    summary_rows = []
    for program_id, data in primary.groupby("program_id"):
        target = data[data["contrast"].eq("alcohol_cirrhosis_vs_alcohol_hepatitis")]
        context_enriched = (
            len(target) == 2
            and target["robust_ci95_low"].gt(0).all()
            and target["above_random_95th_percentile"].astype(bool).all()
        )
        summary_rows.append(
            {
                "program_id": program_id,
                "post_lock_label": (
                    "CIRRHOSIS_CONTEXT_ENRICHED_POST_LOCK"
                    if context_enriched
                    else "NOT_CIRRHOSIS_CONTEXT_ENRICHED_POST_LOCK"
                ),
                "alcohol_cirrhosis_vs_hepatitis_dual_score_positive_ci": bool(
                    len(target) == 2 and target["robust_ci95_low"].gt(0).all()
                ),
                "alcohol_cirrhosis_vs_hepatitis_dual_score_random_specific": bool(
                    len(target) == 2
                    and target["above_random_95th_percentile"].astype(bool).all()
                ),
            }
        )
    summary = pd.DataFrame(summary_rows)
    summary.to_csv(output / "alcohol_context_primary_summary.csv", index=False)
    run = {
        "contrasts": list(CONTRASTS),
        "effect_rows": len(effects),
        "random_benchmark_rows": len(benchmark),
        "primary_programs": sorted(PRIMARY_PROGRAMS),
        "context_enriched_programs": summary.loc[
            summary["post_lock_label"].eq("CIRRHOSIS_CONTEXT_ENRICHED_POST_LOCK"),
            "program_id",
        ].tolist(),
    }
    (output / "alcohol_context_run_summary.json").write_text(
        json.dumps(run, indent=2) + "\n", encoding="utf-8"
    )
    print(summary.to_string(index=False))
    print(json.dumps(run, indent=2))


if __name__ == "__main__":
    main()
