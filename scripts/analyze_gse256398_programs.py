from __future__ import annotations

import hashlib
import itertools
import json
from pathlib import Path

import numpy as np
import pandas as pd
from scipy.io import mmread
from scipy.stats import rankdata, spearmanr

from analyze_gse202379_programs import effect_statistics, singscore_up, stable_seed
from audit_gse256398_gates import CONTRASTS


def benjamini_hochberg(values: pd.Series) -> pd.Series:
    values = values.astype(float)
    order = np.argsort(values.to_numpy())
    ranked = values.to_numpy()[order]
    adjusted = ranked * len(ranked) / np.arange(1, len(ranked) + 1)
    adjusted = np.minimum.accumulate(adjusted[::-1])[::-1]
    result = np.empty_like(adjusted)
    result[order] = np.minimum(adjusted, 1.0)
    return pd.Series(result, index=values.index)


def exact_ordinal_permutation_p(values: np.ndarray, labels: np.ndarray) -> tuple[float, str]:
    observed = float(spearmanr(labels, values).statistic)
    n = len(labels)
    n_zero = int((labels == 0).sum())
    n_one = int((labels == 1).sum())
    extreme = 0
    allocations = 0
    indices = set(range(n))
    for zero_indices in itertools.combinations(range(n), n_zero):
        remaining = sorted(indices - set(zero_indices))
        for one_indices in itertools.combinations(remaining, n_one):
            permuted = np.full(n, 2.0)
            permuted[list(zero_indices)] = 0
            permuted[list(one_indices)] = 1
            statistic = float(spearmanr(permuted, values).statistic)
            extreme += abs(statistic) >= abs(observed) - 1e-15
            allocations += 1
    return extreme / allocations, f"exact_{allocations}_allocations"


def main() -> None:
    repo = Path(__file__).resolve().parents[1]
    interim = repo / "data" / "interim" / "GSE256398"
    genes = pd.read_csv(interim / "genes.csv")["gene"].astype(str)
    manifest = pd.read_csv(interim / "donor_lineage_manifest.csv")
    manifest["eligible_30"] = manifest["n_cells"].ge(30)
    manifest["eligible_20"] = manifest["n_cells"].ge(20)
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
    coverage = pd.read_csv(repo / "results" / "qc" / "gse256398_program_coverage.csv")
    coverage_lookup = coverage.set_index("program_id").to_dict("index")
    score_rows: list[dict[str, object]] = []
    score_arrays: dict[tuple[str, str], np.ndarray] = {}

    for lineage, lineage_manifest in manifest.groupby("harmonized_lineage", sort=False):
        columns = lineage_manifest.index.to_numpy()
        reference_columns = lineage_manifest.index[lineage_manifest["eligible_20"]].to_numpy()
        means = log_cpm[:, reference_columns].mean(axis=1)
        standard_deviations = log_cpm[:, reference_columns].std(axis=1, ddof=1)
        invariant = standard_deviations == 0
        safe_standard_deviations = standard_deviations.copy()
        safe_standard_deviations[invariant] = 1
        standardized = (log_cpm[:, columns] - means[:, None]) / safe_standard_deviations[:, None]
        standardized[invariant, :] = 0
        for program_id, rows in programs[programs["cell_lineage"].eq(lineage)].groupby("program_id"):
            coverage_info = coverage_lookup[program_id]
            if coverage_info["coverage_tier"] == "not_evaluated":
                continue
            measured_genes = sorted(set(rows["gene_symbol"].astype(str).str.upper()) & set(gene_to_index))
            gene_indices = np.array([gene_to_index[gene] for gene in measured_genes], dtype=int)
            sing = singscore_up(ranks[:, columns], gene_indices)
            zmean = standardized[gene_indices, :].mean(axis=0)
            for method, values in (("singscore", sing), ("standardized_mean", zmean)):
                score_arrays[(program_id, method)] = np.full(len(manifest), np.nan)
                score_arrays[(program_id, method)][columns] = values
                for local_index, manifest_index in enumerate(columns):
                    score_rows.append(
                        {
                            "dataset_id": "GSE256398_human",
                            "program_id": program_id,
                            "lineage": lineage,
                            "score_method": method,
                            "score": float(values[local_index]),
                            "measured_program_genes": len(measured_genes),
                            "program_coverage": coverage_info["coverage"],
                            "coverage_tier": coverage_info["coverage_tier"],
                            **manifest.loc[manifest_index].to_dict(),
                        }
                    )

    output = repo / "results" / "phase3"
    output.mkdir(parents=True, exist_ok=True)
    scores = pd.DataFrame(score_rows)
    score_path = output / "gse256398_donor_program_scores.csv.gz"
    scores.to_csv(score_path, index=False, compression="gzip")

    effect_rows: list[dict[str, object]] = []
    for contrast, (controls, cases) in CONTRASTS.items():
        group = manifest["disease_group"].map(
            {**{value: "control" for value in controls}, **{value: "case" for value in cases}}
        ).fillna("excluded")
        for program_id, program_rows in programs.groupby("program_id"):
            lineage = str(program_rows["cell_lineage"].iloc[0])
            coverage_info = coverage_lookup[program_id]
            selected = (
                manifest["harmonized_lineage"].eq(lineage)
                & manifest["eligible_30"]
                & group.ne("excluded")
            )
            selected_indices = manifest.index[selected]
            selected_groups = group.loc[selected_indices]
            for score_method in ("singscore", "standardized_mean"):
                values = score_arrays[(program_id, score_method)][selected_indices]
                control = values[selected_groups.eq("control").to_numpy()]
                case = values[selected_groups.eq("case").to_numpy()]
                if len(control) < 3 or len(case) < 3:
                    raise RuntimeError("prequalified 30-cell contrast lost donor support")
                effect_rows.append(
                    {
                        "dataset_id": "GSE256398_human",
                        "contrast": contrast,
                        "program_id": program_id,
                        "lineage": lineage,
                        "score_method": score_method,
                        "analysis_tier": "post_lock_external",
                        "cell_gate": 30,
                        "coverage_tier": coverage_info["coverage_tier"],
                        "program_coverage": coverage_info["coverage"],
                        "expected_direction": "positive_case_higher",
                        **effect_statistics(
                            control,
                            case,
                            stable_seed("GSE256398_human", contrast, program_id, score_method),
                        ),
                    }
                )
    effects = pd.DataFrame(effect_rows).sort_values(["contrast", "lineage", "program_id", "score_method"])
    effect_path = output / "gse256398_program_effects.csv"
    effects.to_csv(effect_path, index=False)

    trend_rows: list[dict[str, object]] = []
    metabolic = manifest["metabolic_order"].notna() & manifest["eligible_30"]
    for program_id, program_rows in programs.groupby("program_id"):
        lineage = str(program_rows["cell_lineage"].iloc[0])
        selected_indices = manifest.index[metabolic & manifest["harmonized_lineage"].eq(lineage)]
        labels = manifest.loc[selected_indices, "metabolic_order"].to_numpy(dtype=float)
        if list(pd.Series(labels).value_counts().sort_index()) != [3, 4, 4]:
            raise RuntimeError("metabolic ordinal analysis does not contain frozen 3/4/4 donors")
        for score_method in ("singscore", "standardized_mean"):
            values = score_arrays[(program_id, score_method)][selected_indices]
            rho = float(spearmanr(labels, values).statistic)
            p_value, mode = exact_ordinal_permutation_p(values, labels)
            trend_rows.append(
                {
                    "dataset_id": "GSE256398_human",
                    "trend": "masld_f0_to_mash_fibrosis_to_mash_cirrhosis",
                    "program_id": program_id,
                    "lineage": lineage,
                    "score_method": score_method,
                    "analysis_tier": "post_lock_external",
                    "cell_gate": 30,
                    "n_f0": int((labels == 0).sum()),
                    "n_fibrosis": int((labels == 1).sum()),
                    "n_cirrhosis": int((labels == 2).sum()),
                    "spearman_rho": rho,
                    "permutation_p_two_sided": p_value,
                    "permutation_mode": mode,
                }
            )
    trends = pd.DataFrame(trend_rows)
    trends["permutation_fdr_within_lineage_method"] = trends.groupby(
        ["lineage", "score_method"], group_keys=False
    )["permutation_p_two_sided"].apply(benjamini_hochberg)
    trend_path = output / "gse256398_metabolic_ordinal_trends.csv"
    trends.sort_values(["lineage", "program_id", "score_method"]).to_csv(trend_path, index=False)

    summary = {
        "dataset_id": "GSE256398_human",
        "normalization": "log2(CPM + 1) within donor-lineage pseudobulk",
        "score_methods": ["direction-aware singscore", "standardized signed mean"],
        "standardization_reference": "all lineage donors passing 20-cell gate; outcome-blind",
        "binary_effect_rows": len(effects),
        "ordinal_trend_rows": len(trends),
        "ordinal_permutation": "exact unique 3/4/4 label allocations",
        "seed": 20260830,
        "score_sha256": hashlib.sha256(score_path.read_bytes()).hexdigest().upper(),
        "effect_sha256": hashlib.sha256(effect_path.read_bytes()).hexdigest().upper(),
        "trend_sha256": hashlib.sha256(trend_path.read_bytes()).hexdigest().upper(),
    }
    (repo / "results" / "logs" / "gse256398_analysis_run.json").write_text(
        json.dumps(summary, indent=2) + "\n", encoding="utf-8"
    )
    print(effects.groupby("contrast").agg(rows=("program_id", "size"), positive_g=("hedges_g", lambda x: int((x > 0).sum())), ci_excludes_zero=("robust_ci95_low", lambda x: int((x > 0).sum()))).to_string())
    print(trends.groupby(["lineage", "score_method"]).agg(rows=("program_id", "size"), positive_rho=("spearman_rho", lambda x: int((x > 0).sum())), fdr_lt_005=("permutation_fdr_within_lineage_method", lambda x: int((x < 0.05).sum()))).to_string())
    print(json.dumps(summary, indent=2))


if __name__ == "__main__":
    main()
