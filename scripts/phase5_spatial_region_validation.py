from __future__ import annotations

import hashlib
import json
from pathlib import Path

import numpy as np
import pandas as pd
from scipy.stats import rankdata


PRIMARY = {
    "Endo2_primary": ["TFF3", "TSPAN5", "PPDPF", "EFEMP1", "NTS", "ADIRF", "LGALS3"],
    "SAEndo1_primary": ["GSN", "RBP7", "PLPP1", "PLVAP", "VWA1"],
}
SECONDARY = {
    "Endo2_secondary": [
        "TFF3", "TSPAN5", "PPDPF", "EFEMP1", "NTS", "ADIRF", "LGALS3",
        "LAPTM5", "TMSB10", "S100A6", "VIM", "S100A10", "CALD1", "ANXA2",
        "GUK1", "C4ORF48", "SNCG",
    ],
    "SAEndo1_secondary": ["GSN", "RBP7", "PLPP1", "PLVAP", "VWA1", "TMEM88"],
}
ITERATIONS = 10_000


def sha256(path: Path) -> str:
    digest = hashlib.sha256()
    with path.open("rb") as handle:
        for block in iter(lambda: handle.read(1024 * 1024), b""):
            digest.update(block)
    return digest.hexdigest().upper()


def parse_samples(columns: list[str]) -> pd.DataFrame:
    rows = []
    for column in columns:
        prefix = column.removesuffix(" Average")
        donor, region = prefix.split("-", 1)
        rows.append(
            {
                "column": column,
                "donor_id": donor,
                "region": {"fibrosis": "scar", "parenchyma": "parenchyma"}[region],
            }
        )
    result = pd.DataFrame(rows)
    counts = result.groupby("donor_id")["region"].nunique()
    if len(result) != 16 or len(counts) != 8 or not counts.eq(2).all():
        raise RuntimeError("expected eight donors with paired scar and parenchymal profiles")
    return result


def matched_pool(log_expression: pd.DataFrame, excluded: set[str]) -> tuple[pd.Series, pd.Series]:
    detection = log_expression.gt(0).sum(axis=1)
    abundance = log_expression.mean(axis=1)
    eligible = (~log_expression.index.isin(excluded)) & detection.gt(0) & abundance.notna()
    detection_bin = detection.astype(int).astype(str)
    abundance_bin = pd.Series("", index=log_expression.index, dtype=object)
    abundance_bin.loc[eligible] = pd.qcut(
        abundance.loc[eligible], q=10, duplicates="drop", labels=False
    ).astype(str)
    return detection_bin, abundance_bin


def choose_matched_indices(
    rng: np.random.Generator,
    target_genes: list[str],
    exact_pools: dict[str, np.ndarray],
    widened_pools: dict[str, np.ndarray],
    fallback_pool: np.ndarray,
) -> np.ndarray:
    selected: list[int] = []
    used: set[int] = set()
    for gene in target_genes:
        candidates = exact_pools[gene]
        if len(candidates) <= len(used):
            # Deterministic widening only when an exact stratum is exhausted.
            candidates = widened_pools[gene]
        if len(candidates) <= len(used):
            candidates = fallback_pool
        choice = int(rng.choice(candidates))
        while choice in used:
            choice = int(rng.choice(candidates))
        selected.append(choice)
        used.add(choice)
    return np.asarray(selected, dtype=int)


def paired_effects(scores: np.ndarray, samples: pd.DataFrame) -> pd.DataFrame:
    rows = []
    for donor, donor_samples in samples.groupby("donor_id", sort=True):
        scar_index = int(donor_samples.index[donor_samples["region"].eq("scar")][0])
        par_index = int(donor_samples.index[donor_samples["region"].eq("parenchyma")][0])
        rows.append(
            {
                "donor_id": donor,
                "scar_score": float(scores[scar_index]),
                "parenchyma_score": float(scores[par_index]),
                "scar_minus_parenchyma": float(scores[scar_index] - scores[par_index]),
            }
        )
    return pd.DataFrame(rows)


def main() -> None:
    repo = Path(__file__).resolve().parents[1]
    source = repo / "data" / "external" / "phase5_spatial" / "HEP4-6-2538-s001.xlsx"
    output = repo / "results" / "phase5"
    output.mkdir(parents=True, exist_ok=True)

    raw = pd.read_excel(source, sheet_name="gene expr")
    raw["Gene"] = raw["Gene"].astype(str).str.upper()
    raw = raw.drop_duplicates("Gene").set_index("Gene")
    expression = raw.apply(pd.to_numeric, errors="coerce").fillna(0.0)
    samples = parse_samples(expression.columns.tolist())
    samples.index = np.arange(len(samples))
    donor_order = sorted(samples["donor_id"].unique())
    scar_indices = np.asarray(
        [samples.index[(samples["donor_id"].eq(donor)) & (samples["region"].eq("scar"))][0] for donor in donor_order],
        dtype=int,
    )
    parenchyma_indices = np.asarray(
        [samples.index[(samples["donor_id"].eq(donor)) & (samples["region"].eq("parenchyma"))][0] for donor in donor_order],
        dtype=int,
    )
    log_expression = np.log2(expression + 1)
    values = log_expression.to_numpy(float)

    gene_mean = values.mean(axis=1)
    gene_sd = values.std(axis=1, ddof=1)
    safe_sd = np.where(np.isfinite(gene_sd) & (gene_sd > 0), gene_sd, 1.0)
    standardized = (values - gene_mean[:, None]) / safe_sd[:, None]
    standardized[~np.isfinite(gene_sd) | (gene_sd == 0), :] = 0.0
    normalized_ranks = (rankdata(values, axis=0, method="average") - 1) / (len(raw) - 1) - 0.5

    all_sets = {**PRIMARY, **SECONDARY}
    excluded = {gene for genes in all_sets.values() for gene in genes}
    detection_bin, abundance_bin = matched_pool(log_expression, excluded)
    gene_to_index = {gene: index for index, gene in enumerate(log_expression.index)}
    rng = np.random.default_rng(20260901)

    coverage_rows = []
    score_rows = []
    benchmark_rows = []
    summary_rows = []
    for set_id, requested in all_sets.items():
        measured = [gene for gene in requested if gene in gene_to_index]
        coverage = len(measured) / len(requested)
        required = 0.80
        coverage_rows.append(
            {
                "spatial_study": "Chung_2022_PMC9426406",
                "gene_set_id": set_id,
                "requested_genes": len(requested),
                "measured_genes": len(measured),
                "coverage_fraction": coverage,
                "coverage_gate": required,
                "coverage_pass": coverage >= required,
                "measured_gene_symbols": ";".join(measured),
                "missing_gene_symbols": ";".join(sorted(set(requested) - set(measured))),
            }
        )
        if coverage < required:
            continue
        indices = np.asarray([gene_to_index[gene] for gene in measured], dtype=int)
        eligible_background = ~log_expression.index.isin(excluded)
        exact_pools = {
            gene: np.flatnonzero(
                (
                    detection_bin.eq(detection_bin.loc[gene])
                    & abundance_bin.eq(abundance_bin.loc[gene])
                    & eligible_background
                ).to_numpy()
            )
            for gene in measured
        }
        widened_pools = {
            gene: np.flatnonzero(
                (detection_bin.eq(detection_bin.loc[gene]) & eligible_background).to_numpy()
            )
            for gene in measured
        }
        fallback_pool = np.flatnonzero(eligible_background)
        real_scores = {
            "standardized_mean": standardized[indices].mean(axis=0),
            "singscore": normalized_ranks[indices].mean(axis=0),
        }
        null_medians = {
            "standardized_mean": np.empty(ITERATIONS, dtype=float),
            "singscore": np.empty(ITERATIONS, dtype=float),
        }
        for iteration in range(ITERATIONS):
            random_indices = choose_matched_indices(
                rng,
                measured,
                exact_pools,
                widened_pools,
                fallback_pool,
            )
            for method, matrix in (
                ("standardized_mean", standardized),
                ("singscore", normalized_ranks),
            ):
                random_scores = matrix[random_indices].mean(axis=0)
                null_medians[method][iteration] = float(
                    np.median(random_scores[scar_indices] - random_scores[parenchyma_indices])
                )

        method_passes = []
        for method, scores in real_scores.items():
            effects = paired_effects(scores, samples)
            effects.insert(0, "score_method", method)
            effects.insert(0, "gene_set_id", set_id)
            score_rows.append(effects)
            observed = float(effects["scar_minus_parenchyma"].median())
            positive_fraction = float(effects["scar_minus_parenchyma"].gt(0).mean())
            null = null_medians[method]
            q95 = float(np.quantile(null, 0.95))
            passed = observed > 0 and positive_fraction >= 0.75 and observed > q95
            method_passes.append(passed)
            benchmark_rows.append(
                {
                    "spatial_study": "Chung_2022_PMC9426406",
                    "gene_set_id": set_id,
                    "score_method": method,
                    "independent_donors": len(effects),
                    "positive_donors": int(effects["scar_minus_parenchyma"].gt(0).sum()),
                    "positive_donor_fraction": positive_fraction,
                    "observed_median_scar_minus_parenchyma": observed,
                    "random_median": float(np.median(null)),
                    "random_95th_percentile": q95,
                    "observed_percentile": float((null <= observed).mean()),
                    "empirical_p_one_sided": (1 + int((null >= observed).sum())) / (ITERATIONS + 1),
                    "descriptive_pass_if_region_averages_substituted": passed,
                }
            )
        summary_rows.append(
            {
                "spatial_study": "Chung_2022_PMC9426406",
                "gene_set_id": set_id,
                "dual_score_descriptive_pass": bool(all(method_passes)),
                "strict_frozen_spatial_pass": False,
                "strict_status": "SPATIAL_RESOURCE_PARTIALLY_EVALUABLE",
                "reason": (
                    "Author-processed donor-by-region average expression was public, but spot-level "
                    "matrices and labels required by the frozen rule were not public."
                ),
            }
        )

    scores = pd.concat(score_rows, ignore_index=True)
    coverage_table = pd.DataFrame(coverage_rows)
    benchmark = pd.DataFrame(benchmark_rows)
    summary = pd.DataFrame(summary_rows)
    coverage_table.to_csv(output / "spatial_2022_gene_set_coverage.csv", index=False)
    scores.to_csv(output / "spatial_2022_donor_region_scores.csv", index=False)
    benchmark.to_csv(output / "spatial_2022_random_benchmark.csv", index=False)
    summary.to_csv(output / "spatial_2022_evaluability_summary.csv", index=False)
    run = {
        "source": str(source.relative_to(repo)).replace("\\", "/"),
        "source_sha256": sha256(source),
        "matrix_level": "author_processed_donor_by_region_average",
        "independent_donors": 8,
        "regions": ["parenchyma", "scar"],
        "random_iterations": ITERATIONS,
        "strict_frozen_status": "SPATIAL_RESOURCE_PARTIALLY_EVALUABLE",
        "primary_descriptive_dual_score_pass": summary.loc[
            summary["gene_set_id"].isin(PRIMARY),
            ["gene_set_id", "dual_score_descriptive_pass"],
        ].set_index("gene_set_id")["dual_score_descriptive_pass"].to_dict(),
    }
    (output / "spatial_2022_run_summary.json").write_text(
        json.dumps(run, indent=2) + "\n", encoding="utf-8"
    )
    print(summary.to_string(index=False))
    print(benchmark[benchmark["gene_set_id"].isin(PRIMARY)].to_string(index=False))


if __name__ == "__main__":
    main()
