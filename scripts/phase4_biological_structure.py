from __future__ import annotations

import hashlib
import json
import math
from itertools import combinations, product
from pathlib import Path

import numpy as np
import pandas as pd
import statsmodels.api as sm
from scipy.io import mmread
from scipy.optimize import minimize_scalar
from scipy.stats import norm, rankdata, spearmanr

from audit_gse202379_gates import contrast_groups as gse202379_groups


SEED = 20260831
PERMUTATIONS = 10_000
RANDOM_MODULES = 10_000
PROGRAMS_OF_INTEREST = ("RAM2019_ENDO_2", "RAM2019_ENDO_6_SAENDO1")


def stable_seed(*parts: str) -> int:
    payload = "|||".join(parts).encode("utf-8")
    offset = int.from_bytes(hashlib.sha256(payload).digest()[:8], "little")
    return (SEED + offset) % (2**32)


def bh_fdr(values: pd.Series) -> pd.Series:
    p = pd.to_numeric(values, errors="coerce").to_numpy(dtype=float)
    result = np.full(len(p), np.nan)
    finite = np.isfinite(p)
    if not finite.any():
        return pd.Series(result, index=values.index)
    observed = p[finite]
    order = np.argsort(observed)
    ranked = observed[order]
    adjusted = ranked * len(ranked) / np.arange(1, len(ranked) + 1)
    adjusted = np.minimum.accumulate(adjusted[::-1])[::-1]
    restored = np.empty_like(adjusted)
    restored[order] = np.minimum(adjusted, 1.0)
    result[np.flatnonzero(finite)] = restored
    return pd.Series(result, index=values.index)


def hedges_g(control: np.ndarray, case: np.ndarray) -> float:
    control = np.asarray(control, dtype=float)
    case = np.asarray(case, dtype=float)
    n_control, n_case = len(control), len(case)
    df = n_control + n_case - 2
    if n_control < 2 or n_case < 2 or df <= 0:
        return float("nan")
    pooled = ((n_control - 1) * control.var(ddof=1) + (n_case - 1) * case.var(ddof=1)) / df
    if pooled <= 0 or not np.isfinite(pooled):
        return float("nan")
    correction = 1 - 3 / (4 * df - 1)
    return float(correction * (case.mean() - control.mean()) / math.sqrt(pooled))


def vectorized_hedges_g(control: np.ndarray, case: np.ndarray) -> tuple[np.ndarray, np.ndarray]:
    n_control, n_case = control.shape[1], case.shape[1]
    df = n_control + n_case - 2
    pooled = ((n_control - 1) * control.var(axis=1, ddof=1) + (n_case - 1) * case.var(axis=1, ddof=1)) / df
    correction = 1 - 3 / (4 * df - 1)
    with np.errstate(divide="ignore", invalid="ignore"):
        effect = correction * (case.mean(axis=1) - control.mean(axis=1)) / np.sqrt(pooled)
    effect[~np.isfinite(effect)] = np.nan
    variance = (n_control + n_case) / (n_control * n_case) + effect**2 / (2 * df)
    variance[~np.isfinite(variance) | (variance <= 0)] = np.nan
    return effect, variance


def reml_meta(effect: np.ndarray, variance: np.ndarray) -> dict[str, float]:
    effect = np.asarray(effect, dtype=float)
    variance = np.asarray(variance, dtype=float)
    keep = np.isfinite(effect) & np.isfinite(variance) & (variance > 0)
    effect, variance = effect[keep], variance[keep]
    k = len(effect)
    if k < 2:
        return {key: float("nan") for key in (
            "fixed_effect", "fixed_se", "fixed_low", "fixed_high", "fixed_p",
            "random_effect", "random_se", "random_low", "random_high", "tau2", "q", "i2",
        )}
    fixed_weights = 1 / variance
    fixed = float(np.sum(fixed_weights * effect) / np.sum(fixed_weights))
    fixed_se = float(math.sqrt(1 / np.sum(fixed_weights)))
    q = float(np.sum(fixed_weights * (effect - fixed) ** 2))
    i2 = float(max(0.0, (q - (k - 1)) / q)) if q > 0 else 0.0

    def objective(tau2: float) -> float:
        weights = 1 / (variance + tau2)
        mean = np.sum(weights * effect) / np.sum(weights)
        return float(np.sum(np.log(variance + tau2)) + np.log(np.sum(weights)) + np.sum(weights * (effect - mean) ** 2))

    upper = max(1.0, float(np.var(effect, ddof=1) * 20 + np.max(variance)))
    optimized = minimize_scalar(objective, bounds=(0.0, upper), method="bounded")
    tau2 = float(max(0.0, optimized.x if optimized.success else 0.0))
    random_weights = 1 / (variance + tau2)
    random_effect = float(np.sum(random_weights * effect) / np.sum(random_weights))
    random_se = float(math.sqrt(1 / np.sum(random_weights)))
    return {
        "fixed_effect": fixed,
        "fixed_se": fixed_se,
        "fixed_low": fixed - 1.96 * fixed_se,
        "fixed_high": fixed + 1.96 * fixed_se,
        "fixed_p": float(2 * norm.sf(abs(fixed / fixed_se))),
        "random_effect": random_effect,
        "random_se": random_se,
        "random_low": random_effect - 1.96 * random_se,
        "random_high": random_effect + 1.96 * random_se,
        "tau2": tau2,
        "q": q,
        "i2": i2,
    }


def load_gene_context(
    repo: Path,
    dataset_id: str,
    interim_name: str,
    contrast: str,
    group_labels: pd.Series,
) -> pd.DataFrame:
    interim = repo / "data" / "interim" / interim_name
    genes = pd.read_csv(interim / "genes.csv")["gene"].astype(str).str.upper()
    manifest = pd.read_csv(interim / "donor_lineage_manifest.csv")
    counts = mmread(interim / "donor_lineage_raw_counts.mtx").tocsc()
    selected = manifest["harmonized_lineage"].eq("endothelial") & manifest["n_cells"].ge(30) & group_labels.ne("excluded")
    indices = manifest.index[selected].to_numpy()
    labels = group_labels.loc[indices]
    controls = indices[labels.eq("control").to_numpy()]
    cases = indices[labels.eq("case").to_numpy()]
    if len(controls) < 3 or len(cases) < 3:
        raise RuntimeError(f"{dataset_id} {contrast}: insufficient endothelial donors")
    library = np.asarray(counts.sum(axis=0)).ravel()
    log_cpm = np.log2(counts.toarray() / library * 1_000_000 + 1)
    effects, variances = vectorized_hedges_g(log_cpm[:, controls], log_cpm[:, cases])
    detection = np.asarray((counts[:, indices] > 0).mean(axis=1)).ravel()
    frame = pd.DataFrame({
        "gene_symbol": genes,
        f"g__{dataset_id}": effects,
        f"var__{dataset_id}": variances,
        f"detection__{dataset_id}": detection,
        f"n_control__{dataset_id}": len(controls),
        f"n_case__{dataset_id}": len(cases),
    })
    frame["abs_g"] = frame[f"g__{dataset_id}"].abs()
    frame = frame.sort_values(["gene_symbol", "abs_g"], ascending=[True, False]).drop_duplicates("gene_symbol")
    return frame.drop(columns="abs_g")


def analysis_a_member_gene_coherence(repo: Path, output: Path) -> dict[str, object]:
    m202 = pd.read_csv(repo / "data" / "interim" / "GSE202379" / "donor_lineage_manifest.csv")
    m244 = pd.read_csv(repo / "data" / "interim" / "GSE244832" / "donor_lineage_manifest.csv")
    m256 = pd.read_csv(repo / "data" / "interim" / "GSE256398" / "donor_lineage_manifest.csv")
    frames = [
        load_gene_context(
            repo,
            "GSE202379",
            "GSE202379",
            "clinical_cirrhosis_vs_healthy",
            gse202379_groups(m202)["clinical_cirrhosis_vs_healthy"],
        ),
        load_gene_context(
            repo,
            "GSE244832",
            "GSE244832",
            "mash_f2f4_group_vs_normal_sensitivity",
            m244["disease_group"].map({"normal": "control", "MASH": "case"}).fillna("excluded"),
        ),
        load_gene_context(
            repo,
            "GSE256398",
            "GSE256398",
            "mash_cirrhosis_vs_healthy",
            m256["disease_group"].map({"healthy": "control", "mash_cirrhosis": "case"}).fillna("excluded"),
        ),
    ]
    merged = frames[0]
    for frame in frames[1:]:
        merged = merged.merge(frame, on="gene_symbol", how="inner")
    datasets = ("GSE202379", "GSE244832", "GSE256398")
    for dataset in datasets:
        merged = merged[merged[f"detection__{dataset}"].ge(0.20) & merged[f"g__{dataset}"].notna()].copy()
    merged["average_detection"] = merged[[f"detection__{x}" for x in datasets]].mean(axis=1)
    merged["detection_decile"] = pd.qcut(
        merged["average_detection"].rank(method="first"), 10, labels=False
    ).astype(int)
    programs = pd.read_csv(repo / "literature" / "program_inventory.csv")
    inventory = {
        program: set(programs.loc[programs["program_id"].eq(program), "gene_symbol"].astype(str).str.upper())
        for program in PROGRAMS_OF_INTEREST
    }
    all_gene_rows: list[dict[str, object]] = []
    summary_rows: list[dict[str, object]] = []
    loco_rows: list[dict[str, object]] = []
    random_rows: list[dict[str, object]] = []
    rng = np.random.default_rng(stable_seed("phase4", "member_gene_random"))
    gcols = [f"g__{dataset}" for dataset in datasets]
    vcols = [f"var__{dataset}" for dataset in datasets]
    for program_id, members in inventory.items():
        data = merged[merged["gene_symbol"].isin(members)].copy()
        if data.empty:
            raise RuntimeError(f"no eligible genes for {program_id}")
        program_records: list[dict[str, object]] = []
        for row in data.itertuples(index=False):
            effects = np.array([getattr(row, column) for column in gcols], dtype=float)
            variances = np.array([getattr(row, column) for column in vcols], dtype=float)
            meta = reml_meta(effects, variances)
            record = {
                "program_id": program_id,
                "gene_symbol": row.gene_symbol,
                **{column: getattr(row, column) for column in gcols + vcols},
                "average_detection": row.average_detection,
                "coherent_positive_member": bool((effects > 0).all()),
                "minimum_g": float(np.min(effects)),
                "median_g": float(np.median(effects)),
                **meta,
            }
            program_records.append(record)
        program_frame = pd.DataFrame(program_records)
        program_frame["fixed_fdr_within_program"] = bh_fdr(program_frame["fixed_p"])
        program_frame["meta_supported_member"] = (
            program_frame["coherent_positive_member"]
            & program_frame["fixed_low"].gt(0)
            & program_frame["fixed_fdr_within_program"].lt(0.05)
        )
        all_gene_rows.extend(program_frame.to_dict("records"))

        real_count = int(program_frame["coherent_positive_member"].sum())
        real_meta_count = int(program_frame["meta_supported_member"].sum())
        target_deciles = data["detection_decile"].value_counts().sort_index().to_dict()
        random_counts = np.empty(RANDOM_MODULES, dtype=int)
        pool = merged[~merged["gene_symbol"].isin(members)].copy()
        for iteration in range(RANDOM_MODULES):
            sampled_parts = []
            for decile, size in target_deciles.items():
                candidates = pool[pool["detection_decile"].eq(decile)]
                if candidates.empty:
                    candidates = pool
                take = rng.choice(candidates.index.to_numpy(), size=size, replace=len(candidates) < size)
                sampled_parts.append(merged.loc[take])
            sampled = pd.concat(sampled_parts, ignore_index=True)
            random_counts[iteration] = int((sampled[gcols].to_numpy(dtype=float) > 0).all(axis=1).sum())
        empirical_p = float((1 + np.sum(random_counts >= real_count)) / (RANDOM_MODULES + 1))
        summary_rows.append({
            "program_id": program_id,
            "published_program_genes": len(members),
            "eligible_shared_genes": len(program_frame),
            "coherent_positive_genes": real_count,
            "coherent_positive_fraction": real_count / len(program_frame),
            "meta_supported_genes": real_meta_count,
            "random_coherent_median": float(np.median(random_counts)),
            "random_coherent_95th": float(np.quantile(random_counts, 0.95)),
            "random_coherent_99th": float(np.quantile(random_counts, 0.99)),
            "empirical_p_coherent_count": empirical_p,
            "above_random_95th": bool(real_count > np.quantile(random_counts, 0.95)),
        })
        for held_out in datasets:
            discovery = [dataset for dataset in datasets if dataset != held_out]
            selected = data[(data[[f"g__{dataset}" for dataset in discovery]] > 0).all(axis=1)]
            held_values = selected[f"g__{held_out}"].to_numpy(dtype=float)
            loco_rows.append({
                "program_id": program_id,
                "held_out_dataset": held_out,
                "discovery_datasets": ";".join(discovery),
                "genes_positive_in_both_discovery": len(selected),
                "held_out_positive_genes": int((held_values > 0).sum()),
                "held_out_sign_retention": float((held_values > 0).mean()) if len(held_values) else np.nan,
                "held_out_median_g": float(np.median(held_values)) if len(held_values) else np.nan,
            })
        random_rows.extend(
            {"program_id": program_id, "iteration": i + 1, "coherent_positive_genes": int(value)}
            for i, value in enumerate(random_counts)
        )
    gene_frame = pd.DataFrame(all_gene_rows).sort_values(["program_id", "fixed_effect"], ascending=[True, False])
    summary = pd.DataFrame(summary_rows)
    loco = pd.DataFrame(loco_rows)
    random_frame = pd.DataFrame(random_rows)
    gene_frame.to_csv(output / "endothelial_member_gene_coherence.csv", index=False)
    summary.to_csv(output / "endothelial_member_gene_summary.csv", index=False)
    loco.to_csv(output / "endothelial_member_gene_loco.csv", index=False)
    random_frame.to_csv(output / "endothelial_member_gene_random_controls.csv.gz", index=False, compression="gzip")
    return {
        "shared_gene_universe": len(merged),
        "program_summaries": summary.to_dict("records"),
    }


def residualize_program_scores(frame: pd.DataFrame, technical: bool) -> pd.DataFrame:
    rows: list[pd.DataFrame] = []
    for (_, _), data in frame.groupby(["program_id", "score_method"], sort=False):
        data = data.copy()
        y = data["score"].to_numpy(dtype=float)
        groups = pd.get_dummies(data["disease_group"].astype(str), drop_first=True, dtype=float)
        design = groups
        if technical:
            technical_frame = pd.DataFrame({
                "log10_n_cells": np.log10(data["n_cells"].astype(float).to_numpy() + 1),
                "log10_library_size": np.log10(data["library_size"].astype(float).to_numpy() + 1),
            }, index=data.index)
            design = pd.concat([design, technical_frame], axis=1)
        design = sm.add_constant(design.astype(float), has_constant="add")
        fit = sm.OLS(y, design).fit()
        data["residual_score"] = fit.resid
        rows.append(data)
    return pd.concat(rows, ignore_index=True)


def permutation_spearman(x: np.ndarray, y: np.ndarray, seed: int) -> tuple[float, float]:
    observed = float(spearmanr(x, y).statistic)
    xr = rankdata(x).astype(float)
    yr = rankdata(y).astype(float)
    xr -= xr.mean()
    yr -= yr.mean()
    denominator = math.sqrt(float(np.dot(xr, xr) * np.dot(yr, yr)))
    if denominator == 0 or not np.isfinite(observed):
        return float("nan"), float("nan")
    rng = np.random.default_rng(seed)
    extreme = 0
    for _ in range(PERMUTATIONS):
        permuted = rng.permutation(yr)
        value = float(np.dot(xr, permuted) / denominator)
        extreme += abs(value) >= abs(observed) - 1e-15
    return observed, float((extreme + 1) / (PERMUTATIONS + 1))


def cross_lineage_pairs(program_lineages: dict[str, str]) -> list[tuple[str, str]]:
    programs = sorted(program_lineages)
    return [
        (left, right)
        for left, right in combinations(programs, 2)
        if program_lineages[left] != program_lineages[right]
    ]


def analysis_b_cross_lineage_coupling(repo: Path, output: Path) -> dict[str, object]:
    paths = {
        "GSE244832": repo / "results" / "sensitivity" / "gse244832_donor_program_scores.csv.gz",
        "GSE256398": repo / "results" / "phase3" / "gse256398_donor_program_scores.csv.gz",
    }
    programs = pd.read_csv(repo / "literature" / "program_inventory.csv")
    program_lineages = programs.drop_duplicates("program_id").set_index("program_id")["cell_lineage"].to_dict()
    pairs = cross_lineage_pairs(program_lineages)
    correlation_rows: list[dict[str, object]] = []
    for dataset, path in paths.items():
        scores = pd.read_csv(path)
        eligibility_column = "eligible_30"
        scores = scores[scores[eligibility_column].astype(bool)].copy()
        scores["donor_id"] = scores["donor_id"].astype(str)
        for adjustment, technical in (("disease_group_plus_qc", True), ("disease_group_only", False)):
            residuals = residualize_program_scores(scores, technical=technical)
            for method in sorted(residuals["score_method"].unique()):
                method_data = residuals[residuals["score_method"].eq(method)]
                pivot = method_data.pivot_table(index="donor_id", columns="program_id", values="residual_score", aggfunc="first")
                for left, right in pairs:
                    paired = pivot[[left, right]].dropna()
                    if len(paired) < 8:
                        continue
                    if adjustment == "disease_group_plus_qc":
                        rho, p = permutation_spearman(
                            paired[left].to_numpy(dtype=float),
                            paired[right].to_numpy(dtype=float),
                            stable_seed("coupling", dataset, method, left, right),
                        )
                    else:
                        result = spearmanr(paired[left], paired[right])
                        rho, p = float(result.statistic), float(result.pvalue)
                    correlation_rows.append({
                        "dataset_id": dataset,
                        "adjustment": adjustment,
                        "score_method": method,
                        "program_left": left,
                        "lineage_left": program_lineages[left],
                        "program_right": right,
                        "lineage_right": program_lineages[right],
                        "n_donors": len(paired),
                        "spearman_rho": rho,
                        "permutation_p_two_sided": p,
                    })
    correlations = pd.DataFrame(correlation_rows)
    correlations["fdr_within_dataset_method_adjustment"] = correlations.groupby(
        ["dataset_id", "score_method", "adjustment"], group_keys=False
    )["permutation_p_two_sided"].apply(bh_fdr)
    correlations.to_csv(output / "cross_lineage_coupling_correlations.csv", index=False)

    primary = correlations[correlations["adjustment"].eq("disease_group_plus_qc")].copy()
    meta_rows: list[dict[str, object]] = []
    for (left, right, method), data in primary.groupby(["program_left", "program_right", "score_method"]):
        if data["dataset_id"].nunique() != 2:
            continue
        rho = np.clip(data["spearman_rho"].to_numpy(dtype=float), -0.999999, 0.999999)
        weights = data["n_donors"].to_numpy(dtype=float) - 3
        z = np.arctanh(rho)
        combined_z = float(np.sum(weights * z) / np.sum(weights))
        se = float(math.sqrt(1 / np.sum(weights)))
        p = float(2 * norm.sf(abs(combined_z / se)))
        meta_rows.append({
            "program_left": left,
            "lineage_left": data["lineage_left"].iloc[0],
            "program_right": right,
            "lineage_right": data["lineage_right"].iloc[0],
            "score_method": method,
            "cohorts": ";".join(sorted(data["dataset_id"].unique())),
            "meta_spearman_rho": float(np.tanh(combined_z)),
            "meta_z_se": se,
            "meta_p_two_sided": p,
            "same_sign_two_cohorts": bool(np.sign(rho[0]) == np.sign(rho[1])),
            "minimum_absolute_cohort_rho": float(np.min(np.abs(rho))),
        })
    meta = pd.DataFrame(meta_rows)
    meta["meta_fdr_within_method"] = meta.groupby("score_method", group_keys=False)["meta_p_two_sided"].apply(bh_fdr)
    meta.to_csv(output / "cross_lineage_coupling_meta.csv", index=False)

    stable_rows: list[dict[str, object]] = []
    for (left, right), data in primary.groupby(["program_left", "program_right"]):
        if len(data) != 4:
            continue
        meta_pair = meta[(meta["program_left"].eq(left)) & (meta["program_right"].eq(right))]
        if meta_pair["score_method"].nunique() != 2:
            continue
        signs = np.sign(data["spearman_rho"].to_numpy(dtype=float))
        stable = (
            len(set(signs)) == 1
            and int(data["spearman_rho"].abs().ge(0.30).sum()) >= 3
            and meta_pair["meta_fdr_within_method"].lt(0.05).all()
        )
        stable_rows.append({
            "program_left": left,
            "lineage_left": data["lineage_left"].iloc[0],
            "program_right": right,
            "lineage_right": data["lineage_right"].iloc[0],
            "same_sign_all_four": len(set(signs)) == 1,
            "estimates_abs_rho_ge_0_30": int(data["spearman_rho"].abs().ge(0.30).sum()),
            "minimum_absolute_rho": float(data["spearman_rho"].abs().min()),
            "stable_cross_lineage_coupling": bool(stable),
        })
    stable_frame = pd.DataFrame(stable_rows)
    stable_frame.to_csv(output / "cross_lineage_coupling_stability.csv", index=False)
    return {
        "correlation_rows": len(correlations),
        "meta_rows": len(meta),
        "stable_pairs": int(stable_frame["stable_cross_lineage_coupling"].sum()),
    }


def bootstrap_shared_components(data: pd.DataFrame, seed: int) -> tuple[np.ndarray, np.ndarray]:
    groups = {
        group: data[data["disease_group"].eq(group)]["score"].to_numpy(dtype=float)
        for group in ("healthy", "mash_cirrhosis", "alcohol_cirrhosis")
    }
    rng = np.random.default_rng(seed)
    shared = np.empty(PERMUTATIONS, dtype=float)
    divergence = np.empty(PERMUTATIONS, dtype=float)
    for iteration in range(PERMUTATIONS):
        sampled = {
            group: rng.choice(values, size=len(values), replace=True)
            for group, values in groups.items()
        }
        g_mash = hedges_g(sampled["healthy"], sampled["mash_cirrhosis"])
        g_alcohol = hedges_g(sampled["healthy"], sampled["alcohol_cirrhosis"])
        shared[iteration] = (g_mash + g_alcohol) / 2
        divergence[iteration] = (g_mash - g_alcohol) / 2
    return shared[np.isfinite(shared)], divergence[np.isfinite(divergence)]


def analysis_c_etiology_geometry(repo: Path, output: Path) -> dict[str, object]:
    scores = pd.read_csv(repo / "results" / "phase3" / "gse256398_donor_program_scores.csv.gz")
    scores = scores[scores["eligible_30"].astype(bool)].copy()
    effects = pd.read_csv(repo / "results" / "phase3" / "gse256398_program_effects.csv")
    benchmark = pd.read_csv(repo / "results" / "random_controls" / "gse256398_random_module_benchmark.csv")
    effect_lookup = effects.set_index(["program_id", "score_method", "contrast"])["hedges_g"].to_dict()
    random_lookup = benchmark.set_index(["program_id", "score_method", "contrast"])["above_random_95th_percentile"].to_dict()
    rows: list[dict[str, object]] = []
    for (program_id, method), data in scores.groupby(["program_id", "score_method"]):
        relevant = data[data["disease_group"].isin(["healthy", "mash_cirrhosis", "alcohol_cirrhosis"])]
        if relevant["disease_group"].nunique() != 3:
            continue
        g_mash = float(effect_lookup[(program_id, method, "mash_cirrhosis_vs_healthy")])
        g_alcohol = float(effect_lookup[(program_id, method, "alcohol_cirrhosis_vs_healthy")])
        g_direct = float(effect_lookup[(program_id, method, "mash_vs_alcohol_cirrhosis_etiology")])
        shared_observed = (g_mash + g_alcohol) / 2
        divergence_observed = (g_mash - g_alcohol) / 2
        shared_boot, divergence_boot = bootstrap_shared_components(
            relevant, stable_seed("etiology", program_id, method)
        )
        rows.append({
            "program_id": program_id,
            "lineage": relevant["lineage"].iloc[0],
            "score_method": method,
            "g_mash_cirrhosis_vs_healthy": g_mash,
            "g_alcohol_cirrhosis_vs_healthy": g_alcohol,
            "g_mash_vs_alcohol_direct": g_direct,
            "shared_cirrhosis_component": shared_observed,
            "shared_bootstrap_low": float(np.quantile(shared_boot, 0.025)),
            "shared_bootstrap_high": float(np.quantile(shared_boot, 0.975)),
            "etiology_divergence_component": divergence_observed,
            "divergence_bootstrap_low": float(np.quantile(divergence_boot, 0.025)),
            "divergence_bootstrap_high": float(np.quantile(divergence_boot, 0.975)),
            "mash_above_random_95th": bool(random_lookup[(program_id, method, "mash_cirrhosis_vs_healthy")]),
            "alcohol_above_random_95th": bool(random_lookup[(program_id, method, "alcohol_cirrhosis_vs_healthy")]),
            "shared_larger_than_divergence": bool(shared_observed > abs(divergence_observed)),
        })
    method_frame = pd.DataFrame(rows)
    method_frame.to_csv(output / "etiology_shared_divergent_components.csv", index=False)
    summary_rows = []
    for program_id, data in method_frame.groupby("program_id"):
        divergence_nonzero = (
            ((data["divergence_bootstrap_low"] > 0) | (data["divergence_bootstrap_high"] < 0)).all()
            and len(set(np.sign(data["etiology_divergence_component"]))) == 1
        )
        directional = (
            data["g_mash_cirrhosis_vs_healthy"].gt(0).all()
            and data["g_alcohol_cirrhosis_vs_healthy"].gt(0).all()
            and data["shared_larger_than_divergence"].all()
        )
        random_specific = directional and data["mash_above_random_95th"].all() and data["alcohol_above_random_95th"].all()
        summary_rows.append({
            "program_id": program_id,
            "lineage": data["lineage"].iloc[0],
            "shared_directional_backbone": bool(directional),
            "shared_random_specific_backbone": bool(random_specific),
            "etiology_divergent": bool(divergence_nonzero),
            "median_shared_component": float(data["shared_cirrhosis_component"].median()),
            "median_absolute_divergence": float(data["etiology_divergence_component"].abs().median()),
        })
    summary = pd.DataFrame(summary_rows)
    summary.to_csv(output / "etiology_program_classification.csv", index=False)
    return {
        "method_rows": len(method_frame),
        "shared_directional": int(summary["shared_directional_backbone"].sum()),
        "shared_random_specific": int(summary["shared_random_specific_backbone"].sum()),
        "etiology_divergent": int(summary["etiology_divergent"].sum()),
    }


def model_group_effect(data: pd.DataFrame, adjusted: bool) -> dict[str, float]:
    y = data["score"].to_numpy(dtype=float)
    y_sd = float(np.std(y, ddof=1))
    if y_sd <= 0:
        return {"beta": np.nan, "se": np.nan, "low": np.nan, "high": np.nan, "p": np.nan}
    y = (y - np.mean(y)) / y_sd
    design = pd.DataFrame({"case": data["comparison_group"].eq("case").astype(float).to_numpy()}, index=data.index)
    if adjusted:
        design["log10_n_cells"] = np.log10(data["n_cells"].astype(float).to_numpy() + 1)
        design["log10_library_size"] = np.log10(data["library_size"].astype(float).to_numpy() + 1)
    design = sm.add_constant(design, has_constant="add")
    fit = sm.OLS(y, design).fit(cov_type="HC3")
    beta = float(fit.params["case"])
    se = float(fit.bse["case"])
    return {
        "beta": beta,
        "se": se,
        "low": beta - 1.96 * se,
        "high": beta + 1.96 * se,
        "p": float(fit.pvalues["case"]),
    }


def analysis_d_composition(repo: Path, output: Path) -> dict[str, object]:
    specs = [
        (
            "GSE244832",
            repo / "results" / "sensitivity" / "gse244832_donor_program_scores.csv.gz",
            "mash_f2f4_group_vs_normal_sensitivity",
            {"normal": "control", "MASH": "case"},
        ),
        (
            "GSE256398",
            repo / "results" / "phase3" / "gse256398_donor_program_scores.csv.gz",
            "mash_cirrhosis_vs_healthy",
            {"healthy": "control", "mash_cirrhosis": "case"},
        ),
        (
            "GSE256398",
            repo / "results" / "phase3" / "gse256398_donor_program_scores.csv.gz",
            "alcohol_cirrhosis_vs_healthy",
            {"healthy": "control", "alcohol_cirrhosis": "case"},
        ),
        (
            "GSE256398",
            repo / "results" / "phase3" / "gse256398_donor_program_scores.csv.gz",
            "mash_fibrosis_vs_masld_f0",
            {"masld_f0": "control", "mash_fibrosis": "case"},
        ),
    ]
    rows: list[dict[str, object]] = []
    for dataset, path, contrast, mapping in specs:
        scores = pd.read_csv(path)
        scores = scores[scores["eligible_30"].astype(bool)].copy()
        scores["comparison_group"] = scores["disease_group"].map(mapping).fillna("excluded")
        scores = scores[scores["comparison_group"].ne("excluded")]
        for (program_id, method), data in scores.groupby(["program_id", "score_method"]):
            if data.groupby("comparison_group")["donor_id"].nunique().min() < 3:
                continue
            unadjusted = model_group_effect(data, adjusted=False)
            adjusted = model_group_effect(data, adjusted=True)
            ratio = abs(adjusted["beta"] / unadjusted["beta"]) if abs(unadjusted["beta"]) > 1e-12 else np.nan
            rows.append({
                "dataset_id": dataset,
                "contrast": contrast,
                "program_id": program_id,
                "lineage": data["lineage"].iloc[0],
                "score_method": method,
                "n_control": int(data["comparison_group"].eq("control").sum()),
                "n_case": int(data["comparison_group"].eq("case").sum()),
                "unadjusted_standardized_beta": unadjusted["beta"],
                "unadjusted_hc3_low": unadjusted["low"],
                "unadjusted_hc3_high": unadjusted["high"],
                "unadjusted_p": unadjusted["p"],
                "adjusted_standardized_beta": adjusted["beta"],
                "adjusted_hc3_low": adjusted["low"],
                "adjusted_hc3_high": adjusted["high"],
                "adjusted_p": adjusted["p"],
                "sign_retained": bool(np.sign(adjusted["beta"]) == np.sign(unadjusted["beta"])),
                "magnitude_retention_ratio": ratio,
            })
    frame = pd.DataFrame(rows)
    frame["adjusted_fdr_within_context_method"] = frame.groupby(
        ["dataset_id", "contrast", "score_method"], group_keys=False
    )["adjusted_p"].apply(bh_fdr)
    frame.to_csv(output / "composition_adjusted_program_effects.csv", index=False)
    summary_rows = []
    for keys, data in frame.groupby(["dataset_id", "contrast", "program_id"]):
        stable = len(data) == 2 and data["sign_retained"].all() and data["magnitude_retention_ratio"].ge(0.70).all()
        summary_rows.append({
            "dataset_id": keys[0],
            "contrast": keys[1],
            "program_id": keys[2],
            "lineage": data["lineage"].iloc[0],
            "composition_stable": bool(stable),
            "minimum_magnitude_retention": float(data["magnitude_retention_ratio"].min()),
            "both_adjusted_intervals_positive": bool(data["adjusted_hc3_low"].gt(0).all()),
        })
    summary = pd.DataFrame(summary_rows)
    summary.to_csv(output / "composition_stability_summary.csv", index=False)
    return {
        "model_rows": len(frame),
        "context_program_rows": len(summary),
        "composition_stable": int(summary["composition_stable"].sum()),
        "stable_positive_intervals": int((summary["composition_stable"] & summary["both_adjusted_intervals_positive"]).sum()),
    }


def endpoint_class(contrast: str) -> str:
    value = contrast.lower()
    if "mash_vs_alcohol" in value:
        return "etiology_contrast"
    if "cirrhosis_vs_healthy" in value:
        return "cirrhosis_vs_healthy"
    if "advanced_f3f4_vs_f0" in value or "f4_vs_f0" in value:
        return "advanced_vs_f0"
    if "mash_f2f4_group_vs_normal" in value:
        return "mixed_mash_vs_normal"
    if "mash_fibrosis_vs_masld_f0" in value:
        return "fibrosis_vs_f0"
    if "mash_cirrhosis_vs_masld_f0" in value:
        return "cirrhosis_vs_f0"
    if "mixed" in value and "healthy" in value:
        return "mixed_fibrosis_vs_healthy"
    if "cirrhosis_vs_healthy" in value:
        return "cirrhosis_vs_healthy"
    return value


def etiology_class(dataset_id: str, contrast: str) -> str:
    value = contrast.lower()
    if "mash_vs_alcohol" in value:
        return "cross_etiology"
    if "alcohol" in value:
        return "alcohol"
    if "mash" in value or "nash" in value or dataset_id in {"GSE202379", "GSE244832"}:
        return "metabolic"
    return "mixed_or_unresolved"


def two_way_variance(matrix: pd.DataFrame) -> dict[str, float]:
    values = matrix.to_numpy(dtype=float)
    n_programs, n_contexts = values.shape
    grand = float(values.mean())
    program_means = values.mean(axis=1)
    context_means = values.mean(axis=0)
    ss_program = n_contexts * float(np.sum((program_means - grand) ** 2))
    ss_context = n_programs * float(np.sum((context_means - grand) ** 2))
    residual = values - program_means[:, None] - context_means[None, :] + grand
    ss_residual = float(np.sum(residual**2))
    ms_program = ss_program / max(1, n_programs - 1)
    ms_context = ss_context / max(1, n_contexts - 1)
    ms_residual = ss_residual / max(1, (n_programs - 1) * (n_contexts - 1))
    var_program = max(0.0, (ms_program - ms_residual) / n_contexts)
    var_context = max(0.0, (ms_context - ms_residual) / n_programs)
    var_residual = max(0.0, ms_residual)
    total = var_program + var_context + var_residual
    return {
        "n_programs": n_programs,
        "n_contexts": n_contexts,
        "program_variance": var_program,
        "context_variance": var_context,
        "interaction_residual_variance": var_residual,
        "program_fraction": var_program / total if total else np.nan,
        "context_fraction": var_context / total if total else np.nan,
        "interaction_residual_fraction": var_residual / total if total else np.nan,
    }


def descriptor_matrix(contexts: list[str], metadata: pd.DataFrame) -> pd.DataFrame:
    lookup = metadata.set_index("context_id")
    rows = []
    for left, right in combinations(contexts, 2):
        a, b = lookup.loc[left], lookup.loc[right]
        rows.append({
            "context_left": left,
            "context_right": right,
            "same_endpoint": float(a["endpoint_class"] == b["endpoint_class"]),
            "same_etiology": float(a["etiology_class"] == b["etiology_class"]),
            "same_assay": float(a["assay"] == b["assay"]),
            "same_annotation": float(a["annotation_role"] == b["annotation_role"]),
        })
    return pd.DataFrame(rows)


def analysis_e_topology(repo: Path, output: Path) -> dict[str, object]:
    effects = pd.read_csv(repo / "results" / "phase3" / "phase3_cross_cohort_effect_matrix.csv")
    effects = effects[np.isfinite(pd.to_numeric(effects["hedges_g"], errors="coerce"))].copy()
    effects["context_id"] = effects["dataset_id"].astype(str) + "::" + effects["contrast"].astype(str)
    effects["endpoint_class"] = effects["contrast"].map(endpoint_class)
    effects["etiology_class"] = [etiology_class(d, c) for d, c in zip(effects["dataset_id"], effects["contrast"])]
    context_meta = effects.drop_duplicates("context_id")[[
        "context_id", "endpoint_class", "etiology_class", "assay", "annotation_role"
    ]].copy()
    variance_rows: list[dict[str, object]] = []
    topology_rows: list[dict[str, object]] = []
    regression_rows: list[dict[str, object]] = []
    predictors = ["same_endpoint", "same_etiology", "same_assay", "same_annotation"]
    for (lineage, method), data in effects.groupby(["lineage", "score_method"]):
        pivot = data.pivot_table(index="program_id", columns="context_id", values="hedges_g", aggfunc="first")
        complete_contexts = pivot.columns[pivot.notna().all(axis=0)].tolist()
        complete_programs = pivot.index[pivot[complete_contexts].notna().all(axis=1)].tolist() if complete_contexts else []
        if len(complete_contexts) >= 2 and len(complete_programs) >= 3:
            variance_rows.append({
                "lineage": lineage,
                "score_method": method,
                **two_way_variance(pivot.loc[complete_programs, complete_contexts]),
            })
        contexts = sorted(pivot.columns)
        group_pair_rows = []
        for left, right in combinations(contexts, 2):
            paired = pivot[[left, right]].dropna()
            if len(paired) < 3:
                continue
            rho = float(spearmanr(paired[left], paired[right]).statistic)
            meta_left = context_meta[context_meta["context_id"].eq(left)].iloc[0]
            meta_right = context_meta[context_meta["context_id"].eq(right)].iloc[0]
            row = {
                "lineage": lineage,
                "score_method": method,
                "context_left": left,
                "context_right": right,
                "shared_programs": len(paired),
                "spearman_rho": rho,
                "sign_agreement": float((np.sign(paired[left]) == np.sign(paired[right])).mean()),
                "same_endpoint": bool(meta_left["endpoint_class"] == meta_right["endpoint_class"]),
                "same_etiology": bool(meta_left["etiology_class"] == meta_right["etiology_class"]),
                "same_assay": bool(meta_left["assay"] == meta_right["assay"]),
                "same_annotation": bool(meta_left["annotation_role"] == meta_right["annotation_role"]),
            }
            topology_rows.append(row)
            group_pair_rows.append(row)
        pair_frame = pd.DataFrame(group_pair_rows)
        if len(pair_frame) < 8:
            continue
        x = pair_frame[predictors].astype(float).to_numpy()
        varying = np.ptp(x, axis=0) > 0
        active_predictors = [name for name, keep in zip(predictors, varying) if keep]
        x = x[:, varying]
        x_design = np.column_stack([np.ones(len(x)), x])
        y = pair_frame["spearman_rho"].to_numpy(dtype=float)
        observed_beta = np.linalg.lstsq(x_design, y, rcond=None)[0]
        contexts_for_group = sorted(set(pair_frame["context_left"]) | set(pair_frame["context_right"]))
        metadata_group = (
            context_meta[context_meta["context_id"].isin(contexts_for_group)]
            .set_index("context_id")
            .loc[contexts_for_group]
            .reset_index()
        )
        context_index = {context: index for index, context in enumerate(contexts_for_group)}
        left_indices = pair_frame["context_left"].map(context_index).to_numpy(dtype=int)
        right_indices = pair_frame["context_right"].map(context_index).to_numpy(dtype=int)
        null = np.empty((PERMUTATIONS, len(active_predictors)), dtype=float)
        rng = np.random.default_rng(stable_seed("topology", lineage, method))
        attribute_columns = ["endpoint_class", "etiology_class", "assay", "annotation_role"]
        attributes = metadata_group[attribute_columns].astype(str).to_numpy()
        active_indices = [predictors.index(name) for name in active_predictors]
        for iteration in range(PERMUTATIONS):
            order = rng.permutation(len(attributes))
            permuted_attributes = attributes[order]
            all_descriptors = (
                permuted_attributes[left_indices, :] == permuted_attributes[right_indices, :]
            ).astype(float)
            px = all_descriptors[:, active_indices]
            pdesign = np.column_stack([np.ones(len(px)), px])
            null[iteration, :] = np.linalg.lstsq(pdesign, y, rcond=None)[0][1:]
        for index, predictor in enumerate(active_predictors):
            beta = float(observed_beta[index + 1])
            p = float((1 + np.sum(np.abs(null[:, index]) >= abs(beta) - 1e-15)) / (PERMUTATIONS + 1))
            regression_rows.append({
                "lineage": lineage,
                "score_method": method,
                "predictor": predictor,
                "coefficient_on_spearman_rho": beta,
                "permutation_p_two_sided": p,
                "pairs": len(pair_frame),
            })
    variance = pd.DataFrame(variance_rows)
    topology = pd.DataFrame(topology_rows)
    regression = pd.DataFrame(regression_rows)
    if not regression.empty:
        regression["fdr_across_all_coefficients"] = bh_fdr(regression["permutation_p_two_sided"])
    variance.to_csv(output / "program_context_variance_components.csv", index=False)
    topology.to_csv(output / "context_transport_topology.csv", index=False)
    regression.to_csv(output / "context_descriptor_permutation_regression.csv", index=False)
    return {
        "variance_rows": len(variance),
        "topology_pairs": len(topology),
        "descriptor_rows": len(regression),
        "descriptor_fdr_hits": int(regression.get("fdr_across_all_coefficients", pd.Series(dtype=float)).lt(0.05).sum()),
    }


def main() -> None:
    repo = Path(__file__).resolve().parents[1]
    output = repo / "results" / "phase4"
    output.mkdir(parents=True, exist_ok=True)
    summary = {
        "analysis_freeze": "reports/phase4_biological_structure_plan_2026-08-31.md",
        "seed": SEED,
        "permutations": PERMUTATIONS,
        "random_modules": RANDOM_MODULES,
    }
    summary["member_gene_coherence"] = analysis_a_member_gene_coherence(repo, output)
    summary["cross_lineage_coupling"] = analysis_b_cross_lineage_coupling(repo, output)
    summary["etiology_geometry"] = analysis_c_etiology_geometry(repo, output)
    summary["composition"] = analysis_d_composition(repo, output)
    summary["topology"] = analysis_e_topology(repo, output)
    log_dir = repo / "results" / "logs"
    log_dir.mkdir(exist_ok=True)
    (log_dir / "phase4_biological_structure_run.json").write_text(
        json.dumps(summary, indent=2, sort_keys=True), encoding="utf-8"
    )
    print(json.dumps(summary, indent=2, sort_keys=True))


if __name__ == "__main__":
    main()
