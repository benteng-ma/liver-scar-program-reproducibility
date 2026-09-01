from __future__ import annotations

import json
import math
from itertools import combinations
from pathlib import Path

import numpy as np
import pandas as pd
import statsmodels.api as sm
from scipy.io import mmread
from scipy.stats import rankdata, spearmanr

from analyze_gse202379_programs import singscore_up, stable_seed
from random_module_benchmark_gse202379 import score_modules, vectorized_hedges_g


DATASETS = {
    "GSE202379": {
        "interim": "GSE202379",
        "control": {"Healthy control"},
        "case": {"NASH with cirrhosis"},
        "group_column": "Disease.status",
        "contrast": "clinical_cirrhosis_vs_healthy",
        "random_prefix": "gse202379",
        "broad_scores": "results/primary/gse202379_donor_program_scores.csv.gz",
    },
    "GSE244832": {
        "interim": "GSE244832",
        "control": {"normal"},
        "case": {"MASH"},
        "group_column": "disease_group",
        "contrast": "mash_f2f4_group_vs_normal_sensitivity",
        "random_prefix": "gse244832",
        "broad_scores": "results/sensitivity/gse244832_donor_program_scores.csv.gz",
    },
    "GSE256398_human": {
        "interim": "GSE256398",
        "control": {"healthy"},
        "case": {"mash_cirrhosis"},
        "group_column": "disease_group",
        "contrast": "mash_cirrhosis_vs_healthy",
        "random_prefix": "gse256398",
        "broad_scores": "results/phase3/gse256398_donor_program_scores.csv.gz",
    },
}
PRIMARY_PROGRAMS = {"RAM2019_ENDO_2", "RAM2019_ENDO_6_SAENDO1"}
GATES = ((30, "primary"), (20, "sensitivity"))
BOOTSTRAPS = 10_000


def effect_hc3(control: np.ndarray, case: np.ndarray) -> dict[str, float]:
    control = np.asarray(control, dtype=float)
    case = np.asarray(case, dtype=float)
    n_control, n_case = len(control), len(case)
    df = n_control + n_case - 2
    var_control, var_case = control.var(ddof=1), case.var(ddof=1)
    pooled_variance = ((n_control - 1) * var_control + (n_case - 1) * var_case) / df
    mean_difference = float(case.mean() - control.mean())
    if pooled_variance <= 0 or not np.isfinite(pooled_variance):
        return {
            "n_control": n_control,
            "n_case": n_case,
            "mean_control": float(control.mean()),
            "mean_case": float(case.mean()),
            "mean_difference": mean_difference,
            "pooled_sd": math.nan,
            "hedges_g": math.nan,
            "robust_se_g_hc3": math.nan,
            "robust_ci95_low": math.nan,
            "robust_ci95_high": math.nan,
            "effect_status": "ZERO_POOLED_VARIANCE",
        }
    pooled_sd = math.sqrt(pooled_variance)
    correction = 1 - 3 / (4 * df - 1)
    g = correction * mean_difference / pooled_sd
    values = np.concatenate([control, case])
    labels = np.concatenate([np.zeros(n_control), np.ones(n_case)])
    fit = sm.OLS(values, sm.add_constant(labels)).fit(cov_type="HC3")
    se_g = correction * float(fit.bse[1]) / pooled_sd
    return {
        "n_control": n_control,
        "n_case": n_case,
        "mean_control": float(control.mean()),
        "mean_case": float(case.mean()),
        "mean_difference": mean_difference,
        "pooled_sd": pooled_sd,
        "hedges_g": g,
        "robust_se_g_hc3": se_g,
        "robust_ci95_low": g - 1.96 * se_g,
        "robust_ci95_high": g + 1.96 * se_g,
        "effect_status": "ESTIMATED",
    }


def bh(values: pd.Series) -> pd.Series:
    valid = values.notna()
    result = pd.Series(np.nan, index=values.index)
    if not valid.any():
        return result
    x = values.loc[valid].astype(float).to_numpy()
    order = np.argsort(x)
    adjusted = x[order] * len(x) / np.arange(1, len(x) + 1)
    adjusted = np.minimum.accumulate(adjusted[::-1])[::-1]
    inverse = np.empty_like(order)
    inverse[order] = np.arange(len(order))
    result.loc[valid] = np.minimum(adjusted[inverse], 1.0)
    return result


def marker_profiles(
    genes: pd.Series,
    counts,
    manifest: pd.DataFrame,
    excluded_genes: set[str],
) -> tuple[pd.DataFrame, dict[tuple[str, str], pd.Series]]:
    rows: list[dict[str, object]] = []
    profiles: dict[tuple[str, str], pd.Series] = {}
    for lineage, lineage_manifest in manifest.groupby("harmonized_lineage", sort=False):
        state_names = sorted(lineage_manifest["source_state"].unique())
        state_columns = []
        for state in state_names:
            columns = lineage_manifest.index[lineage_manifest["source_state"].eq(state)].to_numpy()
            state_columns.append(np.asarray(counts[:, columns].sum(axis=1)).ravel())
        pooled = np.column_stack(state_columns)
        library = pooled.sum(axis=0)
        log_cpm = np.log2(pooled / library * 1_000_000 + 1)
        for index, state in enumerate(state_names):
            if len(state_names) == 1:
                effect = log_cpm[:, index]
            else:
                effect = log_cpm[:, index] - np.delete(log_cpm, index, axis=1).mean(axis=1)
            series = pd.Series(effect, index=genes.str.upper())
            series = series.groupby(level=0, sort=False).mean()
            series = series[~series.index.isin(excluded_genes)]
            profiles[(lineage, state)] = series
            top = series[series.gt(0)].sort_values(ascending=False).head(50)
            for rank, (gene, value) in enumerate(top.items(), start=1):
                rows.append(
                    {
                        "lineage": lineage,
                        "source_state": state,
                        "marker_rank": rank,
                        "gene_symbol": gene,
                        "marker_effect": float(value),
                    }
                )
    return pd.DataFrame(rows), profiles


def comparison_labels(manifest: pd.DataFrame, cfg: dict[str, object]) -> pd.Series:
    column = str(cfg["group_column"])
    control = set(cfg["control"])
    case = set(cfg["case"])
    return manifest[column].astype(str).map(
        {
            **{value: "control" for value in control},
            **{value: "case" for value in case},
        }
    ).fillna("excluded")


def decomposition(
    scores: pd.DataFrame,
    manifest: pd.DataFrame,
    labels: pd.Series,
    lineage: str,
    program_id: str,
    method: str,
    seed: int,
) -> tuple[dict[str, object], list[dict[str, object]]]:
    scoped_manifest = manifest[manifest["harmonized_lineage"].eq(lineage)].copy()
    scoped_manifest["comparison_group"] = labels.loc[scoped_manifest.index]
    scoped_manifest = scoped_manifest[scoped_manifest["comparison_group"].ne("excluded")]
    donors = sorted(scoped_manifest["donor_id"].astype(str).unique())
    states = sorted(scoped_manifest["source_state"].unique())
    donor_index = {donor: index for index, donor in enumerate(donors)}
    state_index = {state: index for index, state in enumerate(states)}
    p = np.zeros((len(donors), len(states)), dtype=float)
    mu = np.full((len(donors), len(states)), np.nan, dtype=float)
    groups = np.empty(len(donors), dtype=object)
    for donor in donors:
        donor_rows = scoped_manifest[scoped_manifest["donor_id"].astype(str).eq(donor)]
        total = donor_rows["n_cells"].sum()
        labels_for_donor = donor_rows["comparison_group"].unique()
        if len(labels_for_donor) != 1:
            raise RuntimeError("comparison group differs across donor states")
        groups[donor_index[donor]] = labels_for_donor[0]
        for row in donor_rows.itertuples():
            p[donor_index[donor], state_index[row.source_state]] = row.n_cells / total
    subset = scores[
        scores["harmonized_lineage"].eq(lineage)
        & scores["program_id"].eq(program_id)
        & scores["score_method"].eq(method)
    ]
    for row in subset.itertuples():
        donor = str(row.donor_id)
        if donor in donor_index and row.source_state in state_index:
            mu[donor_index[donor], state_index[row.source_state]] = row.score
    pooled_mu = np.nanmean(mu, axis=0)
    pooled_mu = np.where(np.isfinite(pooled_mu), pooled_mu, 0.0)

    def quiet_nanmean(values: np.ndarray) -> np.ndarray:
        valid = np.isfinite(values)
        count = valid.sum(axis=0)
        total = np.where(valid, values, 0.0).sum(axis=0)
        result = np.full(values.shape[1], np.nan, dtype=float)
        np.divide(total, count, out=result, where=count > 0)
        return result

    def calculate(indices_control: np.ndarray, indices_case: np.ndarray):
        p0, p1 = p[indices_control].mean(axis=0), p[indices_case].mean(axis=0)
        mu0 = quiet_nanmean(mu[indices_control])
        mu1 = quiet_nanmean(mu[indices_case])
        mu0 = np.where(np.isfinite(mu0), mu0, pooled_mu)
        mu1 = np.where(np.isfinite(mu1), mu1, pooled_mu)
        abundance = (p1 - p0) * (mu1 + mu0) / 2
        intensity = (mu1 - mu0) * (p1 + p0) / 2
        total = float(np.sum(p1 * mu1) - np.sum(p0 * mu0))
        return total, float(abundance.sum()), float(intensity.sum()), p0, p1, mu0, mu1, abundance, intensity

    control_indices = np.flatnonzero(groups == "control")
    case_indices = np.flatnonzero(groups == "case")
    observed = calculate(control_indices, case_indices)
    rng = np.random.default_rng(seed)
    boot = np.empty((BOOTSTRAPS, 3), dtype=float)
    for iteration in range(BOOTSTRAPS):
        sampled_control = rng.choice(control_indices, len(control_indices), replace=True)
        sampled_case = rng.choice(case_indices, len(case_indices), replace=True)
        boot[iteration, :3] = calculate(sampled_control, sampled_case)[:3]
    total, abundance, intensity, p0, p1, mu0, mu1, abundance_state, intensity_state = observed
    summary = {
        "lineage": lineage,
        "program_id": program_id,
        "score_method": method,
        "n_control_donors": len(control_indices),
        "n_case_donors": len(case_indices),
        "states": len(states),
        "reconstructed_total_difference": total,
        "abundance_component": abundance,
        "intensity_component": intensity,
        "identity_error": total - abundance - intensity,
        "abundance_share": abundance / total if total != 0 else math.nan,
        "intensity_share": intensity / total if total != 0 else math.nan,
        "total_bootstrap_low": float(np.quantile(boot[:, 0], 0.025)),
        "total_bootstrap_high": float(np.quantile(boot[:, 0], 0.975)),
        "abundance_bootstrap_low": float(np.quantile(boot[:, 1], 0.025)),
        "abundance_bootstrap_high": float(np.quantile(boot[:, 1], 0.975)),
        "intensity_bootstrap_low": float(np.quantile(boot[:, 2], 0.025)),
        "intensity_bootstrap_high": float(np.quantile(boot[:, 2], 0.975)),
        "bootstrap_iterations": BOOTSTRAPS,
    }
    state_rows = []
    for index, state in enumerate(states):
        state_rows.append(
            {
                "lineage": lineage,
                "program_id": program_id,
                "score_method": method,
                "source_state": state,
                "mean_fraction_control": p0[index],
                "mean_fraction_case": p1[index],
                "mean_intensity_control": mu0[index],
                "mean_intensity_case": mu1[index],
                "abundance_component": abundance_state[index],
                "intensity_component": intensity_state[index],
                "combined_component": abundance_state[index] + intensity_state[index],
            }
        )
    return summary, state_rows


def main() -> None:
    repo = Path(__file__).resolve().parents[1]
    output = repo / "results" / "phase5"
    output.mkdir(parents=True, exist_ok=True)
    programs = pd.read_csv(repo / "literature" / "program_inventory.csv")
    program_genes = set(programs["gene_symbol"].astype(str).str.upper())

    all_scores = []
    all_markers = []
    all_intensity_effects = []
    all_abundance_effects = []
    all_decomposition = []
    all_contributions = []
    all_random = []
    profile_store: dict[str, dict[tuple[str, str], pd.Series]] = {}

    for dataset_id, cfg in DATASETS.items():
        interim = repo / "data" / "interim" / str(cfg["interim"])
        genes = pd.read_csv(interim / "phase5_state_genes.csv")["gene"].astype(str)
        gene_to_index = {gene.upper(): index for index, gene in enumerate(genes)}
        manifest = pd.read_csv(interim / "phase5_donor_state_manifest.csv")
        manifest["donor_id"] = manifest["donor_id"].astype(str)
        counts = mmread(interim / "phase5_donor_state_raw_counts.mtx").tocsc()
        if counts.shape != (len(genes), len(manifest)):
            raise RuntimeError(f"{dataset_id}: state matrix mismatch")
        labels = comparison_labels(manifest, cfg)

        marker_table, profiles = marker_profiles(genes, counts, manifest, program_genes)
        marker_table.insert(0, "dataset_id", dataset_id)
        all_markers.append(marker_table)
        profile_store[dataset_id] = profiles

        library = np.asarray(counts.sum(axis=0)).ravel()
        log_cpm = np.log2(counts.toarray() / library * 1_000_000 + 1)
        ranks = rankdata(log_cpm, axis=0, method="average")
        dataset_score_rows: list[dict[str, object]] = []
        score_cache: dict[tuple[str, str], np.ndarray] = {}
        lineage_arrays: dict[str, dict[str, object]] = {}

        for lineage, lineage_manifest in manifest.groupby("harmonized_lineage", sort=False):
            columns = lineage_manifest.index.to_numpy()
            reference = lineage_manifest.index[lineage_manifest["n_cells"].ge(20)].to_numpy()
            if len(reference) < 2:
                reference = columns
            mean = log_cpm[:, reference].mean(axis=1)
            sd = log_cpm[:, reference].std(axis=1, ddof=1)
            invariant = ~np.isfinite(sd) | (sd == 0)
            safe_sd = sd.copy()
            safe_sd[invariant] = 1
            standardized = (log_cpm[:, columns] - mean[:, None]) / safe_sd[:, None]
            standardized[invariant, :] = 0
            lineage_ranks = ranks[:, columns]
            global_to_local = {int(value): index for index, value in enumerate(columns)}
            lineage_arrays[lineage] = {
                "columns": columns,
                "global_to_local": global_to_local,
                "standardized": standardized,
                "ranks": lineage_ranks,
            }
            for program_id, rows in programs[programs["cell_lineage"].eq(lineage)].groupby("program_id"):
                measured = sorted(set(rows["gene_symbol"].astype(str).str.upper()) & set(gene_to_index))
                if len(measured) / len(rows) < 0.4:
                    continue
                indices = np.array([gene_to_index[value] for value in measured], dtype=int)
                sing = singscore_up(lineage_ranks, indices)
                zmean = standardized[indices, :].mean(axis=0)
                for method, values in (("singscore", sing), ("standardized_mean", zmean)):
                    full = np.full(len(manifest), np.nan)
                    full[columns] = values
                    score_cache[(program_id, method)] = full
                    for local, global_index in enumerate(columns):
                        dataset_score_rows.append(
                            {
                                "dataset_id": dataset_id,
                                "contrast": cfg["contrast"],
                                "state_group_id": manifest.loc[global_index, "state_group_id"],
                                "donor_id": str(manifest.loc[global_index, "donor_id"]),
                                "harmonized_lineage": lineage,
                                "source_state": manifest.loc[global_index, "source_state"],
                                "n_cells": int(manifest.loc[global_index, "n_cells"]),
                                "comparison_group": labels.loc[global_index],
                                "program_id": program_id,
                                "score_method": method,
                                "score": float(values[local]),
                                "measured_program_genes": len(measured),
                                "program_coverage": len(measured) / len(rows),
                            }
                        )

        dataset_scores = pd.DataFrame(dataset_score_rows)
        all_scores.append(dataset_scores)

        # State abundance effects are program-independent.
        for lineage, lineage_manifest in manifest.groupby("harmonized_lineage", sort=False):
            scoped = lineage_manifest.assign(comparison_group=labels.loc[lineage_manifest.index])
            scoped = scoped[scoped["comparison_group"].ne("excluded")]
            donors = sorted(scoped["donor_id"].astype(str).unique())
            donor_groups = (
                scoped.groupby("donor_id")["comparison_group"].first().astype(str).to_dict()
            )
            donor_totals = scoped.groupby("donor_id")["n_cells"].sum().to_dict()
            for state in sorted(scoped["source_state"].unique()):
                observed = scoped[scoped["source_state"].eq(state)].set_index("donor_id")["n_cells"].to_dict()
                fractions = pd.Series(
                    {donor: observed.get(donor, 0) / donor_totals[donor] for donor in donors}
                )
                transformed = np.arcsin(np.sqrt(fractions.to_numpy(float)))
                groups = np.array([donor_groups[donor] for donor in fractions.index])
                for gate, tier in GATES:
                    # Broad-lineage totals define eligibility; states can be structurally zero.
                    eligible = np.array([donor_totals[donor] >= gate for donor in fractions.index])
                    control = transformed[eligible & (groups == "control")]
                    case = transformed[eligible & (groups == "case")]
                    if len(control) < 3 or len(case) < 3:
                        continue
                    stats = effect_hc3(control, case)
                    all_abundance_effects.append(
                        {
                            "dataset_id": dataset_id,
                            "contrast": cfg["contrast"],
                            "lineage": lineage,
                            "source_state": state,
                            "cell_gate": gate,
                            "analysis_tier": tier,
                            "mean_fraction_control": float(fractions.to_numpy()[eligible & (groups == "control")].mean()),
                            "mean_fraction_case": float(fractions.to_numpy()[eligible & (groups == "case")].mean()),
                            **stats,
                        }
                    )

        # State-specific program intensity effects.
        for (lineage, state, program_id, method), data in dataset_scores.groupby(
            ["harmonized_lineage", "source_state", "program_id", "score_method"]
        ):
            for gate, tier in GATES:
                selected = data[data["n_cells"].ge(gate) & data["comparison_group"].isin(["control", "case"])]
                control = selected.loc[selected["comparison_group"].eq("control"), "score"].to_numpy(float)
                case = selected.loc[selected["comparison_group"].eq("case"), "score"].to_numpy(float)
                if len(control) < 3 or len(case) < 3:
                    continue
                stats = effect_hc3(control, case)
                all_intensity_effects.append(
                    {
                        "dataset_id": dataset_id,
                        "contrast": cfg["contrast"],
                        "lineage": lineage,
                        "source_state": state,
                        "program_id": program_id,
                        "score_method": method,
                        "cell_gate": gate,
                        "analysis_tier": tier,
                        **stats,
                    }
                )

        # Decomposition for every program and method.
        for (program_id, method), data in dataset_scores.groupby(["program_id", "score_method"]):
            lineage = str(data["harmonized_lineage"].iloc[0])
            summary, contributions = decomposition(
                dataset_scores,
                manifest,
                labels,
                lineage,
                program_id,
                method,
                stable_seed("phase5_decomposition", dataset_id, program_id, method),
            )
            summary.update({"dataset_id": dataset_id, "contrast": cfg["contrast"]})
            for row in contributions:
                row.update({"dataset_id": dataset_id, "contrast": cfg["contrast"]})
            all_decomposition.append(summary)
            all_contributions.extend(contributions)

        # Frozen random modules for primary-program eligible state effects.
        membership_path = interim / "random_modules" / "matched_random_module_membership.csv.gz"
        membership = pd.read_csv(membership_path)
        primary_effects_dataset = [
            row for row in all_intensity_effects
            if row["dataset_id"] == dataset_id
            and row["analysis_tier"] == "primary"
            and row["program_id"] in PRIMARY_PROGRAMS
        ]
        for program_id in sorted(PRIMARY_PROGRAMS):
            program_rows = programs[programs["program_id"].eq(program_id)]
            if program_rows.empty:
                continue
            lineage = str(program_rows["cell_lineage"].iloc[0])
            if lineage not in lineage_arrays:
                continue
            members = membership[
                membership["lineage"].eq(lineage) & membership["program_id"].eq(program_id)
            ]
            modules = []
            for module_id in sorted(members["module_id"].unique()):
                genes_for_module = members.loc[members["module_id"].eq(module_id), "random_gene"].astype(str)
                modules.append([gene_to_index[value.upper()] for value in genes_for_module])
            if len(modules) != 1000:
                raise RuntimeError(f"{dataset_id} {program_id}: frozen module count != 1000")
            module_array = np.asarray(modules, dtype=int)
            arrays = lineage_arrays[lineage]
            random_sing, random_zmean = score_modules(
                module_array, arrays["ranks"], arrays["standardized"]
            )
            global_to_local = arrays["global_to_local"]
            for effect in primary_effects_dataset:
                if effect["program_id"] != program_id:
                    continue
                selected = manifest[
                    manifest["harmonized_lineage"].eq(lineage)
                    & manifest["source_state"].eq(effect["source_state"])
                    & manifest["n_cells"].ge(effect["cell_gate"])
                    & labels.isin(["control", "case"])
                ]
                selected_global = selected.index.to_numpy()
                selected_local = np.array([global_to_local[int(value)] for value in selected_global], dtype=int)
                selected_labels = labels.loc[selected_global]
                control_local = selected_local[selected_labels.eq("control").to_numpy()]
                case_local = selected_local[selected_labels.eq("case").to_numpy()]
                module_scores = random_sing if effect["score_method"] == "singscore" else random_zmean
                random_g = vectorized_hedges_g(module_scores, control_local, case_local)
                valid = random_g[np.isfinite(random_g)]
                real_g = float(effect["hedges_g"])
                q95 = float(np.quantile(valid, 0.95))
                all_random.append(
                    {
                        "dataset_id": dataset_id,
                        "contrast": cfg["contrast"],
                        "lineage": lineage,
                        "source_state": effect["source_state"],
                        "program_id": program_id,
                        "score_method": effect["score_method"],
                        "cell_gate": effect["cell_gate"],
                        "real_hedges_g": real_g,
                        "random_g_median": float(np.median(valid)),
                        "random_g_95th_percentile": q95,
                        "real_effect_percentile": float((valid <= real_g).mean()),
                        "above_random_95th_percentile": real_g > q95,
                    }
                )

    score_table = pd.concat(all_scores, ignore_index=True)
    marker_table = pd.concat(all_markers, ignore_index=True)
    intensity = pd.DataFrame(all_intensity_effects)
    abundance = pd.DataFrame(all_abundance_effects)
    decomposition_table = pd.DataFrame(all_decomposition)
    contribution_table = pd.DataFrame(all_contributions)
    random_table = pd.DataFrame(all_random)

    # Cross-cohort state matching uses only non-program marker profiles.
    match_rows = []
    for first, second in combinations(DATASETS, 2):
        for (lineage_a, state_a), profile_a in profile_store[first].items():
            for (lineage_b, state_b), profile_b in profile_store[second].items():
                if lineage_a != lineage_b:
                    continue
                common = profile_a.index.intersection(profile_b.index)
                x, y = profile_a.loc[common].to_numpy(), profile_b.loc[common].to_numpy()
                variable = np.isfinite(x) & np.isfinite(y)
                rho = float(spearmanr(x[variable], y[variable]).statistic)
                top_a = set(
                    marker_table[
                        marker_table["dataset_id"].eq(first)
                        & marker_table["lineage"].eq(lineage_a)
                        & marker_table["source_state"].eq(state_a)
                    ]["gene_symbol"]
                )
                top_b = set(
                    marker_table[
                        marker_table["dataset_id"].eq(second)
                        & marker_table["lineage"].eq(lineage_b)
                        & marker_table["source_state"].eq(state_b)
                    ]["gene_symbol"]
                )
                shared = sorted(top_a & top_b)
                union = top_a | top_b
                jaccard = len(shared) / len(union) if union else 0.0
                match_rows.append(
                    {
                        "dataset_a": first,
                        "state_a": state_a,
                        "dataset_b": second,
                        "state_b": state_b,
                        "lineage": lineage_a,
                        "common_nonprogram_genes": len(common),
                        "marker_effect_spearman_rho": rho,
                        "top50_shared_markers": len(shared),
                        "top50_marker_jaccard": jaccard,
                        "shared_markers": ";".join(shared),
                        "supported_match": bool(rho > 0 and jaccard >= 0.10 and len(shared) >= 3),
                    }
                )
    matches = pd.DataFrame(match_rows)

    # Frozen post-lock state support rule.
    primary_intensity = intensity[
        intensity["analysis_tier"].eq("primary")
        & intensity["program_id"].isin(PRIMARY_PROGRAMS)
    ]
    positive_states: dict[str, set[tuple[str, str]]] = {}
    for program_id, data in primary_intensity.groupby("program_id"):
        qualifying = set()
        for (dataset_id, state), state_data in data.groupby(["dataset_id", "source_state"]):
            if len(state_data) == 2 and state_data["hedges_g"].gt(0).all():
                qualifying.add((dataset_id, state))
        positive_states[program_id] = qualifying
    support_rows = []
    for program_id in sorted(PRIMARY_PROGRAMS):
        qualifying = positive_states.get(program_id, set())
        supporting_pairs = []
        for row in matches[matches["lineage"].eq("endothelial") & matches["supported_match"]].itertuples():
            if (row.dataset_a, row.state_a) in qualifying and (row.dataset_b, row.state_b) in qualifying:
                supporting_pairs.append(f"{row.dataset_a}:{row.state_a}<->{row.dataset_b}:{row.state_b}")
        support_rows.append(
            {
                "program_id": program_id,
                "qualifying_dataset_state_count": len(qualifying),
                "qualifying_datasets": ";".join(sorted({value[0] for value in qualifying})),
                "supporting_state_pairs": ";".join(supporting_pairs),
                "post_lock_label": (
                    "STATE_SUPPORTED_POST_LOCK"
                    if len({value[0] for value in qualifying}) >= 2 and supporting_pairs
                    else "NOT_STATE_SUPPORTED_POST_LOCK"
                ),
            }
        )
    support = pd.DataFrame(support_rows)

    score_table.to_csv(output / "state_program_scores.csv.gz", index=False, compression="gzip")
    marker_table.to_csv(output / "state_top50_nonprogram_markers.csv.gz", index=False, compression="gzip")
    abundance.to_csv(output / "state_abundance_effects.csv", index=False)
    intensity.to_csv(output / "state_intensity_effects.csv", index=False)
    decomposition_table.to_csv(output / "state_abundance_intensity_decomposition.csv", index=False)
    contribution_table.to_csv(output / "state_component_contributions.csv", index=False)
    matches.to_csv(output / "cross_cohort_state_matches.csv", index=False)
    random_table.to_csv(output / "primary_state_random_benchmark.csv", index=False)
    support.to_csv(output / "primary_state_support_summary.csv", index=False)

    summary = {
        "datasets": list(DATASETS),
        "state_score_rows": len(score_table),
        "state_abundance_effect_rows": len(abundance),
        "state_intensity_effect_rows": len(intensity),
        "decomposition_rows": len(decomposition_table),
        "state_match_rows": len(matches),
        "supported_state_matches": int(matches["supported_match"].sum()),
        "primary_state_support": dict(zip(support["program_id"], support["post_lock_label"])),
        "bootstrap_iterations": BOOTSTRAPS,
    }
    (output / "state_decomposition_run_summary.json").write_text(
        json.dumps(summary, indent=2) + "\n", encoding="utf-8"
    )
    print(support.to_string(index=False))
    print(json.dumps(summary, indent=2))


if __name__ == "__main__":
    main()
