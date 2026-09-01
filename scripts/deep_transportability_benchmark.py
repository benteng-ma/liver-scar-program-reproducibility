from __future__ import annotations

import gzip
import hashlib
import json
import math
from dataclasses import dataclass
from itertools import combinations
from pathlib import Path

import numpy as np
import pandas as pd
from scipy.io import mmread
from scipy.stats import rankdata, spearmanr

from analyze_gse202379_programs import effect_statistics, singscore_up
from random_module_benchmark_gse202379 import (
    N_RANDOM_MODULES,
    deciles,
    matched_modules,
    score_modules,
    vectorized_hedges_g,
)


SEED = 20260831
N_RESAMPLES = 10_000

LINEAGES = (
    "macrophage_monocyte",
    "endothelial",
    "mesenchymal_hsc_myofibroblast",
)

COHORTS = {
    "GSE202379": {
        "interim": "GSE202379",
        "dataset_id": "GSE202379",
        "assay": "snRNA-seq",
        "annotation": "author",
        "identity_evidence": "non-circular author labels",
        "contrast": "advanced_f3f4_vs_f0_non_end_stage",
    },
    "GSE244832": {
        "interim": "GSE244832",
        "dataset_id": "GSE244832",
        "assay": "snRNA-seq",
        "annotation": "reconstructed",
        "identity_evidence": "partly circular cluster-mapping QC",
        "contrast": "mash_f2f4_group_vs_normal_sensitivity",
    },
    "GSE290642": {
        "interim": "GSE290642",
        "dataset_id": "GSE290642_human",
        "assay": "scRNA-seq",
        "annotation": "reconstructed",
        "identity_evidence": "circular reconstructed-label QC",
        "contrast": "f4_vs_f0_reconstructed_label_sensitivity",
    },
    "Watson6": {
        "interim": "GSE210077_Watson6",
        "dataset_id": "GSE210077_Watson6",
        "assay": "snRNA-seq",
        "annotation": "author",
        "identity_evidence": "non-circular author labels",
        "contrast": "mixed_f2f4_fibrosis_vs_healthy_sensitivity",
    },
    "GSE181483": {
        "interim": "GSE181483",
        "dataset_id": "GSE181483_human",
        "assay": "scRNA-seq",
        "annotation": "reconstructed",
        "identity_evidence": "circular reconstructed-label QC",
        "contrast": "cirrhosis_vs_healthy_directional",
    },
}


@dataclass
class CohortData:
    key: str
    dataset_id: str
    assay: str
    annotation: str
    identity_evidence: str
    contrast: str
    genes: pd.Series
    manifest: pd.DataFrame
    counts: object
    log_cpm: np.ndarray
    ranks: np.ndarray
    gene_to_index: dict[str, int]


def stable_seed(*parts: object) -> int:
    payload = "|".join(str(part) for part in parts).encode("utf-8")
    offset = int.from_bytes(hashlib.sha256(payload).digest()[:8], "little")
    return (SEED + offset) % (2**32 - 1)


def sha256(path: Path) -> str:
    digest = hashlib.sha256()
    with path.open("rb") as handle:
        for block in iter(lambda: handle.read(8 * 1024 * 1024), b""):
            digest.update(block)
    return digest.hexdigest().upper()


def bh_adjust(values: pd.Series) -> pd.Series:
    result = pd.Series(np.nan, index=values.index, dtype=float)
    valid = values.dropna().astype(float)
    if valid.empty:
        return result
    order = np.argsort(valid.to_numpy())
    ranked = valid.to_numpy()[order]
    adjusted = ranked * len(ranked) / np.arange(1, len(ranked) + 1)
    adjusted = np.minimum.accumulate(adjusted[::-1])[::-1]
    adjusted = np.clip(adjusted, 0, 1)
    result.loc[valid.index[order]] = adjusted
    return result


def load_cohort(repo: Path, key: str) -> CohortData:
    cfg = COHORTS[key]
    interim = repo / "data" / "interim" / str(cfg["interim"])
    genes = pd.read_csv(interim / "genes.csv")["gene"].astype(str)
    manifest = pd.read_csv(interim / "donor_lineage_manifest.csv")
    manifest["eligible_30"] = manifest["n_cells"].ge(30)
    manifest["eligible_20"] = manifest["n_cells"].ge(20)
    counts = mmread(interim / "donor_lineage_raw_counts.mtx").tocsc()
    if key == "GSE290642":
        shared = set(pd.read_csv(interim / "shared_genes.csv")["gene"].astype(str).str.upper())
        keep = np.array([gene.upper() in shared for gene in genes], dtype=bool)
        genes = genes.loc[keep].reset_index(drop=True)
        counts = counts[keep, :]
    if counts.shape != (len(genes), len(manifest)):
        raise RuntimeError(f"{key}: count matrix does not match gene and manifest dimensions")
    library = np.asarray(counts.sum(axis=0)).ravel()
    if (library <= 0).any():
        raise RuntimeError(f"{key}: zero-library donor-lineage pseudobulk")
    log_cpm = np.log2(counts.toarray() / library * 1_000_000 + 1)
    ranks = rankdata(log_cpm, axis=0, method="average")
    gene_to_index = {gene.upper(): index for index, gene in enumerate(genes)}
    return CohortData(
        key=key,
        dataset_id=str(cfg["dataset_id"]),
        assay=str(cfg["assay"]),
        annotation=str(cfg["annotation"]),
        identity_evidence=str(cfg["identity_evidence"]),
        contrast=str(cfg["contrast"]),
        genes=genes,
        manifest=manifest,
        counts=counts,
        log_cpm=log_cpm,
        ranks=ranks,
        gene_to_index=gene_to_index,
    )


def representative_groups(data: CohortData) -> pd.Series:
    manifest = data.manifest
    group = pd.Series("excluded", index=manifest.index, dtype=object)
    if data.key == "GSE202379":
        stage = pd.to_numeric(manifest["Fibrosis.score..F0.4."], errors="coerce")
        non_end = ~manifest["Disease.status"].eq("end stage")
        group.loc[non_end & stage.eq(0)] = "control"
        group.loc[non_end & stage.ge(3) & stage.le(4)] = "case"
    elif data.key == "GSE244832":
        group.loc[manifest["disease_group"].eq("normal")] = "control"
        group.loc[manifest["disease_group"].eq("MASH")] = "case"
    elif data.key == "GSE290642":
        group.loc[manifest["fibrosis_stage"].eq("F0")] = "control"
        group.loc[manifest["fibrosis_stage"].eq("F4")] = "case"
    elif data.key == "Watson6":
        group.loc[manifest["disease_group"].eq("healthy")] = "control"
        group.loc[manifest["disease_group"].eq("fibrosis")] = "case"
    elif data.key == "GSE181483":
        group.loc[manifest["disease_group"].eq("healthy")] = "control"
        group.loc[manifest["disease_group"].eq("cirrhosis")] = "case"
    return group


def measured_module(data: CohortData, rows: pd.DataFrame) -> tuple[np.ndarray, np.ndarray, list[str], float]:
    frozen = rows.drop_duplicates("gene_symbol").copy()
    measured = frozen[frozen["gene_symbol"].astype(str).str.upper().isin(data.gene_to_index)].copy()
    genes = measured["gene_symbol"].astype(str).str.upper().tolist()
    indices = np.array([data.gene_to_index[gene] for gene in genes], dtype=int)
    signs = measured["direction"].map({"UP": 1.0, "DOWN": -1.0, "UP_IN_STATE": 1.0}).to_numpy(dtype=float)
    return indices, signs, genes, len(genes) / max(frozen["gene_symbol"].nunique(), 1)


def signed_singscore(ranks: np.ndarray, indices: np.ndarray, signs: np.ndarray) -> np.ndarray:
    pieces: list[np.ndarray] = []
    weights: list[int] = []
    up = indices[signs > 0]
    down = indices[signs < 0]
    if len(up):
        pieces.append(singscore_up(ranks, up))
        weights.append(len(up))
    if len(down):
        pieces.append(-singscore_up(ranks, down))
        weights.append(len(down))
    if not pieces:
        raise ValueError("module has no measured genes")
    return np.average(np.vstack(pieces), axis=0, weights=np.asarray(weights))


def module_scores(
    data: CohortData,
    rows: pd.DataFrame,
    selected_columns: np.ndarray,
    reference_columns: np.ndarray,
) -> tuple[dict[str, np.ndarray], list[str], float]:
    indices, signs, genes, coverage = measured_module(data, rows)
    if len(indices) == 0:
        return {}, genes, coverage
    mean = data.log_cpm[:, reference_columns].mean(axis=1)
    sd = data.log_cpm[:, reference_columns].std(axis=1, ddof=1)
    safe_sd = sd.copy()
    invariant = ~np.isfinite(safe_sd) | (safe_sd == 0)
    safe_sd[invariant] = 1
    standardized = (data.log_cpm[:, selected_columns] - mean[:, None]) / safe_sd[:, None]
    standardized[invariant, :] = 0
    zmean = (standardized[indices, :] * signs[:, None]).mean(axis=0)
    sing = signed_singscore(data.ranks[:, selected_columns], indices, signs)
    return {"singscore": sing, "standardized_mean": zmean}, genes, coverage


def hedges_g_only(control: np.ndarray, case: np.ndarray) -> float:
    control = np.asarray(control, dtype=float)
    case = np.asarray(case, dtype=float)
    df = len(control) + len(case) - 2
    if len(control) < 2 or len(case) < 2 or df <= 0:
        return float("nan")
    pooled_var = ((len(control) - 1) * control.var(ddof=1) + (len(case) - 1) * case.var(ddof=1)) / df
    if not np.isfinite(pooled_var) or pooled_var <= 0:
        return float("nan")
    correction = 1 - 3 / (4 * df - 1)
    return float(correction * (case.mean() - control.mean()) / math.sqrt(pooled_var))


def binary_auc(labels: np.ndarray, scores: np.ndarray) -> float:
    labels = np.asarray(labels, dtype=bool)
    scores = np.asarray(scores, dtype=float)
    n_pos = int(labels.sum())
    n_neg = len(labels) - n_pos
    if n_pos == 0 or n_neg == 0:
        return float("nan")
    ranks = rankdata(scores, method="average")
    return float((ranks[labels].sum() - n_pos * (n_pos + 1) / 2) / (n_pos * n_neg))


def lineage_permutation_p(
    manifest: pd.DataFrame,
    scores: np.ndarray,
    target_lineage: str,
    seed: int,
) -> float:
    labels = manifest["harmonized_lineage"].astype(str).to_numpy()
    observed = scores[labels == target_lineage].mean() - scores[labels != target_lineage].mean()
    donor_groups = [np.asarray(indices, dtype=int) for _, indices in manifest.groupby("donor_id", sort=True).indices.items()]
    rng = np.random.default_rng(seed)
    extreme = 0
    for _ in range(N_RESAMPLES):
        permuted = labels.copy()
        for indices in donor_groups:
            permuted[indices] = rng.permutation(permuted[indices])
        difference = scores[permuted == target_lineage].mean() - scores[permuted != target_lineage].mean()
        extreme += difference >= observed - 1e-15
    return (extreme + 1) / (N_RESAMPLES + 1)


def identity_controls(repo: Path, controls: pd.DataFrame, output: Path) -> tuple[pd.DataFrame, pd.DataFrame, pd.DataFrame]:
    modules = controls[controls["control_class"].eq("lineage_identity")]
    score_rows: list[dict[str, object]] = []
    effect_rows: list[dict[str, object]] = []
    performance_rows: list[dict[str, object]] = []
    module_lineage = modules.drop_duplicates("control_id").set_index("control_id")["lineage"].to_dict()
    for key in COHORTS:
        data = load_cohort(repo, key)
        selected = data.manifest.index[data.manifest["eligible_20"]].to_numpy()
        manifest = data.manifest.loc[selected].copy().reset_index().rename(columns={"index": "manifest_index"})
        reference = selected
        cohort_scores: dict[tuple[str, str], np.ndarray] = {}
        for control_id, rows in modules.groupby("control_id", sort=True):
            scores, measured, coverage = module_scores(data, rows, selected, reference)
            evaluable = len(measured) >= 5 and coverage >= 0.60
            if not evaluable:
                continue
            target = str(rows["lineage"].iloc[0])
            for method, values in scores.items():
                cohort_scores[(control_id, method)] = values
                labels = manifest["harmonized_lineage"].eq(target).to_numpy()
                effect_rows.append(
                    {
                        "dataset_id": data.dataset_id,
                        "assay": data.assay,
                        "annotation": data.annotation,
                        "identity_evidence": data.identity_evidence,
                        "control_id": control_id,
                        "target_lineage": target,
                        "score_method": method,
                        "measured_genes": len(measured),
                        "module_coverage": coverage,
                        "n_target": int(labels.sum()),
                        "n_off_target": int((~labels).sum()),
                        "hedges_g_target_vs_off": hedges_g_only(values[~labels], values[labels]),
                        "auc_one_vs_rest": binary_auc(labels, values),
                        "donor_stratified_permutation_p_one_sided": lineage_permutation_p(
                            manifest,
                            values,
                            target,
                            stable_seed(data.dataset_id, control_id, method, "identity"),
                        ),
                    }
                )
                for local, row in manifest.iterrows():
                    score_rows.append(
                        {
                            "dataset_id": data.dataset_id,
                            "manifest_index": int(row["manifest_index"]),
                            "donor_id": row["donor_id"],
                            "observed_lineage": row["harmonized_lineage"],
                            "control_id": control_id,
                            "target_lineage": target,
                            "score_method": method,
                            "score": float(values[local]),
                        }
                    )
        for method in ("singscore", "standardized_mean"):
            method_modules = sorted(control_id for control_id, m in cohort_scores if m == method)
            if len(method_modules) != 3:
                continue
            matrix = np.column_stack([cohort_scores[(control_id, method)] for control_id in method_modules])
            predicted_modules = np.asarray(method_modules, dtype=object)[np.argmax(matrix, axis=1)]
            predicted_lineage = np.asarray([module_lineage[module] for module in predicted_modules])
            macro_auc = float(
                np.nanmean(
                    [
                        binary_auc(
                            manifest["harmonized_lineage"].eq(module_lineage[module]).to_numpy(),
                            cohort_scores[(module, method)],
                        )
                        for module in method_modules
                    ]
                )
            )
            accuracy = float((predicted_lineage == manifest["harmonized_lineage"].to_numpy()).mean())
            performance_rows.append(
                {
                    "dataset_id": data.dataset_id,
                    "assay": data.assay,
                    "annotation": data.annotation,
                    "identity_evidence": data.identity_evidence,
                    "score_method": method,
                    "n_pseudobulks": len(manifest),
                    "top_score_lineage_accuracy": accuracy,
                    "macro_one_vs_rest_auc": macro_auc,
                    "passes_frozen_positive_control_threshold": accuracy >= 0.80 and macro_auc >= 0.90,
                }
            )
        del data
    scores_frame = pd.DataFrame(score_rows)
    effects = pd.DataFrame(effect_rows)
    effects["fdr_bh"] = bh_adjust(effects["donor_stratified_permutation_p_one_sided"])
    performance = pd.DataFrame(performance_rows)
    scores_frame.to_csv(
        output / "identity_control_scores.csv.gz",
        index=False,
        compression={"method": "gzip", "mtime": 0},
    )
    effects.to_csv(output / "identity_control_effects.csv", index=False)
    performance.to_csv(output / "identity_control_performance.csv", index=False)
    return scores_frame, effects, performance


def exact_label_test(control: np.ndarray, case: np.ndarray) -> float:
    control = np.asarray(control, dtype=float)
    case = np.asarray(case, dtype=float)
    values = np.concatenate([control, case])
    n_case = len(case)
    observed = case.mean() - control.mean()
    extreme = 0
    total = 0
    for case_indices in combinations(range(len(values)), n_case):
        mask = np.zeros(len(values), dtype=bool)
        mask[list(case_indices)] = True
        difference = values[mask].mean() - values[~mask].mean()
        extreme += abs(difference) >= abs(observed) - 1e-15
        total += 1
    return extreme / total


def watson_identity_diagnostic(
    repo: Path, identity_scores: pd.DataFrame, output: Path
) -> tuple[pd.DataFrame, pd.DataFrame]:
    scores = identity_scores[identity_scores["dataset_id"].eq("GSE210077_Watson6")].copy()
    module_to_lineage = scores.drop_duplicates("control_id").set_index("control_id")["target_lineage"].to_dict()
    data = load_cohort(repo, "Watson6")
    library = np.asarray(data.counts.sum(axis=0)).ravel()
    detected = np.asarray((data.counts > 0).sum(axis=0)).ravel()
    qc = data.manifest[["donor_id", "harmonized_lineage", "n_cells", "disease_group"]].copy()
    qc["library_size"] = library
    qc["detected_genes"] = detected
    rows: list[dict[str, object]] = []
    for (donor_id, method), values in scores.groupby(["donor_id", "score_method"], sort=True):
        pivot = values.pivot_table(index="observed_lineage", columns="control_id", values="score", aggfunc="first")
        if len(pivot) != 3 or pivot.shape[1] != 3:
            continue
        predictions = pivot.idxmax(axis=1).map(module_to_lineage)
        matched_margins: list[float] = []
        for observed_lineage, score_row in pivot.iterrows():
            matched_module = next(module for module, lineage in module_to_lineage.items() if lineage == observed_lineage)
            off = score_row.drop(labels=[matched_module])
            matched_margins.append(float(score_row[matched_module] - off.max()))
        donor_qc = qc[qc["donor_id"].eq(donor_id)]
        rows.append(
            {
                "donor_id": donor_id,
                "disease_group": donor_qc["disease_group"].iloc[0],
                "score_method": method,
                "lineage_top_score_accuracy": float((predictions.to_numpy() == pivot.index.to_numpy()).mean()),
                "mean_matched_minus_best_off_margin": float(np.mean(matched_margins)),
                "total_target_cells_or_nuclei": int(donor_qc["n_cells"].sum()),
                "median_lineage_library_size": float(donor_qc["library_size"].median()),
                "median_lineage_detected_genes": float(donor_qc["detected_genes"].median()),
            }
        )
    donors = pd.DataFrame(rows)
    donors.to_csv(output / "watson_identity_control_donor_audit.csv", index=False)
    contrast_rows: list[dict[str, object]] = []
    outcomes = [
        "lineage_top_score_accuracy",
        "mean_matched_minus_best_off_margin",
        "total_target_cells_or_nuclei",
        "median_lineage_library_size",
        "median_lineage_detected_genes",
    ]
    for method, values in donors.groupby("score_method", sort=True):
        for outcome in outcomes:
            healthy = values.loc[values["disease_group"].eq("healthy"), outcome].to_numpy(dtype=float)
            fibrosis = values.loc[values["disease_group"].eq("fibrosis"), outcome].to_numpy(dtype=float)
            contrast_rows.append(
                {
                    "score_method": method,
                    "outcome": outcome,
                    "n_healthy": len(healthy),
                    "n_fibrosis": len(fibrosis),
                    "healthy_median": float(np.median(healthy)),
                    "fibrosis_median": float(np.median(fibrosis)),
                    "fibrosis_minus_healthy": float(np.mean(fibrosis) - np.mean(healthy)),
                    "hedges_g_fibrosis_vs_healthy": hedges_g_only(healthy, fibrosis),
                    "exact_permutation_p_two_sided": exact_label_test(healthy, fibrosis),
                    "interpretation": "post-result cohort diagnostic; disease and technical context are confounded",
                }
            )
    contrasts = pd.DataFrame(contrast_rows)
    contrasts.to_csv(output / "watson_identity_control_contrasts.csv", index=False)
    del data
    return donors, contrasts


def representative_gate_lookup(repo: Path) -> dict[tuple[str, str], int]:
    effects = pd.read_csv(repo / "results" / "meta" / "cross_cohort_effect_matrix.csv")
    lookup: dict[tuple[str, str], int] = {}
    for cfg in COHORTS.values():
        subset = effects[
            effects["dataset_id"].eq(cfg["dataset_id"])
            & effects["contrast"].eq(cfg["contrast"])
        ]
        for lineage, rows in subset.groupby("lineage"):
            gates = rows["cell_gate"].dropna().astype(int).unique()
            if len(gates):
                lookup[(str(cfg["dataset_id"]), str(lineage))] = int(gates.max())
    return lookup


def disease_response_controls(
    repo: Path, controls: pd.DataFrame, output: Path
) -> tuple[pd.DataFrame, pd.DataFrame]:
    modules = controls[controls["control_class"].eq("disease_response")]
    gate_lookup = representative_gate_lookup(repo)
    score_rows: list[dict[str, object]] = []
    effect_rows: list[dict[str, object]] = []
    for key in COHORTS:
        data = load_cohort(repo, key)
        group = representative_groups(data)
        for control_id, rows in modules.groupby("control_id", sort=True):
            lineage = str(rows["lineage"].iloc[0])
            gate = gate_lookup.get((data.dataset_id, lineage))
            if gate is None:
                continue
            lineage_columns = data.manifest.index[data.manifest["harmonized_lineage"].eq(lineage)].to_numpy()
            reference = data.manifest.index[
                data.manifest["harmonized_lineage"].eq(lineage) & data.manifest["eligible_20"]
            ].to_numpy()
            scores, measured, coverage = module_scores(data, rows, lineage_columns, reference)
            if len(measured) < 5 or coverage < 0.60:
                continue
            eligibility = data.manifest[f"eligible_{gate}"]
            selected = data.manifest["harmonized_lineage"].eq(lineage) & eligibility & group.ne("excluded")
            selected_global = data.manifest.index[selected].to_numpy()
            global_to_local = {int(index): local for local, index in enumerate(lineage_columns)}
            selected_local = np.asarray([global_to_local[int(index)] for index in selected_global], dtype=int)
            selected_group = group.loc[selected_global]
            for method, values in scores.items():
                chosen = values[selected_local]
                control_values = chosen[selected_group.eq("control").to_numpy()]
                case_values = chosen[selected_group.eq("case").to_numpy()]
                if len(control_values) < 2 or len(case_values) < 2:
                    continue
                try:
                    stats = effect_statistics(
                        control_values,
                        case_values,
                        stable_seed(data.dataset_id, control_id, method, "disease_control"),
                    )
                    status = "evaluated"
                except ValueError as exc:
                    stats = {
                        "n_control": len(control_values),
                        "n_case": len(case_values),
                        "hedges_g": np.nan,
                        "robust_ci95_low": np.nan,
                        "robust_ci95_high": np.nan,
                        "permutation_p_two_sided": np.nan,
                        "permutation_mode": "not_evaluable",
                    }
                    status = f"not evaluable: {exc}"
                effect_rows.append(
                    {
                        "dataset_id": data.dataset_id,
                        "contrast": data.contrast,
                        "assay": data.assay,
                        "annotation": data.annotation,
                        "control_id": control_id,
                        "lineage": lineage,
                        "score_method": method,
                        "cell_gate": gate,
                        "measured_genes": len(measured),
                        "module_coverage": coverage,
                        "status": status,
                        **stats,
                    }
                )
                for local, global_index in zip(selected_local, selected_global):
                    score_rows.append(
                        {
                            "dataset_id": data.dataset_id,
                            "contrast": data.contrast,
                            "control_id": control_id,
                            "lineage": lineage,
                            "score_method": method,
                            "donor_id": data.manifest.loc[global_index, "donor_id"],
                            "group": group.loc[global_index],
                            "score": float(values[local]),
                        }
                    )
        del data
    scores_frame = pd.DataFrame(score_rows)
    effects = pd.DataFrame(effect_rows)
    effects["fdr_bh"] = bh_adjust(effects["permutation_p_two_sided"])
    effects["positive_fdr_005"] = effects["hedges_g"].gt(0) & effects["fdr_bh"].lt(0.05)
    scores_frame.to_csv(
        output / "disease_control_scores.csv.gz",
        index=False,
        compression={"method": "gzip", "mtime": 0},
    )
    effects.to_csv(output / "disease_control_effects.csv", index=False)
    return scores_frame, effects


def representative_effects(repo: Path) -> pd.DataFrame:
    effects = pd.read_csv(repo / "results" / "meta" / "cross_cohort_effect_matrix.csv")
    keep = pd.Series(False, index=effects.index)
    for cfg in COHORTS.values():
        keep |= effects["dataset_id"].eq(cfg["dataset_id"]) & effects["contrast"].eq(cfg["contrast"])
    frame = effects[keep].copy()
    frame["cohort_key"] = frame["dataset_id"].map(
        {str(cfg["dataset_id"]): key for key, cfg in COHORTS.items()}
    )
    return frame


def pairwise_decomposition(repo: Path, output: Path) -> tuple[pd.DataFrame, pd.DataFrame, pd.DataFrame]:
    effects = representative_effects(repo)
    effects["effect_z"] = effects.groupby(["dataset_id", "lineage", "score_method"])["hedges_g"].transform(
        lambda values: (values - values.mean()) / values.std(ddof=1) if values.std(ddof=1) > 0 else np.nan
    )
    rows: list[dict[str, object]] = []
    for (lineage, method, program_id), values in effects.groupby(["lineage", "score_method", "program_id"], sort=True):
        by_dataset = {dataset: row.iloc[0] for dataset, row in values.groupby("dataset_id")}
        for left, right in combinations(sorted(by_dataset), 2):
            a = by_dataset[left]
            b = by_dataset[right]
            endpoint_comparable = {
                str(a["endpoint_family"]),
                str(b["endpoint_family"]),
            }.issubset({"advanced_f3f4_vs_f0", "advanced_f3f4_vs_f0_reconstructed_sensitivity"})
            group_ns = np.asarray([a["n_control"], a["n_case"], b["n_control"], b["n_case"]], dtype=float)
            harmonic = float(len(group_ns) / np.sum(1.0 / group_ns))
            rows.append(
                {
                    "lineage": lineage,
                    "score_method": method,
                    "program_id": program_id,
                    "dataset_left": left,
                    "dataset_right": right,
                    "same_assay": str(a["assay"]) == str(b["assay"]),
                    "both_author_labelled": left in {"GSE202379", "GSE210077_Watson6"} and right in {"GSE202379", "GSE210077_Watson6"},
                    "comparable_advanced_endpoint": endpoint_comparable,
                    "includes_watson": "GSE210077_Watson6" in {left, right},
                    "coverage_abs_difference": abs(float(a["program_coverage"]) - float(b["program_coverage"])),
                    "minimum_program_coverage": min(float(a["program_coverage"]), float(b["program_coverage"])),
                    "harmonic_mean_group_donors": harmonic,
                    "minimum_group_donors": float(group_ns.min()),
                    "left_hedges_g": float(a["hedges_g"]),
                    "right_hedges_g": float(b["hedges_g"]),
                    "sign_discordant": np.sign(float(a["hedges_g"])) != np.sign(float(b["hedges_g"])),
                    "absolute_effect_difference": abs(float(a["hedges_g"]) - float(b["hedges_g"])),
                    "absolute_rank_z_difference": abs(float(a["effect_z"]) - float(b["effect_z"])),
                }
            )
    pairs = pd.DataFrame(rows)
    pairs.to_csv(output / "transfer_failure_pairwise_rows.csv", index=False)

    rng = np.random.default_rng(stable_seed("failure_decomposition_bootstrap"))
    program_ids = sorted(pairs["program_id"].unique())
    blocks = {program: pairs[pairs["program_id"].eq(program)] for program in program_ids}
    summary_rows: list[dict[str, object]] = []
    categorical = ["same_assay", "both_author_labelled", "comparable_advanced_endpoint", "includes_watson"]
    outcomes = ["sign_discordant", "absolute_effect_difference", "absolute_rank_z_difference"]
    for descriptor in categorical:
        for outcome in outcomes:
            positive = pairs.loc[pairs[descriptor].astype(bool), outcome].astype(float)
            negative = pairs.loc[~pairs[descriptor].astype(bool), outcome].astype(float)
            if positive.empty or negative.empty:
                continue
            observed = float(positive.median() - negative.median())
            boot = np.empty(N_RESAMPLES, dtype=float)
            for iteration in range(N_RESAMPLES):
                sampled = rng.choice(program_ids, size=len(program_ids), replace=True)
                frame = pd.concat([blocks[program] for program in sampled], ignore_index=True)
                pos = frame.loc[frame[descriptor].astype(bool), outcome].astype(float)
                neg = frame.loc[~frame[descriptor].astype(bool), outcome].astype(float)
                boot[iteration] = pos.median() - neg.median() if len(pos) and len(neg) else np.nan
            valid = boot[np.isfinite(boot)]
            summary_rows.append(
                {
                    "analysis_type": "categorical_median_contrast",
                    "descriptor": descriptor,
                    "outcome": outcome,
                    "n_rows": len(pairs),
                    "estimate_true_minus_false": observed,
                    "bootstrap_ci95_low": float(np.quantile(valid, 0.025)),
                    "bootstrap_ci95_high": float(np.quantile(valid, 0.975)),
                    "spearman_rho": np.nan,
                    "bootstrap_rho_ci95_low": np.nan,
                    "bootstrap_rho_ci95_high": np.nan,
                }
            )
    continuous = [
        "coverage_abs_difference",
        "minimum_program_coverage",
        "harmonic_mean_group_donors",
        "minimum_group_donors",
    ]
    for descriptor in continuous:
        for outcome in ["absolute_effect_difference", "absolute_rank_z_difference"]:
            observed = float(spearmanr(pairs[descriptor], pairs[outcome]).statistic)
            boot = np.empty(N_RESAMPLES, dtype=float)
            for iteration in range(N_RESAMPLES):
                sampled = rng.choice(program_ids, size=len(program_ids), replace=True)
                frame = pd.concat([blocks[program] for program in sampled], ignore_index=True)
                boot[iteration] = float(spearmanr(frame[descriptor], frame[outcome]).statistic)
            valid = boot[np.isfinite(boot)]
            summary_rows.append(
                {
                    "analysis_type": "continuous_spearman",
                    "descriptor": descriptor,
                    "outcome": outcome,
                    "n_rows": len(pairs),
                    "estimate_true_minus_false": np.nan,
                    "bootstrap_ci95_low": np.nan,
                    "bootstrap_ci95_high": np.nan,
                    "spearman_rho": observed,
                    "bootstrap_rho_ci95_low": float(np.quantile(valid, 0.025)),
                    "bootstrap_rho_ci95_high": float(np.quantile(valid, 0.975)),
                }
            )
    summaries = pd.DataFrame(summary_rows)
    summaries.to_csv(output / "transfer_failure_decomposition_summary.csv", index=False)

    prediction_rows: list[dict[str, object]] = []
    for row in effects.itertuples(index=False):
        train = effects[
            effects["program_id"].eq(row.program_id)
            & effects["lineage"].eq(row.lineage)
            & effects["score_method"].eq(row.score_method)
            & ~effects["dataset_id"].eq(row.dataset_id)
        ]
        if len(train) < 2 or not np.isfinite(row.effect_z):
            continue
        prediction_rows.append(
            {
                "held_out_dataset": row.dataset_id,
                "lineage": row.lineage,
                "score_method": row.score_method,
                "program_id": row.program_id,
                "observed_effect_z": float(row.effect_z),
                "predicted_effect_z_from_other_cohorts": float(train["effect_z"].mean()),
                "n_training_cohorts": int(train["dataset_id"].nunique()),
            }
        )
    predictions = pd.DataFrame(prediction_rows)
    predictions.to_csv(output / "leave_one_cohort_out_predictions.csv", index=False)
    pred_summary_rows: list[dict[str, object]] = []
    for keys, values in [(('ALL', 'ALL', 'ALL'), predictions)] + list(
        predictions.groupby(["held_out_dataset", "lineage", "score_method"], sort=True)
    ):
        if len(values) < 3:
            continue
        observed = values["observed_effect_z"].to_numpy(dtype=float)
        predicted = values["predicted_effect_z_from_other_cohorts"].to_numpy(dtype=float)
        r2 = 1 - np.sum((observed - predicted) ** 2) / np.sum((observed - observed.mean()) ** 2)
        pred_summary_rows.append(
            {
                "held_out_dataset": keys[0],
                "lineage": keys[1],
                "score_method": keys[2],
                "n_program_predictions": len(values),
                "spearman_rho": float(spearmanr(observed, predicted).statistic),
                "predictive_r_squared": float(r2),
            }
        )
    prediction_summary = pd.DataFrame(pred_summary_rows)
    prediction_summary.to_csv(output / "leave_one_cohort_out_prediction_summary.csv", index=False)
    return pairs, summaries, prediction_summary


def transfer_sensitivity_diagnostics(
    repo: Path, pairs: pd.DataFrame, output: Path
) -> tuple[pd.DataFrame, pd.DataFrame]:
    strata = {
        "all_pairs": pd.Series(True, index=pairs.index),
        "exclude_watson": ~pairs["includes_watson"].astype(bool),
        "minimum_coverage_ge_0_80": pairs["minimum_program_coverage"].ge(0.80),
        "minimum_coverage_ge_0_80_exclude_watson": pairs["minimum_program_coverage"].ge(0.80)
        & ~pairs["includes_watson"].astype(bool),
        "comparable_advanced_endpoint": pairs["comparable_advanced_endpoint"].astype(bool),
    }
    rows: list[dict[str, object]] = []
    for name, keep in strata.items():
        frame = pairs[keep].copy()
        pair_rhos: list[float] = []
        pair_sign: list[float] = []
        for _, values in frame.groupby(
            ["dataset_left", "dataset_right", "lineage", "score_method"], sort=True
        ):
            if len(values) < 3:
                continue
            pair_rhos.append(float(spearmanr(values["left_hedges_g"], values["right_hedges_g"]).statistic))
            pair_sign.append(float((~values["sign_discordant"].astype(bool)).mean()))
        rows.append(
            {
                "stratum": name,
                "n_program_pair_rows": len(frame),
                "n_unique_cohort_lineage_score_pairs": len(pair_rhos),
                "program_pair_sign_agreement": float((~frame["sign_discordant"].astype(bool)).mean()) if len(frame) else np.nan,
                "median_pairwise_program_spearman": float(np.median(pair_rhos)) if pair_rhos else np.nan,
                "median_pairwise_sign_agreement": float(np.median(pair_sign)) if pair_sign else np.nan,
                "median_absolute_effect_difference": float(frame["absolute_effect_difference"].median()) if len(frame) else np.nan,
                "median_absolute_rank_z_difference": float(frame["absolute_rank_z_difference"].median()) if len(frame) else np.nan,
            }
        )
    summary = pd.DataFrame(rows)
    summary.to_csv(output / "transfer_failure_stratified_summary.csv", index=False)

    effects = representative_effects(repo)
    effects = effects[~effects["dataset_id"].eq("GSE210077_Watson6")].copy()
    effects["effect_z"] = effects.groupby(["dataset_id", "lineage", "score_method"])["hedges_g"].transform(
        lambda values: (values - values.mean()) / values.std(ddof=1) if values.std(ddof=1) > 0 else np.nan
    )
    prediction_rows: list[dict[str, object]] = []
    for row in effects.itertuples(index=False):
        train = effects[
            effects["program_id"].eq(row.program_id)
            & effects["lineage"].eq(row.lineage)
            & effects["score_method"].eq(row.score_method)
            & ~effects["dataset_id"].eq(row.dataset_id)
        ]
        if len(train) < 2 or not np.isfinite(row.effect_z):
            continue
        prediction_rows.append(
            {
                "held_out_dataset": row.dataset_id,
                "lineage": row.lineage,
                "score_method": row.score_method,
                "program_id": row.program_id,
                "observed_effect_z": float(row.effect_z),
                "predicted_effect_z_from_other_non_watson_cohorts": float(train["effect_z"].mean()),
                "n_training_cohorts": int(train["dataset_id"].nunique()),
            }
        )
    predictions = pd.DataFrame(prediction_rows)
    predictions.to_csv(output / "leave_one_cohort_out_predictions_excluding_watson.csv", index=False)
    prediction_summaries: list[dict[str, object]] = []
    grouped = [(('ALL_NON_WATSON', 'ALL', 'ALL'), predictions)] + list(
        predictions.groupby(["held_out_dataset", "lineage", "score_method"], sort=True)
    )
    for keys, values in grouped:
        if len(values) < 3:
            continue
        observed = values["observed_effect_z"].to_numpy(dtype=float)
        predicted = values["predicted_effect_z_from_other_non_watson_cohorts"].to_numpy(dtype=float)
        denominator = np.sum((observed - observed.mean()) ** 2)
        r2 = 1 - np.sum((observed - predicted) ** 2) / denominator if denominator > 0 else np.nan
        prediction_summaries.append(
            {
                "held_out_dataset": keys[0],
                "lineage": keys[1],
                "score_method": keys[2],
                "n_program_predictions": len(values),
                "spearman_rho": float(spearmanr(observed, predicted).statistic),
                "predictive_r_squared": float(r2),
            }
        )
    prediction_summary = pd.DataFrame(prediction_summaries)
    prediction_summary.to_csv(output / "leave_one_cohort_out_prediction_summary_excluding_watson.csv", index=False)
    return summary, prediction_summary


def gene_effects_discovery(repo: Path, output: Path) -> tuple[pd.DataFrame, pd.DataFrame]:
    programs = pd.read_csv(repo / "literature" / "program_inventory.csv")
    recurrence = (
        programs.groupby(["cell_lineage", "gene_symbol"])["program_id"]
        .nunique()
        .rename("program_recurrence")
        .reset_index()
    )
    recurrence["gene_symbol"] = recurrence["gene_symbol"].astype(str).str.upper()
    recurrence = recurrence[recurrence["program_recurrence"].ge(2)].copy()

    gse202379 = load_cohort(repo, "GSE202379")
    gse244832 = load_cohort(repo, "GSE244832")
    rows: list[dict[str, object]] = []
    for lineage in LINEAGES:
        candidates = recurrence[recurrence["cell_lineage"].eq(lineage)]
        m1 = gse202379.manifest
        stage = pd.to_numeric(m1["Fibrosis.score..F0.4."], errors="coerce")
        selected1 = m1.index[
            m1["harmonized_lineage"].eq(lineage)
            & m1["eligible_30"]
            & ~m1["Disease.status"].eq("end stage")
            & stage.notna()
        ].to_numpy()
        m2 = gse244832.manifest
        group2 = representative_groups(gse244832)
        selected2 = m2.index[
            m2["harmonized_lineage"].eq(lineage) & m2["eligible_30"] & group2.ne("excluded")
        ].to_numpy()
        for candidate in candidates.itertuples(index=False):
            gene = str(candidate.gene_symbol)
            measured1 = gene in gse202379.gene_to_index
            measured2 = gene in gse244832.gene_to_index
            row: dict[str, object] = {
                "lineage": lineage,
                "gene_symbol": gene,
                "program_recurrence": int(candidate.program_recurrence),
                "measured_gse202379": measured1,
                "measured_gse244832": measured2,
            }
            if measured1 and measured2 and len(selected1) >= 10 and len(selected2) >= 4:
                idx1 = gse202379.gene_to_index[gene]
                idx2 = gse244832.gene_to_index[gene]
                values1 = gse202379.log_cpm[idx1, selected1]
                stages = stage.loc[selected1].to_numpy(dtype=float)
                values2 = gse244832.log_cpm[idx2, selected2]
                groups2 = group2.loc[selected2]
                control2 = values2[groups2.eq("control").to_numpy()]
                case2 = values2[groups2.eq("case").to_numpy()]
                row.update(
                    {
                        "gse202379_n_donors": len(selected1),
                        "gse202379_detection_rate": float(np.mean(gse202379.counts[idx1, selected1].toarray().ravel() > 0)),
                        "gse202379_stage_spearman_rho": float(spearmanr(stages, values1).statistic),
                        "gse244832_n_control": len(control2),
                        "gse244832_n_case": len(case2),
                        "gse244832_detection_rate": float(np.mean(gse244832.counts[idx2, selected2].toarray().ravel() > 0)),
                        "gse244832_mash_vs_normal_hedges_g": hedges_g_only(control2, case2),
                    }
                )
            rows.append(row)
    effects = pd.DataFrame(rows)
    effects["passes_measurement"] = (
        effects["measured_gse202379"].astype(bool)
        & effects["measured_gse244832"].astype(bool)
        & effects["gse202379_detection_rate"].ge(0.20)
        & effects["gse244832_detection_rate"].ge(0.20)
    )
    effects["passes_discovery_effects"] = (
        effects["passes_measurement"]
        & effects["gse202379_stage_spearman_rho"].ge(0.10)
        & effects["gse244832_mash_vs_normal_hedges_g"].ge(0.20)
    )
    effects["stage_effect_percentile"] = effects.groupby("lineage")["gse202379_stage_spearman_rho"].rank(pct=True)
    effects["mash_effect_percentile"] = effects.groupby("lineage")["gse244832_mash_vs_normal_hedges_g"].rank(pct=True)
    effects["mean_discovery_effect_percentile"] = effects[["stage_effect_percentile", "mash_effect_percentile"]].mean(axis=1)
    effects.to_csv(output / "core_discovery_gene_effects.csv", index=False)

    core_rows: list[dict[str, object]] = []
    for lineage, values in effects[effects["passes_discovery_effects"]].groupby("lineage", sort=True):
        ranked = values.sort_values(
            ["mean_discovery_effect_percentile", "gene_symbol"], ascending=[False, True]
        )
        if len(ranked) < 5:
            continue
        for rank, row in enumerate(ranked.head(5).itertuples(index=False), start=1):
            core_rows.append(
                {
                    "core_id": f"EXPLORATORY_{lineage.upper()}_CORE5",
                    "lineage": lineage,
                    "gene_symbol": row.gene_symbol,
                    "direction": "UP",
                    "core_rank": rank,
                    "program_recurrence": row.program_recurrence,
                    "gse202379_stage_spearman_rho": row.gse202379_stage_spearman_rho,
                    "gse244832_mash_vs_normal_hedges_g": row.gse244832_mash_vs_normal_hedges_g,
                    "mean_discovery_effect_percentile": row.mean_discovery_effect_percentile,
                }
            )
    cores = pd.DataFrame(core_rows)
    cores.to_csv(output / "minimal_core_membership.csv", index=False)
    del gse202379, gse244832
    return effects, cores


def validate_cores(
    repo: Path, cores: pd.DataFrame, output: Path
) -> tuple[pd.DataFrame, pd.DataFrame, pd.DataFrame]:
    if cores.empty:
        empty = pd.DataFrame()
        empty.to_csv(output / "minimal_core_validation_effects.csv", index=False)
        empty.to_csv(output / "minimal_core_random_benchmark.csv", index=False)
        empty.to_csv(output / "minimal_core_validation_summary.csv", index=False)
        return empty, empty, empty
    gate_lookup = representative_gate_lookup(repo)
    effect_rows: list[dict[str, object]] = []
    benchmark_rows: list[dict[str, object]] = []
    random_rows: list[dict[str, object]] = []
    membership_rows: list[dict[str, object]] = []
    held_out = ["GSE290642", "Watson6", "GSE181483"]
    for key in held_out:
        data = load_cohort(repo, key)
        group = representative_groups(data)
        for core_id, rows in cores.groupby("core_id", sort=True):
            lineage = str(rows["lineage"].iloc[0])
            gate = gate_lookup.get((data.dataset_id, lineage))
            if gate is None:
                continue
            module_rows = rows[["gene_symbol", "direction"]].copy()
            lineage_columns = data.manifest.index[data.manifest["harmonized_lineage"].eq(lineage)].to_numpy()
            reference = data.manifest.index[
                data.manifest["harmonized_lineage"].eq(lineage) & data.manifest["eligible_20"]
            ].to_numpy()
            scores, measured, coverage = module_scores(data, module_rows, lineage_columns, reference)
            if len(measured) < 4 or coverage < 0.80:
                continue
            selected = (
                data.manifest["harmonized_lineage"].eq(lineage)
                & data.manifest[f"eligible_{gate}"]
                & group.ne("excluded")
            )
            selected_global = data.manifest.index[selected].to_numpy()
            global_to_local = {int(index): local for local, index in enumerate(lineage_columns)}
            selected_local = np.asarray([global_to_local[int(index)] for index in selected_global], dtype=int)
            selected_group = group.loc[selected_global]
            control_local = selected_local[selected_group.eq("control").to_numpy()]
            case_local = selected_local[selected_group.eq("case").to_numpy()]
            if len(control_local) < 2 or len(case_local) < 2:
                continue
            real_effects: dict[str, float] = {}
            for method, values in scores.items():
                stats = effect_statistics(
                    values[control_local],
                    values[case_local],
                    stable_seed(data.dataset_id, core_id, method, "core_validation"),
                )
                real_effects[method] = float(stats["hedges_g"])
                effect_rows.append(
                    {
                        "dataset_id": data.dataset_id,
                        "contrast": data.contrast,
                        "assay": data.assay,
                        "annotation": data.annotation,
                        "core_id": core_id,
                        "lineage": lineage,
                        "score_method": method,
                        "cell_gate": gate,
                        "measured_genes": len(measured),
                        "core_coverage": coverage,
                        **stats,
                    }
                )

            indices, _, _, _ = measured_module(data, module_rows)
            mean_expression = data.log_cpm[:, reference].mean(axis=1)
            detection = np.asarray((data.counts[:, reference] > 0).mean(axis=1)).ravel()
            expr_decile = deciles(mean_expression)
            det_decile = deciles(detection)
            sd = data.log_cpm[:, reference].std(axis=1, ddof=1)
            invariant = ~np.isfinite(sd) | (sd == 0)
            safe_sd = sd.copy()
            safe_sd[invariant] = 1
            standardized = (data.log_cpm[:, lineage_columns] - mean_expression[:, None]) / safe_sd[:, None]
            standardized[invariant, :] = 0
            rng = np.random.default_rng(stable_seed(data.dataset_id, core_id, "matched_random"))
            modules, distances, target_positions = matched_modules(indices, expr_decile, det_decile, rng)
            random_sing, random_zmean = score_modules(modules, data.ranks[:, lineage_columns], standardized)
            for module_index in range(N_RANDOM_MODULES):
                for position in range(modules.shape[1]):
                    membership_rows.append(
                        {
                            "dataset_id": data.dataset_id,
                            "core_id": core_id,
                            "module_id": f"R{module_index + 1:04d}",
                            "target_gene": data.genes.iloc[indices[target_positions[module_index, position]]],
                            "random_gene": data.genes.iloc[modules[module_index, position]],
                            "bin_distance": int(distances[module_index, position]),
                        }
                    )
            for method, random_scores in (("singscore", random_sing), ("standardized_mean", random_zmean)):
                random_g = vectorized_hedges_g(random_scores, control_local, case_local)
                real_g = real_effects[method]
                q95 = float(np.quantile(random_g, 0.95))
                benchmark_rows.append(
                    {
                        "dataset_id": data.dataset_id,
                        "contrast": data.contrast,
                        "core_id": core_id,
                        "lineage": lineage,
                        "score_method": method,
                        "real_hedges_g": real_g,
                        "random_g_median": float(np.median(random_g)),
                        "random_g_95th_percentile": q95,
                        "real_effect_percentile": float(np.mean(random_g <= real_g)),
                        "empirical_p_one_sided": (1 + int(np.sum(random_g >= real_g))) / (N_RANDOM_MODULES + 1),
                        "above_random_95th_percentile": real_g > q95,
                        "random_modules": N_RANDOM_MODULES,
                    }
                )
                for module_index, random_g_value in enumerate(random_g, start=1):
                    random_rows.append(
                        {
                            "dataset_id": data.dataset_id,
                            "core_id": core_id,
                            "score_method": method,
                            "module_id": f"R{module_index:04d}",
                            "random_hedges_g": float(random_g_value),
                        }
                    )
        del data
    effects = pd.DataFrame(effect_rows)
    benchmarks = pd.DataFrame(benchmark_rows)
    effects.to_csv(output / "minimal_core_validation_effects.csv", index=False)
    benchmarks.to_csv(output / "minimal_core_random_benchmark.csv", index=False)
    pd.DataFrame(random_rows).to_csv(
        output / "minimal_core_random_effects.csv.gz",
        index=False,
        compression={"method": "gzip", "mtime": 0},
    )
    pd.DataFrame(membership_rows).to_csv(
        output / "minimal_core_random_membership.csv.gz",
        index=False,
        compression={"method": "gzip", "mtime": 0},
    )
    summary_rows: list[dict[str, object]] = []
    for core_id, values in effects.groupby("core_id", sort=True):
        positive = (
            values.assign(positive=values["hedges_g"].gt(0))
            .pivot_table(index="dataset_id", columns="score_method", values="positive", aggfunc="first", fill_value=False)
            .reindex(columns=["singscore", "standardized_mean"], fill_value=False)
            .all(axis=1)
        )
        random_values = benchmarks[benchmarks["core_id"].eq(core_id)]
        random_specific = (
            random_values.pivot_table(
                index="dataset_id",
                columns="score_method",
                values="above_random_95th_percentile",
                aggfunc="first",
                fill_value=False,
            )
            .reindex(columns=["singscore", "standardized_mean"], fill_value=False)
            .all(axis=1)
        )
        passes = int(positive.sum()) >= 2 and bool(random_specific.any())
        summary_rows.append(
            {
                "core_id": core_id,
                "lineage": values["lineage"].iloc[0],
                "held_out_datasets_evaluable": int(values["dataset_id"].nunique()),
                "held_out_datasets_positive_both_scores": int(positive.sum()),
                "positive_both_score_datasets": ";".join(positive.index[positive].astype(str)),
                "held_out_datasets_random_specific_both_scores": int(random_specific.sum()),
                "random_specific_both_score_datasets": ";".join(random_specific.index[random_specific].astype(str)),
                "held_out_directional_core": passes,
            }
        )
    summaries = pd.DataFrame(summary_rows)
    summaries.to_csv(output / "minimal_core_validation_summary.csv", index=False)
    return effects, benchmarks, summaries


def build_report_card(repo: Path, output: Path) -> pd.DataFrame:
    effects = representative_effects(repo)
    programs = pd.read_csv(repo / "literature" / "program_lineage_matrix.csv")[["program_id", "cell_lineage"]].drop_duplicates()
    programs = programs.rename(columns={"cell_lineage": "lineage"})
    random_frames = [pd.read_csv(path) for path in sorted((repo / "results" / "random_controls").glob("*_random_module_benchmark.csv"))]
    random = pd.concat(random_frames, ignore_index=True, sort=False)
    representative_pairs = {(str(cfg["dataset_id"]), str(cfg["contrast"])) for cfg in COHORTS.values()}
    random = random[
        [
            (str(dataset), str(contrast)) in representative_pairs
            for dataset, contrast in zip(random["dataset_id"], random["contrast"])
        ]
    ].copy()
    evidence = pd.read_csv(repo / "results" / "exploratory" / "program_evidence_attrition_matrix.csv")
    trends = pd.read_csv(repo / "results" / "exploratory" / "gse202379_stage_trends.csv")
    trends = trends[trends["analysis_set"].eq("non_end_stage_primary_exploratory")]
    classifications = pd.read_csv(repo / "results" / "meta" / "program_classification_table.csv")
    rows: list[dict[str, object]] = []
    for program in programs.itertuples(index=False):
        subset = effects[effects["program_id"].eq(program.program_id)].copy()
        coverage_by_cohort = subset.groupby("dataset_id")["program_coverage"].first()
        median_coverage = float(coverage_by_cohort.median()) if len(coverage_by_cohort) else 0.0
        minimum_coverage = float(coverage_by_cohort.min()) if len(coverage_by_cohort) else 0.0
        method_pivot = subset.pivot_table(
            index=["dataset_id", "lineage"], columns="score_method", values="hedges_g", aggfunc="first"
        ).dropna()
        method_reliability = float(
            (np.sign(method_pivot.get("singscore")) == np.sign(method_pivot.get("standardized_mean"))).mean()
        ) if len(method_pivot) and {"singscore", "standardized_mean"}.issubset(method_pivot.columns) else 0.0
        directional = float(subset["hedges_g"].gt(0).mean()) if len(subset) else 0.0
        random_subset = random[random["program_id"].eq(program.program_id)]
        random_pivot = random_subset.pivot_table(
            index=["dataset_id", "lineage"],
            columns="score_method",
            values="above_random_95th_percentile",
            aggfunc="first",
            fill_value=False,
        ).reindex(columns=["singscore", "standardized_mean"], fill_value=False)
        random_specificity = float(random_pivot.all(axis=1).mean()) if len(random_pivot) else 0.0
        evidence_row = evidence[evidence["program_id"].eq(program.program_id)].iloc[0]
        trend_row = trends[trends["program_id"].eq(program.program_id)]
        class_row = classifications[classifications["program_id"].eq(program.program_id)].iloc[0]
        stage_pass = bool(trend_row["dual_score_positive_fdr_005"].astype(bool).any()) if len(trend_row) else False
        meta_pass = bool(class_row["advanced_sensitivity_fixed_ci_excludes_zero_both_scores"])
        endpoint_components = {
            "positive_formal_primary_interval_both_scores": bool(evidence_row["positive_primary_interval_both_scores"]),
            "dual_score_fdr_positive_stage_trend": stage_pass,
            "advanced_meta_positive_interval_both_scores": meta_pass,
            "assay_robust": bool(evidence_row["assay_robust"]),
        }
        measurement_score = 10 * np.clip(median_coverage, 0, 1) + 10 * np.clip(minimum_coverage, 0, 1)
        method_score = 20 * method_reliability
        directional_score = 20 * directional
        random_score = 20 * random_specificity
        endpoint_score = 5 * sum(endpoint_components.values())
        rows.append(
            {
                "program_id": program.program_id,
                "lineage": program.lineage,
                "representative_cohorts_evaluable": int(subset["dataset_id"].nunique()),
                "median_program_coverage": median_coverage,
                "minimum_program_coverage": minimum_coverage,
                "score_method_sign_agreement_fraction": method_reliability,
                "expected_positive_effect_fraction": directional,
                "random_specific_both_scores_context_fraction": random_specificity,
                **endpoint_components,
                "measurement_domain_0_20": measurement_score,
                "score_method_domain_0_20": method_score,
                "directional_transfer_domain_0_20": directional_score,
                "matched_random_specificity_domain_0_20": random_score,
                "endpoint_evidence_domain_0_20": endpoint_score,
                "transportability_readiness_total_0_100": measurement_score
                + method_score
                + directional_score
                + random_score
                + endpoint_score,
                "interpretation": "exploratory prioritization index; not a clinical score or replication definition",
            }
        )
    report_card = pd.DataFrame(rows).sort_values(
        ["transportability_readiness_total_0_100", "program_id"], ascending=[False, True]
    )
    report_card.to_csv(output / "program_transportability_report_card.csv", index=False)
    return report_card


def write_report(
    repo: Path,
    identity_performance: pd.DataFrame,
    disease_effects: pd.DataFrame,
    decomposition: pd.DataFrame,
    prediction_summary: pd.DataFrame,
    transfer_strata: pd.DataFrame,
    prediction_without_watson: pd.DataFrame,
    watson_contrasts: pd.DataFrame,
    discovery: pd.DataFrame,
    cores: pd.DataFrame,
    core_summary: pd.DataFrame,
    report_card: pd.DataFrame,
) -> None:
    noncircular = identity_performance[identity_performance["annotation"].eq("author")]
    pass_by_cohort = noncircular.groupby("dataset_id")["passes_frozen_positive_control_threshold"].all()
    disease_dual = (
        disease_effects.pivot_table(
            index=["dataset_id", "control_id"],
            columns="score_method",
            values="positive_fdr_005",
            aggfunc="first",
            fill_value=False,
        )
        .reindex(columns=["singscore", "standardized_mean"], fill_value=False)
        .all(axis=1)
    )
    overall_prediction = prediction_summary[prediction_summary["held_out_dataset"].eq("ALL")]
    lines = [
        "# Deep transportability benchmark results",
        "",
        "Date: 2026-08-31",
        "",
        "All results below are exploratory and preserve the frozen primary classifications.",
        "",
        "## A. Positive-negative control ladder",
        "",
        f"- Non-circular author-label identity cohorts evaluated: {len(pass_by_cohort)}; cohorts passing the frozen dual-method accuracy/AUROC threshold: {int(pass_by_cohort.sum())}/{len(pass_by_cohort)}.",
        f"- Across the non-circular identity rows, median top-score accuracy was {noncircular['top_score_lineage_accuracy'].median():.3f} and median macro AUROC was {noncircular['macro_one_vs_rest_auc'].median():.3f}.",
        f"- Canonical disease-response cohort-module comparisons positive with both scores after FDR control: {int(disease_dual.sum())}/{len(disease_dual)}.",
        "- Matched random modules, direction randomization and donor-label permutations remain negative-control references; a positive identity result demonstrates technical sensitivity, not disease-program transportability.",
        "",
        "## B. Quantitative failure decomposition",
        "",
    ]
    if len(overall_prediction):
        row = overall_prediction.iloc[0]
        lines.extend(
            [
                f"- Leave-one-cohort-out prediction of held-out program ordering: Spearman rho {row['spearman_rho']:.3f}; predictive R-squared {row['predictive_r_squared']:.3f} across {int(row['n_program_predictions'])} predictions.",
                "- Descriptor contrasts are descriptive because assay, endpoint, annotation and etiology are partly confounded. Bootstrap intervals are provided in `results/deep_benchmark/transfer_failure_decomposition_summary.csv`.",
            ]
        )
    non_watson_prediction = prediction_without_watson[
        prediction_without_watson["held_out_dataset"].eq("ALL_NON_WATSON")
    ]
    if len(non_watson_prediction):
        row = non_watson_prediction.iloc[0]
        lines.append(
            f"- After fully excluding Watson from both training and testing, leave-one-cohort-out Spearman rho was {row['spearman_rho']:.3f} and predictive R-squared was {row['predictive_r_squared']:.3f} across {int(row['n_program_predictions'])} predictions."
        )
    high_coverage = transfer_strata[
        transfer_strata["stratum"].eq("minimum_coverage_ge_0_80_exclude_watson")
    ]
    if len(high_coverage):
        row = high_coverage.iloc[0]
        lines.append(
            f"- In program pairs with at least 80% coverage after excluding Watson, sign agreement was {row['program_pair_sign_agreement']:.3f} and median program-rank rho was {row['median_pairwise_program_spearman']:.3f}."
        )
    watson_margin = watson_contrasts[
        watson_contrasts["outcome"].eq("mean_matched_minus_best_off_margin")
        & watson_contrasts["score_method"].eq("standardized_mean")
    ]
    if len(watson_margin):
        row = watson_margin.iloc[0]
        lines.append(
            f"- Triggered Watson audit: the standardized-mean matched-lineage margin changed by {row['fibrosis_minus_healthy']:.3f} in fibrotic versus healthy donors (exact P={row['exact_permutation_p_two_sided']:.3f}); disease and technical context are inseparable in this 3+3 cohort."
        )
    strongest = decomposition.copy()
    strongest["magnitude"] = strongest["estimate_true_minus_false"].abs().fillna(strongest["spearman_rho"].abs())
    if len(strongest):
        row = strongest.sort_values("magnitude", ascending=False).iloc[0]
        value = row["estimate_true_minus_false"] if pd.notna(row["estimate_true_minus_false"]) else row["spearman_rho"]
        lines.append(f"- Largest prespecified descriptor association by absolute estimate: {row['descriptor']} with {row['outcome']} ({float(value):.3f}).")
    lines.extend(["", "## C. Frozen discovery and held-out minimal cores", ""])
    for lineage in LINEAGES:
        subset = discovery[discovery["lineage"].eq(lineage)]
        selected = subset[subset["passes_discovery_effects"]]
        core = cores[cores["lineage"].eq(lineage)] if not cores.empty else pd.DataFrame()
        if core.empty:
            lines.append(f"- {lineage}: {len(selected)} genes passed the frozen discovery gates; fewer than five were available, so no core was declared.")
        else:
            genes = ", ".join(core.sort_values("core_rank")["gene_symbol"])
            validation = core_summary[core_summary["lineage"].eq(lineage)] if not core_summary.empty else pd.DataFrame()
            if len(validation):
                status = "passed" if bool(validation.iloc[0]["held_out_directional_core"]) else "did not pass"
                lines.append(f"- {lineage}: frozen Core5 = {genes}; it {status} the held-out directional/random-specific rule.")
            else:
                lines.append(f"- {lineage}: frozen Core5 = {genes}; no held-out validation row was evaluable.")
    lines.extend(
        [
            "",
            "## D. Program transportability report card",
            "",
            "The report card is a transparent prioritization resource, not a clinical score. The five highest totals were:",
            "",
        ]
    )
    for row in report_card.head(5).itertuples(index=False):
        lines.append(
            f"- {row.program_id}: {row.transportability_readiness_total_0_100:.1f}/100 "
            f"(measurement {row.measurement_domain_0_20:.1f}, method {row.score_method_domain_0_20:.1f}, "
            f"direction {row.directional_transfer_domain_0_20:.1f}, random specificity {row.matched_random_specificity_domain_0_20:.1f}, "
            f"endpoint evidence {row.endpoint_evidence_domain_0_20:.1f})."
        )
    lines.extend(
        [
            "",
            "## Interpretation boundary",
            "",
            "These analyses can convert the project from a simple null report into a calibrated transportability benchmark only to the extent that technical positive controls pass and failure persists in held-out prediction. They do not establish a causal mechanism, universal scar state, diagnostic biomarker or treatment target.",
            "",
        ]
    )
    (repo / "reports" / "deep_transportability_benchmark_results.md").write_text("\n".join(lines), encoding="utf-8")


def main() -> None:
    repo = Path(__file__).resolve().parents[1]
    output = repo / "results" / "deep_benchmark"
    output.mkdir(parents=True, exist_ok=True)
    controls = pd.read_csv(repo / "config" / "deep_benchmark_control_programs.csv")

    identity_scores, identity_effects, identity_performance = identity_controls(repo, controls, output)
    _, watson_contrasts = watson_identity_diagnostic(repo, identity_scores, output)
    _, disease_effects = disease_response_controls(repo, controls, output)
    pairs, decomposition, prediction_summary = pairwise_decomposition(repo, output)
    transfer_strata, prediction_without_watson = transfer_sensitivity_diagnostics(repo, pairs, output)
    discovery, cores = gene_effects_discovery(repo, output)
    core_effects, core_benchmarks, core_summary = validate_cores(repo, cores, output)
    report_card = build_report_card(repo, output)
    write_report(
        repo,
        identity_performance,
        disease_effects,
        decomposition,
        prediction_summary,
        transfer_strata,
        prediction_without_watson,
        watson_contrasts,
        discovery,
        cores,
        core_summary,
        report_card,
    )

    files = sorted(path for path in output.iterdir() if path.is_file())
    summary = {
        "frozen_plan": "reports/deep_transportability_benchmark_plan_2026-08-31.md",
        "seed": SEED,
        "resamples": N_RESAMPLES,
        "identity_effect_rows": len(identity_effects),
        "identity_performance_rows": len(identity_performance),
        "noncircular_identity_cohorts_passing_both_methods": int(
            identity_performance[identity_performance["annotation"].eq("author")]
            .groupby("dataset_id")["passes_frozen_positive_control_threshold"]
            .all()
            .sum()
        ),
        "disease_control_effect_rows": len(disease_effects),
        "pairwise_failure_rows": len(pairs),
        "core_discovery_gene_rows": len(discovery),
        "lineage_cores_declared": int(cores["core_id"].nunique()) if not cores.empty else 0,
        "held_out_directional_cores": int(core_summary["held_out_directional_core"].sum()) if not core_summary.empty else 0,
        "report_card_programs": len(report_card),
        "core_validation_effect_rows": len(core_effects),
        "core_random_benchmark_rows": len(core_benchmarks),
        "output_sha256": {path.name: sha256(path) for path in files},
        "interpretation_boundary": "exploratory only; frozen primary program classifications unchanged",
    }
    log_path = repo / "results" / "logs" / "deep_transportability_benchmark_run.json"
    log_path.write_text(json.dumps(summary, indent=2) + "\n", encoding="utf-8")
    print(json.dumps(summary, indent=2))


if __name__ == "__main__":
    main()
