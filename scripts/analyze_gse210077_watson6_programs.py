from __future__ import annotations

import hashlib
import json
from pathlib import Path

import numpy as np
import pandas as pd
from scipy.io import mmread
from scipy.stats import rankdata

from analyze_gse202379_programs import effect_statistics, singscore_up, stable_seed
from audit_gse210077_watson6_gates import contrast_groups


def main() -> None:
    repo = Path(__file__).resolve().parents[1]
    interim = repo / "data" / "interim" / "GSE210077_Watson6"
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
    coverage = pd.read_csv(repo / "results" / "qc" / "gse210077_watson6_program_coverage.csv")
    coverage_lookup = coverage.set_index("program_id").to_dict("index")
    score_rows: list[dict[str, object]] = []
    score_arrays: dict[tuple[str, str], np.ndarray] = {}

    for lineage, lineage_manifest in manifest.groupby("harmonized_lineage", sort=False):
        columns = lineage_manifest.index.to_numpy()
        reference_columns = lineage_manifest.index[lineage_manifest["eligible_20"]].to_numpy()
        means = log_cpm[:, reference_columns].mean(axis=1)
        standard_deviations = log_cpm[:, reference_columns].std(axis=1, ddof=1)
        invariant = standard_deviations == 0
        safe_sd = standard_deviations.copy()
        safe_sd[invariant] = 1
        standardized = (log_cpm[:, columns] - means[:, None]) / safe_sd[:, None]
        standardized[invariant, :] = 0
        for program_id, rows in programs[programs["cell_lineage"].eq(lineage)].groupby("program_id"):
            coverage_info = coverage_lookup[program_id]
            if coverage_info["coverage_tier"] == "not_evaluated":
                continue
            measured_genes = sorted(set(rows["gene_symbol"].astype(str).str.upper()) & set(gene_to_index))
            gene_indices = np.array([gene_to_index[gene] for gene in measured_genes], dtype=int)
            sing = singscore_up(ranks[:, columns], gene_indices)
            zmean = standardized[gene_indices, :].mean(axis=0)
            for score_method, values in (("singscore", sing), ("standardized_mean", zmean)):
                full = np.full(len(manifest), np.nan)
                full[columns] = values
                score_arrays[(program_id, score_method)] = full
                for local_index, manifest_index in enumerate(columns):
                    score_rows.append(
                        {
                            "dataset_id": "GSE210077_Watson6",
                            "program_id": program_id,
                            "lineage": lineage,
                            "score_method": score_method,
                            "score": float(values[local_index]),
                            "measured_program_genes": len(measured_genes),
                            "program_coverage": coverage_info["coverage"],
                            "coverage_tier": coverage_info["coverage_tier"],
                            **manifest.loc[manifest_index].to_dict(),
                            "library_size": int(library_sizes[manifest_index]),
                        }
                    )
    scores = pd.DataFrame(score_rows)
    sensitivity_dir = repo / "results" / "sensitivity"
    score_path = sensitivity_dir / "gse210077_watson6_donor_program_scores.csv.gz"
    scores.to_csv(score_path, index=False, compression="gzip")

    gates = pd.read_csv(repo / "results" / "qc" / "gse210077_watson6_donor_gate_summary.csv")
    gate_status = gates.drop_duplicates(["contrast", "lineage"]).set_index(["contrast", "lineage"])
    effect_rows: list[dict[str, object]] = []
    for contrast, groups in contrast_groups(manifest).items():
        for program_id, program_rows in programs.groupby("program_id"):
            lineage = str(program_rows["cell_lineage"].iloc[0])
            coverage_info = coverage_lookup[program_id]
            if coverage_info["coverage_tier"] == "not_evaluated":
                continue
            status = gate_status.loc[(contrast, lineage)]
            if status["formal_30_cell_gate"] == "PASS":
                eligibility_column = "eligible_30"
                cell_gate = 30
            elif status["formal_20_cell_gate"] == "PASS":
                eligibility_column = "eligible_20"
                cell_gate = 20
            else:
                continue
            selected = manifest["harmonized_lineage"].eq(lineage) & manifest[eligibility_column] & groups.ne("excluded")
            selected_indices = manifest.index[selected]
            selected_groups = groups.loc[selected_indices]
            for score_method in ("singscore", "standardized_mean"):
                values = score_arrays[(program_id, score_method)][selected_indices]
                control = values[selected_groups.eq("control").to_numpy()]
                case = values[selected_groups.eq("case").to_numpy()]
                if len(control) != 3 or len(case) != 3:
                    raise RuntimeError("Watson effect must use exactly three donors per group")
                stats = effect_statistics(
                    control,
                    case,
                    stable_seed("GSE210077_Watson6", contrast, program_id, score_method),
                )
                effect_rows.append(
                    {
                        "dataset_id": "GSE210077_Watson6",
                        "contrast": contrast,
                        "program_id": program_id,
                        "lineage": lineage,
                        "score_method": score_method,
                        "analysis_tier": "sensitivity",
                        "cell_gate": cell_gate,
                        "coverage_tier": coverage_info["coverage_tier"],
                        "program_coverage": coverage_info["coverage"],
                        "expected_direction": "positive_case_higher",
                        "endpoint_limitation": "case donors are one each F2, F3, F4; not a cirrhosis cohort",
                        **stats,
                    }
                )
    effects = pd.DataFrame(effect_rows).sort_values(["lineage", "program_id", "score_method"])
    effect_path = sensitivity_dir / "gse210077_watson6_sensitivity_effects.csv"
    effects.to_csv(effect_path, index=False)
    run_summary = {
        "dataset_id": "GSE210077_Watson6",
        "normalization": "log2(CPM + 1) within donor-lineage pseudobulk",
        "score_methods": ["direction-aware singscore", "standardized signed mean"],
        "standardization_reference": "all six lineage donors passing 20-cell gate; outcome-blind",
        "feature_space": f"{len(genes)} genes measured in all six donors",
        "endpoint": "three healthy versus F2/F3/F4 fibrosis; sensitivity only",
        "effect_rows": len(effects),
        "effect_file_sha256": hashlib.sha256(effect_path.read_bytes()).hexdigest().upper(),
        "score_file_sha256": hashlib.sha256(score_path.read_bytes()).hexdigest().upper(),
    }
    (repo / "results" / "logs" / "gse210077_watson6_analysis_run.json").write_text(json.dumps(run_summary, indent=2) + "\n", encoding="utf-8")
    print(effects[["lineage", "program_id", "score_method", "n_control", "n_case", "hedges_g", "robust_ci95_low", "robust_ci95_high", "permutation_p_two_sided"]].to_string(index=False))
    print(json.dumps(run_summary, indent=2))


if __name__ == "__main__":
    main()
