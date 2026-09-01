from __future__ import annotations

import hashlib
import json
import math
from itertools import combinations
from pathlib import Path

import gseapy as gp
import numpy as np
import pandas as pd
from scipy.io import mmread
from scipy.stats import spearmanr

from analyze_gse202379_programs import stable_seed
from audit_gse202379_gates import contrast_groups as gse202379_groups


PERMUTATIONS = 1_000
MIN_PATHWAY_SIZE = 15
MAX_PATHWAY_SIZE = 500
MIN_GENE_DETECTION = 0.20


def sha256(path: Path) -> str:
    digest = hashlib.sha256()
    with path.open("rb") as handle:
        for block in iter(lambda: handle.read(1024 * 1024), b""):
            digest.update(block)
    return digest.hexdigest().upper()


def vectorized_hedges_g(control: np.ndarray, case: np.ndarray) -> np.ndarray:
    n_control = control.shape[1]
    n_case = case.shape[1]
    degrees_freedom = n_control + n_case - 2
    pooled = (
        (n_control - 1) * control.var(axis=1, ddof=1)
        + (n_case - 1) * case.var(axis=1, ddof=1)
    ) / degrees_freedom
    correction = 1 - 3 / (4 * degrees_freedom - 1)
    with np.errstate(divide="ignore", invalid="ignore"):
        effect = correction * (case.mean(axis=1) - control.mean(axis=1)) / np.sqrt(pooled)
    effect[~np.isfinite(effect)] = np.nan
    return effect


def contexts(repo: Path) -> list[dict[str, object]]:
    result: list[dict[str, object]] = []

    manifest = pd.read_csv(repo / "data" / "interim" / "GSE202379" / "donor_lineage_manifest.csv")
    groups = gse202379_groups(manifest)
    result.append(
        {
            "dataset_id": "GSE202379",
            "interim": "GSE202379",
            "contrast": "advanced_f3f4_vs_f0_non_end_stage",
            "manifest": manifest,
            "groups": groups["advanced_f3f4_vs_f0_non_end_stage"],
        }
    )

    manifest = pd.read_csv(repo / "data" / "interim" / "GSE244832" / "donor_lineage_manifest.csv")
    result.append(
        {
            "dataset_id": "GSE244832",
            "interim": "GSE244832",
            "contrast": "mash_f2f4_group_vs_normal_sensitivity",
            "manifest": manifest,
            "groups": manifest["disease_group"].map({"normal": "control", "MASH": "case"}).fillna("excluded"),
        }
    )

    manifest = pd.read_csv(repo / "data" / "interim" / "GSE256398" / "donor_lineage_manifest.csv")
    for contrast, controls, cases in (
        ("mash_cirrhosis_vs_healthy", {"healthy"}, {"mash_cirrhosis"}),
        ("alcohol_cirrhosis_vs_healthy", {"healthy"}, {"alcohol_cirrhosis"}),
        ("mash_fibrosis_vs_masld_f0", {"masld_f0"}, {"mash_fibrosis"}),
    ):
        result.append(
            {
                "dataset_id": "GSE256398_human",
                "interim": "GSE256398",
                "contrast": contrast,
                "manifest": manifest,
                "groups": manifest["disease_group"].map(
                    {**{value: "control" for value in controls}, **{value: "case" for value in cases}}
                ).fillna("excluded"),
            }
        )
    return result


def gene_effects_and_enrichment(repo: Path, output: Path, gmt: Path) -> tuple[pd.DataFrame, pd.DataFrame]:
    gene_rows: list[pd.DataFrame] = []
    enrichment_rows: list[pd.DataFrame] = []
    for context in contexts(repo):
        interim = repo / "data" / "interim" / str(context["interim"])
        genes = pd.read_csv(interim / "genes.csv")["gene"].astype(str)
        counts = mmread(interim / "donor_lineage_raw_counts.mtx").tocsc()
        manifest = context["manifest"].copy()
        groups = context["groups"]
        manifest["eligible_30"] = manifest["n_cells"].ge(30)
        library_sizes = np.asarray(counts.sum(axis=0)).ravel()
        log_cpm = np.log2(counts.toarray() / library_sizes * 1_000_000 + 1)
        for lineage, lineage_manifest in manifest.groupby("harmonized_lineage", sort=True):
            selected = manifest["harmonized_lineage"].eq(lineage) & manifest["eligible_30"] & groups.ne("excluded")
            selected_indices = manifest.index[selected].to_numpy()
            selected_groups = groups.loc[selected_indices]
            control_indices = selected_indices[selected_groups.eq("control").to_numpy()]
            case_indices = selected_indices[selected_groups.eq("case").to_numpy()]
            if len(control_indices) < 3 or len(case_indices) < 3:
                continue
            effects = vectorized_hedges_g(log_cpm[:, control_indices], log_cpm[:, case_indices])
            detection = np.asarray((counts[:, selected_indices] > 0).mean(axis=1)).ravel()
            frame = pd.DataFrame(
                {
                    "dataset_id": context["dataset_id"],
                    "contrast": context["contrast"],
                    "lineage": lineage,
                    "gene_symbol": genes.str.upper(),
                    "n_control": len(control_indices),
                    "n_case": len(case_indices),
                    "detection_fraction": detection,
                    "hedges_g": effects,
                }
            )
            frame = frame[frame["hedges_g"].notna() & frame["detection_fraction"].ge(MIN_GENE_DETECTION)].copy()
            frame["absolute_effect"] = frame["hedges_g"].abs()
            frame = frame.sort_values(["gene_symbol", "absolute_effect"], ascending=[True, False]).drop_duplicates("gene_symbol")
            frame = frame.drop(columns="absolute_effect")
            # Gene-level effects can be exactly tied, especially at zero.  Add a
            # fixed sub-numerical tie breaker so preranking never depends on row
            # order while leaving all reported Hedges-g values unchanged.
            tie_breaker = frame["gene_symbol"].map(
                lambda gene: int.from_bytes(hashlib.sha256(gene.encode("utf-8")).digest()[:8], "little") / 2**64
            )
            frame["prerank_score"] = frame["hedges_g"] + (tie_breaker - 0.5) * 1e-10
            frame = frame.sort_values("prerank_score", ascending=False)
            gene_rows.append(frame)
            ranking = frame[["gene_symbol", "prerank_score"]].rename(columns={"gene_symbol": "gene", "prerank_score": "score"})
            pre = gp.prerank(
                rnk=ranking,
                gene_sets=str(gmt),
                min_size=MIN_PATHWAY_SIZE,
                max_size=MAX_PATHWAY_SIZE,
                permutation_num=PERMUTATIONS,
                weight=1,
                seed=stable_seed("reactome", context["dataset_id"], context["contrast"], lineage),
                threads=4,
                outdir=None,
                no_plot=True,
                verbose=False,
            )
            enriched = pre.res2d.copy()
            rename = {
                "Term": "pathway_name",
                "ES": "enrichment_score",
                "NES": "normalized_enrichment_score",
                "NOM p-val": "nominal_p",
                "FDR q-val": "fdr_q",
                "FWER p-val": "fwer_p",
                "Tag %": "tag_fraction",
                "Gene %": "gene_fraction",
                "Lead_genes": "leading_edge_genes",
            }
            enriched = enriched.rename(columns=rename)
            if "Name" in enriched.columns:
                enriched = enriched.drop(columns="Name")
            enriched.insert(0, "lineage", lineage)
            enriched.insert(0, "contrast", context["contrast"])
            enriched.insert(0, "dataset_id", context["dataset_id"])
            enrichment_rows.append(enriched)
            print(f"{context['dataset_id']} {context['contrast']} {lineage}: {len(frame)} genes; {len(enriched)} pathways", flush=True)
    genes = pd.concat(gene_rows, ignore_index=True)
    enrichment = pd.concat(enrichment_rows, ignore_index=True)
    genes.to_csv(output / "reactome_gene_effects.csv.gz", index=False, compression="gzip")
    enrichment.to_csv(output / "reactome_preranked_enrichment.csv.gz", index=False, compression="gzip")
    return genes, enrichment


def transfer_summary(enrichment: pd.DataFrame, output: Path) -> tuple[pd.DataFrame, pd.DataFrame]:
    context_column = enrichment["dataset_id"].astype(str) + "::" + enrichment["contrast"].astype(str)
    enrichment = enrichment.copy()
    enrichment["context_id"] = context_column
    rows: list[dict[str, object]] = []
    for lineage, lineage_data in enrichment.groupby("lineage", sort=True):
        contexts = sorted(lineage_data["context_id"].unique())
        for left, right in combinations(contexts, 2):
            left_data = lineage_data[lineage_data["context_id"].eq(left)].set_index("pathway_name")
            right_data = lineage_data[lineage_data["context_id"].eq(right)].set_index("pathway_name")
            shared = left_data[["normalized_enrichment_score"]].join(
                right_data[["normalized_enrichment_score"]], how="inner", lsuffix="_left", rsuffix="_right"
            ).dropna()
            if len(shared) < 50:
                continue
            left_nes = shared["normalized_enrichment_score_left"].astype(float)
            right_nes = shared["normalized_enrichment_score_right"].astype(float)
            rows.append(
                {
                    "lineage": lineage,
                    "context_left": left,
                    "context_right": right,
                    "shared_pathways": len(shared),
                    "spearman_rho_nes": float(spearmanr(left_nes, right_nes).statistic),
                    "sign_agreement_nes": float((np.sign(left_nes) == np.sign(right_nes)).mean()),
                    "median_absolute_nes_difference": float(np.median(np.abs(left_nes - right_nes))),
                }
            )
    pairwise = pd.DataFrame(rows).sort_values(["lineage", "context_left", "context_right"])
    pairwise.to_csv(output / "reactome_pathway_transfer_pairwise.csv", index=False)

    enrichment["fdr_q"] = pd.to_numeric(enrichment["fdr_q"], errors="coerce")
    enrichment["normalized_enrichment_score"] = pd.to_numeric(enrichment["normalized_enrichment_score"], errors="coerce")
    recurrence = (
        enrichment[enrichment["fdr_q"].lt(0.05)]
        .assign(direction=lambda x: np.where(x["normalized_enrichment_score"].gt(0), "positive", "negative"))
        .groupby(["lineage", "pathway_name", "direction"], as_index=False)
        .agg(
            significant_contexts=("context_id", "nunique"),
            contexts=("context_id", lambda x: ";".join(sorted(set(x)))),
            median_nes=("normalized_enrichment_score", "median"),
            minimum_absolute_nes=("normalized_enrichment_score", lambda x: float(np.min(np.abs(x)))),
        )
        .sort_values(["significant_contexts", "minimum_absolute_nes"], ascending=[False, False])
    )
    recurrence.to_csv(output / "reactome_recurrent_pathways.csv", index=False)
    return pairwise, recurrence


def main() -> None:
    repo = Path(__file__).resolve().parents[1]
    output = repo / "results" / "phase3"
    output.mkdir(parents=True, exist_ok=True)
    gmt = repo / "data" / "external" / "reactome_2026-06" / "ReactomePathways.gmt"
    genes, enrichment = gene_effects_and_enrichment(repo, output, gmt)
    pairwise, recurrence = transfer_summary(enrichment, output)
    summary = {
        "reactome_snapshot": "official Reactome current download retrieved 2026-08-31; server files dated 2026-06",
        "reactome_gmt_sha256": sha256(gmt),
        "gene_effect_rows": len(genes),
        "enrichment_rows": len(enrichment),
        "context_lineage_analyses": int(enrichment[["dataset_id", "contrast", "lineage"]].drop_duplicates().shape[0]),
        "pairwise_transfer_rows": len(pairwise),
        "recurrent_fdr_pathways_two_or_more_contexts": int((recurrence["significant_contexts"] >= 2).sum()),
        "permutations": PERMUTATIONS,
        "gene_detection_minimum": MIN_GENE_DETECTION,
        "pathway_size_range": [MIN_PATHWAY_SIZE, MAX_PATHWAY_SIZE],
        "interpretation": "post-lock pathway convergence cannot relabel a failed published program",
    }
    (repo / "results" / "logs" / "phase3_reactome_run.json").write_text(
        json.dumps(summary, indent=2) + "\n", encoding="utf-8"
    )
    print(json.dumps(summary, indent=2))


if __name__ == "__main__":
    main()
