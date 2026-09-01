from __future__ import annotations

import hashlib
import json
from pathlib import Path

import numpy as np
import openpyxl
import pandas as pd


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


def random_matched_set(
    rng: np.random.Generator,
    measured: list[str],
    gene_index: pd.Index,
    pools: dict[str, np.ndarray],
    wider: dict[str, np.ndarray],
    fallback: np.ndarray,
) -> list[str]:
    selected: list[int] = []
    used: set[int] = set()
    for gene in measured:
        pool = pools[gene] if len(pools[gene]) > len(used) else wider[gene]
        if len(pool) <= len(used):
            pool = fallback
        choice = int(rng.choice(pool))
        while choice in used:
            choice = int(rng.choice(pool))
        selected.append(choice)
        used.add(choice)
    return gene_index[selected].tolist()


def recurrence_metrics(gene_set: set[str], ranked_lists: dict[str, list[str]]) -> tuple[int, int, float]:
    hit_count = 0
    positive_samples = 0
    rank_weight = 0.0
    for genes in ranked_lists.values():
        hits = [(index + 1, gene) for index, gene in enumerate(genes) if gene in gene_set]
        hit_count += len(hits)
        positive_samples += int(bool(hits))
        rank_weight += sum((101 - rank) / 100 for rank, _ in hits)
    return hit_count, positive_samples, rank_weight


def read_spot_sheet(path: Path, sheet: str, header_row: int) -> pd.DataFrame:
    return pd.read_excel(path, sheet_name=sheet, header=header_row, index_col=0)


def main() -> None:
    repo = Path(__file__).resolve().parents[1]
    supplement = repo / "data" / "external" / "phase5_spatial" / "PMC12162837_supplement"
    ev2 = supplement / "44321_2025_230_MOESM5_ESM.xlsx"
    ev6 = supplement / "44321_2025_230_MOESM9_ESM.xlsx"
    ev7 = supplement / "44321_2025_230_MOESM10_ESM.xlsx"
    output = repo / "results" / "phase5"
    output.mkdir(parents=True, exist_ok=True)

    # Dataset EV6 provides ranked top-100 scar-cluster genes for seven samples (A-G).
    workbook = openpyxl.load_workbook(ev6, read_only=True, data_only=True)
    sheet = workbook["Summary Scar Clusters"]
    top_rows = list(sheet.iter_rows(min_row=4, max_row=103, values_only=True))
    ranked_lists = {
        sample: [str(row[column]).upper() for row in top_rows if row[column] is not None]
        for sample, column in zip("ABCDEFG", [0, 3, 6, 9, 12, 15, 18])
    }
    workbook.close()

    # Dataset EV2 supplies a complete sample-B feature universe and cluster-level abundance
    # for expression/detection matching; it does not replace unavailable spot expression.
    background = pd.read_excel(
        ev2,
        sheet_name="Data Set EV2 - all genes",
        header=None,
        skiprows=2,
    )
    genes = background.iloc[:, 1].astype(str).str.upper()
    average_columns = [2, 5, 8, 11, 14, 17]
    averages = background.iloc[:, average_columns].apply(pd.to_numeric, errors="coerce").fillna(0)
    duplicate = genes.duplicated(keep="first")
    genes = genes.loc[~duplicate].reset_index(drop=True)
    averages = averages.loc[~duplicate].reset_index(drop=True)
    gene_index = pd.Index(genes)
    detection = averages.gt(0).sum(axis=1)
    abundance = np.log1p(averages.mean(axis=1))
    abundance_bin = pd.Series(pd.qcut(abundance, q=10, duplicates="drop", labels=False).astype(str), index=gene_index)
    detection_bin = pd.Series(detection.astype(int).astype(str).to_numpy(), index=gene_index)

    all_sets = {**PRIMARY, **SECONDARY}
    excluded = {gene for values in all_sets.values() for gene in values}
    rng = np.random.default_rng(20260901)
    recurrence_rows = []
    hit_rows = []
    for set_id, requested in all_sets.items():
        measured = [gene for gene in requested if gene in gene_index]
        observed = recurrence_metrics(set(measured), ranked_lists)
        eligible = ~gene_index.isin(excluded)
        pools = {
            gene: np.flatnonzero(
                (
                    eligible
                    & detection_bin.eq(detection_bin.loc[gene])
                    & abundance_bin.eq(abundance_bin.loc[gene])
                ).to_numpy()
            )
            for gene in measured
        }
        wider = {
            gene: np.flatnonzero(
                (eligible & detection_bin.eq(detection_bin.loc[gene])).to_numpy()
            )
            for gene in measured
        }
        fallback = np.flatnonzero(eligible)
        for sample, ranked in ranked_lists.items():
            hits = [(index + 1, gene) for index, gene in enumerate(ranked) if gene in measured]
            hit_rows.append(
                {
                    "spatial_study": "Hammond_2025_PMC12162837",
                    "gene_set_id": set_id,
                    "sample_id": sample,
                    "scar_cluster": dict(zip("ABCDEFG", ["A6", "B6", "C5", "D3", "E2", "F2", "G1"]))[sample],
                    "top100_hit_count": len(hits),
                    "hit_genes": ";".join(gene for _, gene in hits),
                    "hit_ranks": ";".join(str(rank) for rank, _ in hits),
                }
            )
        null = np.empty((ITERATIONS, 3), dtype=float)
        for iteration in range(ITERATIONS):
            random_genes = random_matched_set(
                rng,
                measured,
                gene_index,
                pools,
                wider,
                fallback,
            )
            null[iteration] = recurrence_metrics(set(random_genes), ranked_lists)
        recurrence_rows.append(
            {
                "spatial_study": "Hammond_2025_PMC12162837",
                "gene_set_id": set_id,
                "requested_genes": len(requested),
                "measured_genes": len(measured),
                "coverage_fraction": len(measured) / len(requested),
                "observed_total_top100_hits": observed[0],
                "observed_samples_with_hit": observed[1],
                "observed_rank_weight": observed[2],
                "random_total_hits_95th": float(np.quantile(null[:, 0], 0.95)),
                "random_positive_samples_95th": float(np.quantile(null[:, 1], 0.95)),
                "random_rank_weight_95th": float(np.quantile(null[:, 2], 0.95)),
                "total_hits_percentile": float((null[:, 0] <= observed[0]).mean()),
                "positive_samples_percentile": float((null[:, 1] <= observed[1]).mean()),
                "rank_weight_percentile": float((null[:, 2] <= observed[2]).mean()),
                "descriptive_recurrence_pass": bool(
                    observed[1] > np.quantile(null[:, 1], 0.95)
                    and observed[2] > np.quantile(null[:, 2], 0.95)
                ),
            }
        )

    # Dataset EV7 exposes spot-level author cell2location output. Scar-cluster labels are
    # auditable for samples A and B through EV6 (A6 and B6); sample C/D1 lacks a deposited
    # region_cluster field and is retained as not evaluable for this comparison.
    colocation_rows = []
    for donor, sheet_name, header in (
        ("A", "A1 10factors_sd8", 2),
        ("B", "B1 10factors_sd8", 0),
    ):
        spots = read_spot_sheet(ev7, sheet_name, header)
        scar = spots["region_cluster"].eq(6)
        scar_ec = spots[
            [
                "mean_spot_factors14 Endothelial_4",
                "mean_spot_factors15 Endothelial_5",
                "mean_spot_factors17 Endothelial_6",
            ]
        ].sum(axis=1)
        sinusoidal_ec = spots[
            [
                "mean_spot_factors1 Endothelial_1",
                "mean_spot_factors7 Endothelial_2",
                "mean_spot_factors8 Endothelial_3",
            ]
        ].sum(axis=1)
        fraction = scar_ec / (scar_ec + sinusoidal_ec).replace(0, np.nan)
        for metric, values in (
            ("EC4_6_abundance", scar_ec),
            ("EC1_3_abundance", sinusoidal_ec),
            ("EC4_6_fraction_of_EC1_6", fraction),
        ):
            scar_median = float(values.loc[scar].median())
            parenchyma_median = float(values.loc[~scar].median())
            colocation_rows.append(
                {
                    "spatial_study": "Hammond_2025_PMC12162837",
                    "donor_id": donor,
                    "scar_cluster": f"{donor}6",
                    "metric": metric,
                    "scar_spots": int(scar.sum()),
                    "other_spots": int((~scar).sum()),
                    "scar_median": scar_median,
                    "other_median": parenchyma_median,
                    "scar_minus_other": scar_median - parenchyma_median,
                }
            )
    colocation_rows.append(
        {
            "spatial_study": "Hammond_2025_PMC12162837",
            "donor_id": "C",
            "scar_cluster": "C5",
            "metric": "NOT_EVALUABLE",
            "scar_spots": np.nan,
            "other_spots": np.nan,
            "scar_median": np.nan,
            "other_median": np.nan,
            "scar_minus_other": np.nan,
        }
    )

    recurrence = pd.DataFrame(recurrence_rows)
    hits = pd.DataFrame(hit_rows)
    colocation = pd.DataFrame(colocation_rows)
    recurrence.to_csv(output / "spatial_2025_scar_top100_random_benchmark.csv", index=False)
    hits.to_csv(output / "spatial_2025_scar_top100_hits.csv", index=False)
    colocation.to_csv(output / "spatial_2025_endothelial_colocation.csv", index=False)
    summary = {
        "source_files": {
            path.name: sha256(path) for path in (ev2, ev6, ev7)
        },
        "experimental_independent_donors": 3,
        "validation_independent_donors": 4,
        "scar_top100_samples": 7,
        "spot_level_deconvolution_samples": 3,
        "spot_region_label_evaluable_samples": 2,
        "raw_fastq_size_gib": {
            "E-MTAB-13132": 92.72457328997552,
            "E-MTAB-14960": 25.156081312336028,
        },
        "strict_frozen_status": "SPATIAL_RESOURCE_PARTIALLY_EVALUABLE",
        "strict_reason": (
            "Public supplements provide scar-cluster top genes and spot-level cell2location output, "
            "but not a spot-by-gene expression matrix with complete scar labels for all donors."
        ),
        "primary_descriptive_recurrence": recurrence.loc[
            recurrence["gene_set_id"].isin(PRIMARY),
            ["gene_set_id", "descriptive_recurrence_pass"],
        ].set_index("gene_set_id")["descriptive_recurrence_pass"].to_dict(),
        "spatially_replicated_post_lock": False,
    }
    (output / "spatial_2025_run_summary.json").write_text(
        json.dumps(summary, indent=2) + "\n", encoding="utf-8"
    )
    print(recurrence.to_string(index=False))
    print(colocation.to_string(index=False))


if __name__ == "__main__":
    main()
