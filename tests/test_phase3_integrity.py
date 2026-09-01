from __future__ import annotations

from pathlib import Path

import pandas as pd


ROOT = Path(__file__).resolve().parents[1]


def test_gse256398_design_and_lineage_eligibility_are_complete() -> None:
    donors = pd.read_csv(ROOT / "metadata" / "gse256398_donor_manifest.csv")
    eligibility = pd.read_csv(
        ROOT / "metadata" / "gse256398_donor_lineage_eligibility.csv"
    )

    assert donors["donor_id"].nunique() == 26
    assert len(eligibility) == 26 * 3
    assert eligibility.groupby("donor_id")["harmonized_lineage"].nunique().eq(3).all()
    assert eligibility["n_cells"].gt(0).all()
    assert int(eligibility["eligible_30"].sum()) == 75
    assert int(eligibility["eligible_20"].sum()) == 76

    gates = pd.read_csv(ROOT / "results" / "qc" / "gse256398_donor_gate_summary.csv")
    assert len(gates) == 5 * 3 * 2
    assert gates["formal_30_cell_gate"].eq("PASS").all()
    assert gates["formal_20_cell_gate"].eq("PASS").all()

    healthy = donors.loc[donors["disease_group"].eq("healthy")]
    assert len(healthy) == 6
    assert healthy["metabolic_order"].isna().all()


def test_gse256398_program_coverage_and_effect_grid_are_locked() -> None:
    coverage = pd.read_csv(ROOT / "results" / "qc" / "gse256398_program_coverage.csv")
    effects = pd.read_csv(ROOT / "results" / "phase3" / "gse256398_program_effects.csv")
    trends = pd.read_csv(
        ROOT / "results" / "phase3" / "gse256398_metabolic_ordinal_trends.csv"
    )

    assert len(coverage) == 19
    assert coverage["coverage"].ge(0.80).all()
    assert coverage["coverage_tier"].eq("primary").all()
    assert len(effects) == 5 * 19 * 2
    assert effects["cell_gate"].eq(30).all()
    assert effects["coverage_tier"].eq("primary").all()
    assert len(trends) == 19 * 2


def test_gse256398_random_module_benchmark_is_complete() -> None:
    effects = pd.read_csv(ROOT / "results" / "phase3" / "gse256398_program_effects.csv")
    benchmark = pd.read_csv(
        ROOT / "results" / "random_controls" / "gse256398_random_module_benchmark.csv"
    )

    assert len(benchmark) == len(effects) == 190
    assert benchmark["random_modules"].eq(1000).all()
    assert int(benchmark["above_random_95th_percentile"].sum()) == 33


def test_phase3_cross_cohort_synthesis_and_meta_support_are_locked() -> None:
    matrix = pd.read_csv(
        ROOT / "results" / "phase3" / "phase3_cross_cohort_effect_matrix.csv"
    )
    summary = pd.read_csv(
        ROOT
        / "results"
        / "phase3"
        / "phase3_mash_cirrhosis_meta_program_summary.csv"
    )

    assert len(matrix) == 464
    assert int(summary["fixed_ci_positive_both_scores"].sum()) == 2
    assert int(summary["random_reml_ci_positive_both_scores"].sum()) == 0
    supported = set(
        summary.loc[summary["fixed_ci_positive_both_scores"], "program_id"]
    )
    assert supported == {"RAM2019_ENDO_2", "RAM2019_ENDO_6_SAENDO1"}
    assert not summary["formal_replication_eligible"].any()


def test_reactome_output_dimensions_and_pairwise_grid_are_locked() -> None:
    enrichment = pd.read_csv(
        ROOT / "results" / "phase3" / "reactome_preranked_enrichment.csv.gz"
    )
    gene_effects = pd.read_csv(
        ROOT / "results" / "phase3" / "reactome_gene_effects.csv.gz"
    )
    pairwise = pd.read_csv(
        ROOT / "results" / "phase3" / "reactome_pathway_transfer_pairwise.csv"
    )

    assert len(enrichment) == 18_688
    assert len(gene_effects) == 302_163
    assert len(pairwise) == 26
    assert pairwise["shared_pathways"].gt(0).all()


def test_phase3_figure_set_exists_in_both_formats() -> None:
    figures = ROOT / "results" / "figures"
    stems = [
        "figure_7_post_lock_external_enrichment",
        "supplementary_figure_1_gse256398_qc",
        "supplementary_figure_2_redundancy_weights",
        "supplementary_figure_3_histology_progression",
        "supplementary_figure_4_reactome_transfer",
        "supplementary_figure_5_threshold_precision",
        "supplementary_figure_6_bootstrap_stability",
    ]
    for stem in stems:
        for suffix in (".png", ".pdf"):
            path = figures / f"{stem}{suffix}"
            assert path.exists()
            assert path.stat().st_size > 10_000
