from __future__ import annotations

import json
import unittest
from pathlib import Path

import pandas as pd


ROOT = Path(__file__).resolve().parents[1]


class Phase2IntegrityTests(unittest.TestCase):
    def test_phase2_environment_verification_passed(self) -> None:
        report = json.loads(
            (ROOT / "results" / "logs" / "phase2_environment_verification.json").read_text()
        )
        self.assertIs(report["verified"], True)
        self.assertIs(report["isolated_from_system_site_packages"], True)

    def test_program_freeze_checksum_is_unchanged(self) -> None:
        freeze = (ROOT / "metadata" / "program_freeze.sha256").read_text()
        expected = "CF05D2E589EAB4C5BCAEC4005D868D3E3EB3F9F7E5E7A7D15D8365C2EE88F386"
        self.assertIn(expected, freeze)

    def test_gse202379_primary_rows_respect_gate_and_coverage(self) -> None:
        effects = pd.read_csv(ROOT / "results" / "primary" / "gse202379_primary_effects.csv")
        gates = (
            pd.read_csv(ROOT / "results" / "qc" / "gse202379_donor_gate_summary.csv")
            .drop_duplicates(["contrast", "lineage"])
            .set_index(["contrast", "lineage"])
        )
        coverage = pd.read_csv(ROOT / "results" / "qc" / "gse202379_program_coverage.csv").set_index("program_id")
        for row in effects.itertuples(index=False):
            self.assertEqual(gates.loc[(row.contrast, row.lineage), "formal_primary_gate"], "PASS")
            self.assertEqual(coverage.loc[row.program_id, "coverage_tier"], "primary")
            self.assertEqual(row.cell_gate, 30)
            self.assertGreaterEqual(row.n_control, 3)
            self.assertGreaterEqual(row.n_case, 3)

    def test_gse202379_random_benchmark_is_complete(self) -> None:
        effects = pd.concat(
            [
                pd.read_csv(ROOT / "results" / "primary" / "gse202379_primary_effects.csv"),
                pd.read_csv(ROOT / "results" / "sensitivity" / "gse202379_sensitivity_effects.csv"),
            ],
            ignore_index=True,
        )
        benchmark = pd.read_csv(ROOT / "results" / "random_controls" / "gse202379_random_module_benchmark.csv")
        self.assertEqual(len(benchmark), len(effects))
        self.assertTrue(benchmark["random_modules"].eq(1000).all())

    def test_gse290642_mapping_does_not_override_automatic_rule(self) -> None:
        mapping = pd.read_csv(ROOT / "metadata" / "gse290642_cluster_mapping.csv").set_index("cluster")
        audit = pd.read_csv(ROOT / "results" / "qc" / "gse290642_cluster_marker_audit.csv").set_index("cluster")
        included = mapping[mapping["target_included"].eq("yes")]
        self.assertEqual(len(mapping), 40)
        self.assertTrue(audit.loc[included.index, "passes_automatic_rule"].all())
        self.assertTrue(
            (
                audit.loc[included.index, "provisional_label"]
                == included["harmonized_lineage"]
            ).all()
        )

    def test_gse290642_uses_only_shared_features(self) -> None:
        summary = json.loads(
            (ROOT / "results" / "logs" / "gse290642_analysis_run.json").read_text()
        )
        gate_summary = json.loads(
            (ROOT / "results" / "qc" / "gse290642_gate_audit_summary.json").read_text()
        )
        n_shared = gate_summary["genes_shared_all_24_donors"]
        self.assertEqual(summary["feature_space"], f"{n_shared} genes measured in all 24 donors")
        self.assertLess(n_shared, gate_summary["genes"])

    def test_gse290642_effects_are_gate_eligible_and_sensitivity_only(self) -> None:
        effects = pd.read_csv(ROOT / "results" / "sensitivity" / "gse290642_sensitivity_effects.csv")
        gates = (
            pd.read_csv(ROOT / "results" / "qc" / "gse290642_donor_gate_summary.csv")
            .drop_duplicates(["contrast", "lineage"])
            .set_index(["contrast", "lineage"])
        )
        self.assertTrue(effects["analysis_tier"].eq("sensitivity").all())
        for row in effects.itertuples(index=False):
            column = "formal_30_cell_gate" if row.cell_gate == 30 else "formal_20_cell_gate"
            self.assertEqual(gates.loc[(row.contrast, row.lineage), column], "PASS")
            self.assertGreaterEqual(row.n_control, 3)
            self.assertGreaterEqual(row.n_case, 3)

    def test_gse290642_random_modules_match_every_effect(self) -> None:
        effects = pd.read_csv(ROOT / "results" / "sensitivity" / "gse290642_sensitivity_effects.csv")
        benchmark = pd.read_csv(ROOT / "results" / "random_controls" / "gse290642_random_module_benchmark.csv")
        matching = pd.read_csv(ROOT / "results" / "random_controls" / "gse290642_random_module_matching_qc.csv")
        self.assertEqual(len(effects), len(benchmark))
        self.assertTrue(benchmark["random_modules"].eq(1000).all())
        self.assertTrue(matching["exact_bin_match_fraction"].eq(1.0).all())
        self.assertTrue(matching["within_module_replacement"].eq(False).all())

    def test_gse244832_is_mixed_stage_sensitivity_only(self) -> None:
        effects = pd.read_csv(ROOT / "results" / "sensitivity" / "gse244832_sensitivity_effects.csv")
        gates = pd.read_csv(ROOT / "results" / "qc" / "gse244832_donor_gate_summary.csv")
        benchmark = pd.read_csv(ROOT / "results" / "random_controls" / "gse244832_random_module_benchmark.csv")
        self.assertTrue(effects["analysis_tier"].eq("sensitivity").all())
        self.assertTrue(effects["endpoint_limitation"].str.contains("F2-F4").all())
        self.assertTrue(gates["formal_30_cell_gate"].eq("PASS").all())
        self.assertEqual(len(effects), len(benchmark))
        self.assertTrue(benchmark["random_modules"].eq(1000).all())

    def test_watson_uses_shared_features_and_has_no_expected_random_exceedance(self) -> None:
        gate = json.loads(
            (ROOT / "results" / "qc" / "gse210077_watson6_gate_audit_summary.json").read_text()
        )
        analysis = json.loads(
            (ROOT / "results" / "logs" / "gse210077_watson6_analysis_run.json").read_text()
        )
        effects = pd.read_csv(
            ROOT / "results" / "sensitivity" / "gse210077_watson6_sensitivity_effects.csv"
        )
        benchmark = pd.read_csv(
            ROOT / "results" / "random_controls" / "gse210077_watson6_random_module_benchmark.csv"
        )
        self.assertEqual(
            analysis["feature_space"],
            f"{gate['genes_shared_all_six_donors']} genes measured in all six donors",
        )
        self.assertTrue(effects["hedges_g"].lt(0).all())
        self.assertFalse(benchmark["above_random_95th_percentile"].any())

    def test_gse181483_mapping_and_directional_boundary(self) -> None:
        mapping = pd.read_csv(ROOT / "metadata" / "gse181483_cluster_mapping.csv").set_index("cluster")
        audit = pd.read_csv(ROOT / "results" / "qc" / "gse181483_cluster_marker_audit.csv").set_index("cluster")
        included = mapping[mapping["target_included"].eq("yes")]
        self.assertEqual(len(mapping), 25)
        self.assertTrue(audit.loc[included.index, "passes_automatic_rule"].all())
        self.assertTrue(
            (audit.loc[included.index, "provisional_label"] == included["harmonized_lineage"]).all()
        )
        effects = pd.read_csv(ROOT / "results" / "sensitivity" / "gse181483_directional_effects.csv")
        gates = pd.read_csv(ROOT / "results" / "qc" / "gse181483_donor_gate_summary.csv")
        self.assertTrue(effects["analysis_tier"].eq("directional").all())
        self.assertNotIn("permutation_p_two_sided", effects.columns)
        self.assertNotIn("robust_ci95_low", effects.columns)
        self.assertTrue(gates["formal_three_donor_gate"].eq("FAIL").all())

    def test_cross_cohort_labels_do_not_overreach(self) -> None:
        classifications = pd.read_csv(ROOT / "results" / "meta" / "program_classification_table.csv")
        meta = pd.read_csv(ROOT / "results" / "meta" / "advanced_endothelial_sensitivity_meta.csv")
        self.assertEqual(len(classifications), 19)
        self.assertFalse(classifications["within_cell_state_replicated"].any())
        self.assertFalse(classifications["pan_cirrhotic_transportable"].any())
        self.assertTrue(classifications["assay_transfer_class"].eq("UNRESOLVED").all())
        self.assertFalse(meta["formal_replication_eligible"].any())
        self.assertTrue(meta["k"].eq(2).all())

    def test_validation_random_benchmarks_have_one_thousand_modules(self) -> None:
        for prefix, effect_file in (
            ("gse181483", "gse181483_directional_effects.csv"),
            ("gse244832", "gse244832_sensitivity_effects.csv"),
            ("gse210077_watson6", "gse210077_watson6_sensitivity_effects.csv"),
        ):
            effects = pd.read_csv(ROOT / "results" / "sensitivity" / effect_file)
            benchmark = pd.read_csv(
                ROOT / "results" / "random_controls" / f"{prefix}_random_module_benchmark.csv"
            )
            self.assertEqual(len(effects), len(benchmark))
            self.assertTrue(benchmark["random_modules"].eq(1000).all())


if __name__ == "__main__":
    unittest.main()
