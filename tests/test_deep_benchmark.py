from __future__ import annotations

import csv
import unittest
from pathlib import Path


ROOT = Path(__file__).resolve().parents[1]
RESULTS = ROOT / "results" / "deep_benchmark"


def read_rows(name: str) -> list[dict[str, str]]:
    with (RESULTS / name).open(encoding="utf-8-sig", newline="") as handle:
        return list(csv.DictReader(handle))


class DeepTransportabilityBenchmarkTests(unittest.TestCase):
    def test_amendment_and_control_configuration_exist(self) -> None:
        self.assertTrue((ROOT / "reports" / "deep_transportability_benchmark_plan_2026-08-31.md").is_file())
        controls = ROOT / "config" / "deep_benchmark_control_programs.csv"
        self.assertTrue(controls.is_file())
        with controls.open(encoding="utf-8-sig", newline="") as handle:
            rows = list(csv.DictReader(handle))
        self.assertGreaterEqual(len(rows), 40)
        self.assertEqual({row["control_class"] for row in rows}, {"lineage_identity", "disease_response"})

    def test_noncircular_identity_calibration_is_frozen(self) -> None:
        rows = read_rows("identity_control_performance.csv")
        author = [row for row in rows if row["annotation"] == "author"]
        self.assertEqual(len(author), 4)
        gse = [row for row in author if row["dataset_id"] == "GSE202379"]
        watson = [row for row in author if row["dataset_id"] == "GSE210077_Watson6"]
        self.assertEqual(len(gse), 2)
        self.assertEqual(len(watson), 2)
        field = "passes_frozen_positive_control_threshold"
        self.assertTrue(all(row[field] == "True" for row in gse))
        self.assertTrue(all(row[field] == "False" for row in watson))
        self.assertTrue(all(float(row["top_score_lineage_accuracy"]) == 1.0 for row in gse))
        self.assertTrue(all(abs(float(row["top_score_lineage_accuracy"]) - 2 / 3) < 1e-12 for row in watson))

    def test_endpoint_alignment_result_is_frozen(self) -> None:
        rows = {row["stratum"]: row for row in read_rows("transfer_failure_stratified_summary.csv")}
        all_pairs = rows["all_pairs"]
        matched = rows["comparable_advanced_endpoint"]
        self.assertAlmostEqual(float(all_pairs["program_pair_sign_agreement"]), 0.5)
        self.assertAlmostEqual(float(all_pairs["median_pairwise_program_spearman"]), 0.4464285714285715)
        self.assertAlmostEqual(float(matched["program_pair_sign_agreement"]), 11 / 12)
        self.assertAlmostEqual(float(matched["median_pairwise_program_spearman"]), 0.6857142857142857)

    def test_held_out_magnitude_failure_is_frozen(self) -> None:
        all_rows = read_rows("leave_one_cohort_out_prediction_summary.csv")
        overall = next(row for row in all_rows if row["held_out_dataset"] == "ALL")
        self.assertEqual(int(overall["n_program_predictions"]), 152)
        self.assertAlmostEqual(float(overall["spearman_rho"]), 0.3635139660604569)
        self.assertAlmostEqual(float(overall["predictive_r_squared"]), 0.024403662615125388)
        no_watson = read_rows("leave_one_cohort_out_prediction_summary_excluding_watson.csv")
        overall_no_watson = next(row for row in no_watson if row["held_out_dataset"] == "ALL_NON_WATSON")
        self.assertEqual(int(overall_no_watson["n_program_predictions"]), 90)
        self.assertAlmostEqual(float(overall_no_watson["predictive_r_squared"]), 0.009989907026054734)

    def test_core5_membership_and_failure_are_frozen(self) -> None:
        membership = read_rows("minimal_core_membership.csv")
        by_lineage: dict[str, set[str]] = {}
        for row in membership:
            by_lineage.setdefault(row["lineage"], set()).add(row["gene_symbol"])
        self.assertEqual(by_lineage["endothelial"], {"TFF3", "PLPP1", "FTL", "CPE", "FTH1"})
        self.assertEqual(
            by_lineage["mesenchymal_hsc_myofibroblast"],
            {"MDK", "CST3", "TM4SF1", "SERPINF1", "TMSB10"},
        )
        self.assertNotIn("macrophage_monocyte", by_lineage)
        summary = read_rows("minimal_core_validation_summary.csv")
        self.assertTrue(all(row["held_out_directional_core"] == "False" for row in summary))

    def test_report_card_is_prioritization_not_replication(self) -> None:
        rows = read_rows("program_transportability_report_card.csv")
        self.assertEqual(len(rows), 19)
        ordered = sorted(rows, key=lambda row: float(row["transportability_readiness_total_0_100"]), reverse=True)
        self.assertEqual(ordered[0]["program_id"], "RAM2019_MAC_SIG_B_SAM")
        self.assertEqual(ordered[1]["program_id"], "RAM2019_MAC_SIG_A_SAM")
        self.assertTrue(all(float(row["transportability_readiness_total_0_100"]) <= 100 for row in rows))
        self.assertFalse(any(float(row["endpoint_evidence_domain_0_20"]) == 20 for row in rows))


if __name__ == "__main__":
    unittest.main()
