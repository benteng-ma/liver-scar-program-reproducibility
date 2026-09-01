from __future__ import annotations

import json
import subprocess
import unittest
from pathlib import Path

import pandas as pd


ROOT = Path(__file__).resolve().parents[1]


class Phase0IntegrityTests(unittest.TestCase):
    @classmethod
    def setUpClass(cls):
        cls.datasets = pd.read_csv(ROOT / "metadata" / "dataset_manifest.csv")
        cls.donors = pd.read_csv(ROOT / "metadata" / "donor_manifest.csv")
        cls.programs = pd.read_csv(ROOT / "literature" / "program_inventory.csv")
        cls.coverage = pd.read_csv(ROOT / "metadata" / "program_assay_coverage.csv")

    def test_project_isolation(self):
        actual = subprocess.check_output(
            ["git", "rev-parse", "--show-toplevel"], cwd=ROOT, text=True
        ).strip()
        self.assertEqual(Path(actual).resolve(), ROOT.resolve())

    def test_donor_ids_are_nonempty(self):
        self.assertFalse(self.donors["donor_id"].isna().any())
        self.assertFalse(self.donors["sample_id"].isna().any())

    def test_technical_fractions_are_not_donors(self):
        gse = self.donors[self.donors.dataset_id.eq("GSE136103")]
        self.assertEqual(len(gse), 20)
        self.assertEqual(gse.donor_id.nunique(), 10)

    def test_regions_are_not_donors(self):
        gse = self.donors[self.donors.dataset_id.eq("GSE202379")]
        self.assertEqual(len(gse), 59)
        self.assertEqual(gse.donor_id.nunique(), 47)

    def test_e_mtab_is_normal_and_reuse_flagged(self):
        gse = self.donors[self.donors.dataset_id.eq("E-MTAB-10553_unique")]
        self.assertEqual(gse.donor_id.nunique(), 6)
        self.assertEqual(set(gse.disease_group.str.lower()), {"normal"})
        note = self.datasets.loc[
            self.datasets.dataset_id.eq("E-MTAB-10553_unique"), "limitations"
        ].iloc[0]
        self.assertIn("GSE136103", note)

    def test_human_gse181483_is_four_donors(self):
        gse = self.donors[self.donors.dataset_id.eq("GSE181483_human")]
        self.assertEqual(gse.donor_id.nunique(), 4)
        self.assertEqual(gse.disease_group.value_counts().to_dict(), {"healthy": 2, "cirrhosis": 2})

    def test_programs_are_traceable_and_directional(self):
        self.assertEqual(self.programs.program_id.nunique(), 19)
        self.assertGreaterEqual(self.programs.program_id.nunique(), 15)
        self.assertFalse(self.programs.direction.isna().any())
        self.assertFalse(self.programs.source_table.isna().any())
        self.assertTrue((self.programs.verified.astype(str).str.upper() == "TRUE").all())
        self.assertEqual(self.programs.duplicated(["program_id", "gene_symbol"]).sum(), 0)

    def test_discovery_data_cannot_be_external_validation(self):
        self.assertEqual(set(self.programs.discovery_dataset), {"GSE136103"})
        self.assertEqual(set(self.programs.data_lineage_class), {"REUSED"})

    def test_missing_genes_are_not_imputed_in_coverage(self):
        calculated = self.coverage.n_detected / self.coverage.n_program_genes
        self.assertTrue((calculated.round(6) == self.coverage.coverage.round(6)).all())
        self.assertTrue((self.coverage.n_detected <= self.coverage.n_program_genes).all())

    def test_spatial_failure_is_coverage_not_negativity(self):
        spatial = self.coverage[self.coverage.assay_resource.eq("MERFISH_GSE210077_panel")]
        self.assertEqual((spatial.coverage_tier == "spatial_evaluable").sum(), 0)
        decision = json.loads((ROOT / "reports" / "phase0_decision.json").read_text(encoding="utf-8"))
        self.assertEqual(decision["decision"], "CONDITIONAL_GO_NO_SPATIAL")
        self.assertFalse(decision["spatial_main_analysis"])

    def test_smoke_tests_did_not_run_disease_effects(self):
        result = json.loads((ROOT / "results" / "qc" / "phase0_smoke_test.json").read_text(encoding="utf-8"))
        self.assertIn("No disease effect", result["prohibition"])
        self.assertEqual(result["coverage_summary"]["scRNA_primary_coverage_programs"], 19)
        self.assertEqual(result["coverage_summary"]["snRNA_primary_coverage_programs"], 19)


if __name__ == "__main__":
    unittest.main()
