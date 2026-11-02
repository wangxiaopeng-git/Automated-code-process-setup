import json
import random
import tempfile
import unittest
from pathlib import Path
from batch_primer_design.core import DesignConfig, TemplateRecord


class WorkflowTests(unittest.TestCase):
    def test_pairing_includes_reverse_amplicon(self):
        from batch_primer_design.specificity import Hit, pair_hits

        hits = [
            Hit("L", "chr", 10, 30, "+"),
            Hit("R", "chr", 110, 130, "-"),
            Hit("R", "chr", 210, 230, "+"),
            Hit("L", "chr", 310, 330, "-"),
        ]
        products = pair_hits(hits, 100, 150)
        self.assertEqual(
            [(p["start"], p["end"], p["strand"]) for p in products],
            [(10, 130, "+"), (210, 330, "-")],
        )

    def test_wrong_orientation_does_not_amplify(self):
        from batch_primer_design.specificity import Hit, pair_hits

        self.assertEqual(
            pair_hits(
                [Hit("L", "c", 10, 30, "-"), Hit("R", "c", 110, 130, "+")], 100, 150
            ),
            [],
        )

    def test_exact_scan_repeated_offtargets(self):
        from batch_primer_design.specificity import exact_hits, pair_hits

        hits = exact_hits({"c": "AAAACCCCGGGGTTTTAAAACCCCGGGGTTTT"}, "AAAA", "CCCC")
        products = pair_hits(hits, 8, 16)
        self.assertGreaterEqual(len(products), 2)

    def test_real_design_and_exports(self):
        from batch_primer_design.workflow import run_design
        from batch_primer_design.reference import Target

        random.seed(41)
        seq = "".join(random.choices("ACGT", k=600))
        with tempfile.TemporaryDirectory() as d:
            root = Path(d) / "results"
            report = run_design(
                [Target(TemplateRecord("a", seq, 280, 5))],
                DesignConfig(num_return=2),
                root,
            )
            self.assertEqual(report["summary"]["successful_records"], 1)
            self.assertEqual(report["results"][0]["specificity_status"], "not_run")
            for name in (
                "results.csv",
                "results.tsv",
                "results.xlsx",
                "results.json",
                "primers.fasta",
                "amplicons.fasta",
                "primers.bed",
                "report.html",
                "manifest.json",
            ):
                self.assertTrue((root / name).is_file(), name)
            row = report["results"][0]
            a = row["left_binding_start_1based"]
            b = row["right_binding_end_1based"]
            self.assertIn(seq[a - 1 : b], (root / "amplicons.fasta").read_text())
            with self.assertRaises(Exception):
                run_design([Target(TemplateRecord("a", seq))], DesignConfig(), root)

    def test_duplicate_records_counted_once_per_input(self):
        from batch_primer_design.workflow import run_design
        from batch_primer_design.reference import Target

        with tempfile.TemporaryDirectory() as d:
            r = run_design(
                [
                    Target(TemplateRecord("a", "ACGT" * 100)),
                    Target(TemplateRecord("a", "ACGT" * 100)),
                ],
                DesignConfig(),
                Path(d) / "r",
            )
            self.assertEqual(r["summary"]["input_records"], 2)
            self.assertEqual(
                r["summary"]["successful_records"] + r["summary"]["failed_records"], 2
            )
            self.assertEqual(r["results"][-1]["status"], "error")
