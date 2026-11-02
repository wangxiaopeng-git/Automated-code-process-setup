import json
import tempfile
import unittest
from pathlib import Path


class AdvancedCliTests(unittest.TestCase):
    def test_design_and_reference_commands(self):
        from batch_primer_design.advanced_cli import main

        with tempfile.TemporaryDirectory() as d:
            d = Path(d)
            genome = d / "g.fa"
            genome.write_text(">c\n" + "ACGT" * 300 + "\n")
            self.assertEqual(
                main(
                    [
                        "reference",
                        "import",
                        "--genome",
                        str(genome),
                        "--taxid",
                        "123",
                        "--assembly",
                        "local-test",
                        "-o",
                        str(d / "ref"),
                    ]
                ),
                0,
            )
            self.assertEqual(main(["reference", "verify", str(d / "ref")]), 0)
            bed = d / "a.bed"
            bed.write_text("c\t490\t500\ttarget\n")
            code = main(
                [
                    "design",
                    str(bed),
                    "--reference",
                    str(d / "ref"),
                    "--specificity",
                    "exact",
                    "--flank",
                    "200",
                    "-o",
                    str(d / "out"),
                ]
            )
            self.assertIn(code, [0, 1])
            self.assertTrue((d / "out" / "results.json").exists())

    def test_invalid_config_does_not_write_output(self):
        from batch_primer_design.advanced_cli import main

        with tempfile.TemporaryDirectory() as d:
            d = Path(d)
            p = d / "a.fa"
            p.write_text(">a\nACGT\n")
            self.assertEqual(
                main(
                    [
                        "design",
                        str(p),
                        "--config",
                        str(d / "missing.json"),
                        "-o",
                        str(d / "out"),
                    ]
                ),
                2,
            )
            self.assertFalse((d / "out").exists())
