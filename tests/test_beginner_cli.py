"""User workflows: use local files directly, no manual bundle setup."""

import contextlib
import io
import json
import os
import random
import shutil
import tempfile
import unittest
from pathlib import Path
from unittest.mock import patch
from batch_primer_design.advanced_cli import main


class BeginnerCliTests(unittest.TestCase):
    def setUp(self):
        self.tmp = tempfile.TemporaryDirectory()
        self.addCleanup(self.tmp.cleanup)
        self.root = Path(self.tmp.name)
        self.project = self.root / "data" / "project1"
        self.project.mkdir(parents=True)
        rng = random.Random(42)
        seq = "".join(rng.choices("ACGT", k=1200))
        self.genome = self.project / "genome.fa"
        self.genome.write_text(">chr1\n" + seq + "\n")
        self.fasta = self.project / "targets.fa"
        self.fasta.write_text(">target\n" + seq[300:900] + "\n")
        self.gff = self.project / "genes.gff3"
        self.gff.write_text("chr1\tx\tgene\t601\t610\t.\t+\t.\tID=gene1\n")
        self.gtf = self.project / "genes.gtf"
        self.gtf.write_text(
            'chr1\tx\texon\t601\t610\t.\t-\t.\tgene_id "gene1"; transcript_id "tx1";\n'
        )
        self.ids = self.project / "genes.txt"
        self.ids.write_text("gene1\n")

    def run_cli(self, args):
        with contextlib.redirect_stdout(io.StringIO()), contextlib.redirect_stderr(
            io.StringIO()
        ) as errors:
            result = main(args)
            self.last_error = errors.getvalue()
            return result

    def test_no_reference_defaults_to_project_results_and_no_check(self):
        old = Path.cwd()
        try:
            os.chdir(self.root)
            self.assertEqual(self.run_cli(["design", "data/project1/targets.fa"]), 0)
            report = json.loads(Path("results/project1/results.json").read_text())
            self.assertEqual(report["results"][0]["specificity_status"], "not_run")
            self.assertIsNone(report["reference"])
        finally:
            os.chdir(old)

    def test_raw_fasta_and_gff_auto_import_and_reuse_cache(self):
        for n in range(2):
            out = self.root / f"results{n}"
            self.assertEqual(
                self.run_cli(
                    [
                        "design",
                        str(self.ids),
                        "--reference",
                        str(self.genome),
                        "--gff",
                        str(self.gff),
                        "--specificity",
                        "exact",
                        "-o",
                        str(out),
                    ]
                ),
                0,
            )
            report = json.loads((out / "results.json").read_text())
            self.assertEqual(
                report["results"][0]["specificity_status"], "single_product_exact"
            )
            self.assertIsNone(report["reference"]["taxid"])
        manifests = list((self.project / ".primer-cache").rglob("manifest.json"))
        self.assertEqual(len(manifests), 1)
        self.assertEqual(self.genome.read_text().splitlines()[0], ">chr1")

    def test_gtf_alias_resolves_negative_strand(self):
        out = self.root / "out"
        self.assertEqual(
            self.run_cli(
                [
                    "design",
                    str(self.ids),
                    "--reference",
                    str(self.genome),
                    "--gtf",
                    str(self.gtf),
                    "--specificity",
                    "exact",
                    "-o",
                    str(out),
                ]
            ),
            0,
        )
        report = json.loads((out / "results.json").read_text())
        self.assertEqual(report["results"][0]["template_strand"], "-")

    @unittest.skipUnless(shutil.which("blastn"), "BLAST+ not installed")
    def test_reference_enables_blast_without_extra_flag(self):
        out = self.root / "out"
        self.assertEqual(
            self.run_cli(
                [
                    "design",
                    str(self.fasta),
                    "--reference",
                    str(self.genome),
                    "-o",
                    str(out),
                ]
            ),
            0,
        )
        report = json.loads((out / "results.json").read_text())
        self.assertEqual(report["results"][0]["specificity_method"], "blast")
        self.assertTrue(report["results"][0]["predicted_products"])

    def test_annotation_requires_reference(self):
        out = self.root / "out"
        self.assertEqual(
            self.run_cli(
                ["design", str(self.fasta), "--gtf", str(self.gtf), "-o", str(out)]
            ),
            2,
        )
        self.assertFalse(out.exists())

    def test_missing_blast_is_actionable_and_does_not_make_output(self):
        out = self.root / "out"
        with patch("shutil.which", return_value=None):
            self.assertEqual(
                self.run_cli(
                    [
                        "design",
                        str(self.fasta),
                        "--reference",
                        str(self.genome),
                        "-o",
                        str(out),
                    ]
                ),
                2,
            )
        self.assertFalse(out.exists())

    def test_changed_reference_content_gets_new_cache(self):
        for n in range(2):
            if n:
                self.genome.write_text(self.genome.read_text() + "A\n")
            self.assertEqual(
                self.run_cli(
                    [
                        "design",
                        str(self.fasta),
                        "--reference",
                        str(self.genome),
                        "--specificity",
                        "exact",
                        "-o",
                        str(self.root / f"out{n}"),
                    ]
                ),
                0,
            )
        self.assertEqual(
            len(list((self.project / ".primer-cache").glob("*/manifest.json"))), 2
        )

    def test_annotation_without_reference_does_not_create_cache(self):
        self.assertEqual(
            self.run_cli(
                [
                    "design",
                    str(self.fasta),
                    "--gff",
                    str(self.gff),
                    "-o",
                    str(self.root / "out"),
                ]
            ),
            2,
        )
        self.assertFalse((self.project / ".primer-cache").exists())

    def test_explicit_none_skips_blast_requirement(self):
        with patch("shutil.which", return_value=None):
            out = self.root / "out"
            self.assertEqual(
                self.run_cli(
                    [
                        "design",
                        str(self.fasta),
                        "--reference",
                        str(self.genome),
                        "--specificity",
                        "none",
                        "-o",
                        str(out),
                    ]
                ),
                0,
            )
        r = json.loads((out / "results.json").read_text())
        self.assertEqual(r["results"][0]["specificity_status"], "not_run")

    def test_gtf_and_gff_cannot_be_supplied_together(self):
        with self.assertRaises(SystemExit) as caught:
            self.run_cli(
                [
                    "design",
                    str(self.fasta),
                    "--reference",
                    str(self.genome),
                    "--gtf",
                    str(self.gtf),
                    "--gff",
                    str(self.gff),
                ]
            )
        self.assertEqual(caught.exception.code, 2)

    def test_bare_input_uses_current_project_folder_name(self):
        old = Path.cwd()
        try:
            os.chdir(self.project)
            self.assertEqual(self.run_cli(["design", "targets.fa"]), 0)
            self.assertTrue(Path("results/project1/results.json").is_file())
        finally:
            os.chdir(old)
