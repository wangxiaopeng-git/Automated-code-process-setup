import hashlib
import io
import shutil
import tempfile
import unittest
import zipfile
from pathlib import Path
from unittest.mock import patch
from types import SimpleNamespace
from batch_primer_design.core import InputError
from batch_primer_design.reference import Target
from batch_primer_design.core import TemplateRecord


class ReviewRegressionTests(unittest.TestCase):
    def test_unequal_primers_short_offtarget(self):
        from batch_primer_design.specificity import check_pair, Hit

        row = {
            "left_sequence_5to3": "A" * 20,
            "right_sequence_5to3": "C" * 25,
            "left_binding_start_1based": 1,
            "right_binding_end_1based": 100,
        }
        target = Target(TemplateRecord("a", "A" * 100), "c", 0, 100, "+", "test")
        bundle = SimpleNamespace(
            manifest={"contigs": {"c": 400}, "assembly": "test"}, sequences={}
        )
        hits = [
            Hit("L", "c", 0, 20, "+"),
            Hit("R", "c", 75, 100, "-"),
            Hit("L", "c", 300, 320, "+"),
            Hit("R", "c", 320, 345, "-"),
        ]
        with patch("batch_primer_design.specificity.exact_hits", return_value=hits):
            r = check_pair(row, target, bundle, "exact", max_product=200)
        self.assertEqual(r["specificity_status"], "multiple_products")
        self.assertEqual([p["size_bp"] for p in r["predicted_products"]], [100, 45])

    def test_empty_md5_manifest_not_verified(self):
        from batch_primer_design.ncbi import safe_extract

        with tempfile.TemporaryDirectory() as d:
            p = Path(d) / "a.zip"
            with zipfile.ZipFile(p, "w") as z:
                z.writestr("md5sum.txt", "")
                z.writestr("genomic.fna", ">a\nACGT\n")
            self.assertFalse(safe_extract(p, Path(d) / "out", 10000))

    def test_partial_md5_manifest_not_verified(self):
        from batch_primer_design.ncbi import safe_extract

        with tempfile.TemporaryDirectory() as d:
            p = Path(d) / "a.zip"
            with zipfile.ZipFile(p, "w") as z:
                z.writestr(
                    "md5sum.txt", hashlib.md5(b"info").hexdigest() + "  info.txt\n"
                )
                z.writestr("info.txt", "info")
                z.writestr("genomic.fna", ">a\nACGT\n")
            self.assertFalse(safe_extract(p, Path(d) / "out", 10000))

    @unittest.skipUnless(shutil.which("blastn"), "BLAST+ not installed")
    def test_real_blast_recovers_terminal_and_two_mismatch_hits(self):
        from batch_primer_design.specificity import blast_hits

        primer = "ACGTCAGTGCATGACTGA"
        terminal = "T" + primer[1:]
        two = list(primer)
        two[5] = "T"
        two[11] = "C"
        two = "".join(two)
        with tempfile.TemporaryDirectory() as d:
            p = Path(d) / "g.fa"
            p.write_text(
                ">terminal\n"
                + "G" * 30
                + terminal
                + "G" * 30
                + "\n>two\n"
                + "G" * 30
                + two
                + "G" * 30
                + "\n"
            )
            hits = blast_hits(p, primer, "TTGCCAGTACGATCGTAC", max_mismatches=2)
            matching = [
                h for h in hits if h.primer == "L" and h.start == 30 and h.strand == "+"
            ]
            self.assertEqual({h.contig for h in matching}, {"terminal", "two"})
            self.assertEqual(
                next(h.mismatches for h in matching if h.contig == "terminal"), (1,)
            )
