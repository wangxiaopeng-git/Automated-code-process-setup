import random
import shutil
import tempfile
import unittest
from pathlib import Path
from batch_primer_design.core import DesignConfig
from batch_primer_design.reference import ReferenceBundle, import_reference
from batch_primer_design.workflow import run_design
from Bio.Seq import Seq


class ReferenceIntegrationTests(unittest.TestCase):
    def exercise(self, method):
        with tempfile.TemporaryDirectory() as d:
            d = Path(d)
            rng = random.Random(123)
            sequence = "".join(rng.choices("ACGT", k=2400))
            genome = d / "g.fa"
            genome.write_text(">c\n" + sequence + "\n")
            gff = d / "a.gff"
            gff.write_text(
                "c\tx\tgene\t801\t810\t.\t+\t.\tID=plus\nc\tx\tgene\t1601\t1610\t.\t-\t.\tID=minus\n"
            )
            import_reference(
                genome, gff, d / "ref", taxid="32630", assembly="synthetic"
            )
            with ReferenceBundle(d / "ref") as bundle:
                report = run_design(
                    [bundle.resolve_gene("plus"), bundle.resolve_gene("minus")],
                    DesignConfig(num_return=1),
                    d / "out",
                    bundle=bundle,
                    specificity=method,
                    max_product=400,
                )
            self.assertEqual(report["summary"]["successful_records"], 2)
            self.assertEqual(report["summary"]["specificity_errors"], 0)
            for row in report["results"]:
                self.assertEqual(row["specificity_status"], "single_product_" + method)
                self.assertTrue(row["predicted_products"][0]["is_expected"])
                for side in ("left", "right"):
                    a = row[side + "_genomic_start_0based"]
                    b = row[side + "_genomic_end_exclusive"]
                    actual = sequence[a:b]
                    if row[side + "_genomic_strand"] == "-":
                        actual = str(Seq(actual).reverse_complement())
                    self.assertEqual(actual, row[side + "_sequence_5to3"])

    def test_exact_expected_product_both_strands(self):
        self.exercise("exact")

    @unittest.skipUnless(shutil.which("blastn"), "BLAST+ not installed")
    def test_blast_expected_product_both_strands(self):
        self.exercise("blast")
