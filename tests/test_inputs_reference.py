import gzip
import tempfile
import unittest
from pathlib import Path
from batch_primer_design.core import InputError


class InputsReferenceTests(unittest.TestCase):
    def setUp(self):
        self.tmp = tempfile.TemporaryDirectory()
        self.addCleanup(self.tmp.cleanup)
        self.root = Path(self.tmp.name)

    def test_gzip_fasta(self):
        from batch_primer_design.inputs import load_templates

        p = self.root / "a.fa.gz"
        with gzip.open(p, "wt") as h:
            h.write(">a\nacgtacgt\n")
        self.assertEqual(load_templates(p)[0].sequence, "ACGTACGT")

    def test_xlsx(self):
        from openpyxl import Workbook
        from batch_primer_design.inputs import load_templates

        w = Workbook()
        w.active.append(["id", "sequence"])
        w.active.append(["a", "acgt"])
        p = self.root / "a.xlsx"
        w.save(p)
        self.assertEqual(load_templates(p)[0].sequence, "ACGT")

    def test_genbank(self):
        from Bio.Seq import Seq
        from Bio.SeqRecord import SeqRecord
        from Bio import SeqIO
        from batch_primer_design.inputs import load_templates

        p = self.root / "a.gb"
        SeqIO.write(
            SeqRecord(Seq("ACGT"), id="gene1", annotations={"molecule_type": "DNA"}),
            p,
            "genbank",
        )
        self.assertEqual(load_templates(p)[0].sample_id, "gene1")

    def bundle(self, annotation):
        from batch_primer_design.reference import import_reference

        fa = self.root / "g.fa"
        fa.write_text(">chr1\nAAAACCCCGGGGTTTT\n")
        gff = self.root / "a.gff"
        gff.write_text(annotation)
        return import_reference(
            fa, gff, self.root / "bundle", taxid="123", assembly="local-v1"
        )

    def test_negative_gene_and_bed_mapping(self):
        from batch_primer_design.reference import ReferenceBundle

        self.bundle("chr1\tx\tgene\t5\t8\t.\t-\t.\tID=g1;Name=GeneA\n")
        with ReferenceBundle(self.root / "bundle") as b:
            target = b.resolve_gene("GeneA", flank=2)
            self.assertEqual(target.record.sequence, "CCGGGGTT")
            self.assertEqual(
                (target.record.target_start, target.record.target_length), (3, 4)
            )
            self.assertEqual(target.genomic_interval(1, 4), (6, 10))
            p = self.root / "a.bed"
            p.write_text("chr1\t4\t8\ta\t0\t+\n")
            self.assertEqual(b.read_bed(p, flank=0)[0].record.sequence, "CCCC")

    def test_mismatched_contig_rejected_without_partial_bundle(self):
        with self.assertRaises(InputError):
            self.bundle("chr2\tx\tgene\t5\t8\t.\t+\t.\tID=g1\n")
        self.assertFalse((self.root / "bundle").exists())

    def test_gtf_aggregates_gene(self):
        from batch_primer_design.reference import ReferenceBundle

        self.bundle(
            'chr1\tx\texon\t3\t5\t.\t+\t.\tgene_id "g1"; transcript_id "t1"; gene_name "A";\nchr1\tx\texon\t9\t12\t.\t+\t.\tgene_id "g1"; transcript_id "t1";\n'
        )
        with ReferenceBundle(self.root / "bundle") as b:
            t = b.resolve_gene("g1", flank=0)
            self.assertEqual(t.record.sequence, "AACCCCGGGG")

    def test_ambiguous_alias_rejected(self):
        from batch_primer_design.reference import ReferenceBundle

        self.bundle(
            "chr1\tx\tgene\t1\t4\t.\t+\t.\tID=g1;Name=A\nchr1\tx\tgene\t9\t12\t.\t+\t.\tID=g2;Name=A\n"
        )
        with ReferenceBundle(self.root / "bundle") as b:
            with self.assertRaisesRegex(InputError, "ambiguous"):
                b.resolve_gene("A")
