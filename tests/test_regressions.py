import tempfile
import unittest
from pathlib import Path
from batch_primer_design.core import (
    DesignConfig,
    InputError,
    read_records,
    design_record,
    TemplateRecord,
)
from batch_primer_design.cli import main


class RegressionTests(unittest.TestCase):
    def test_empty_fasta_title_is_input_error(self):
        with tempfile.TemporaryDirectory() as d:
            p = Path(d) / "bad.fa"
            p.write_text(">\nACGT\n")
            with self.assertRaises(InputError):
                list(read_records(p))

    def test_empty_table_not_success(self):
        with tempfile.TemporaryDirectory() as d:
            p = Path(d) / "empty.csv"
            p.write_text("id,sequence\n")
            self.assertEqual(main([str(p), "-o", str(Path(d) / "out.csv")]), 2)

    def test_invalid_chemistry_rejected(self):
        for kw in (
            {"salt_monovalent_mm": -1},
            {"primer_dna_nm": float("nan")},
            {"dntp_mm": float("inf")},
            {"primer_dna_nm": 0},
        ):
            with self.subTest(kw=kw), self.assertRaises(InputError):
                DesignConfig(**kw).validate()

    def test_engine_error_is_not_no_pair(self):
        row = design_record(
            TemplateRecord("a", "ACGT" * 50),
            DesignConfig(),
            lambda **kw: {"PRIMER_ERROR": "bad chemistry"},
        )[0]
        self.assertEqual(row["status"], "error")
        self.assertIn("bad chemistry", row["message"])


class CountAndExportTests(unittest.TestCase):
    def test_legacy_batch_has_record_identity(self):
        from batch_primer_design.core import design_batch

        rows = list(
            design_batch(
                [TemplateRecord("a", "ACGT" * 50)] * 3,
                DesignConfig(),
                lambda **kw: {"PRIMER_PAIR_NUM_RETURNED": 0},
            )
        )
        self.assertEqual([r["record_id"] for r in rows], [1, 2, 3])

    def test_spreadsheet_output_does_not_execute_sample_id(self):
        import csv
        from batch_primer_design.core import write_results, _blank_result

        with tempfile.TemporaryDirectory() as d:
            p = Path(d) / "out.csv"
            write_results([_blank_result("=1+1", "error", "bad")], p)
            with p.open(encoding="utf-8-sig") as h:
                self.assertEqual(next(csv.DictReader(h))["sample_id"], "'=1+1")
