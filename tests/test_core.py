import sys
import tempfile
import unittest
from pathlib import Path

sys.path.insert(0, str(Path(__file__).resolve().parents[1] / "src"))

from batch_primer_design.core import (
    DesignConfig,
    InputError,
    TemplateRecord,
    design_batch,
    design_record,
    read_records,
)


def fake_engine(seq_args, global_args):
    assert seq_args["SEQUENCE_TARGET"] == [9, 12]
    assert global_args["PRIMER_PRODUCT_SIZE_RANGE"] == [[100, 300]]
    return {
        "PRIMER_PAIR_NUM_RETURNED": 1,
        "PRIMER_LEFT_0": [4, 20],
        "PRIMER_RIGHT_0": [153, 21],
        "PRIMER_LEFT_0_SEQUENCE": "ACGTACGTACGTACGTACGT",
        "PRIMER_RIGHT_0_SEQUENCE": "TGCATGCATGCATGCATGCAT",
        "PRIMER_LEFT_0_TM": 59.9876,
        "PRIMER_RIGHT_0_TM": 60.1234,
        "PRIMER_LEFT_0_GC_PERCENT": 50.0,
        "PRIMER_RIGHT_0_GC_PERCENT": 47.619,
        "PRIMER_PAIR_0_PRODUCT_SIZE": 150,
        "PRIMER_PAIR_0_PENALTY": 0.12345678,
    }


class CoreTests(unittest.TestCase):
    def test_design_record_converts_coordinates(self):
        record = TemplateRecord("sample", "ACGT" * 50, 10, 12)
        row = design_record(record, DesignConfig(), fake_engine)[0]
        self.assertEqual(row["status"], "ok")
        self.assertEqual(row["left_binding_start_1based"], 5)
        self.assertEqual(row["left_binding_end_1based"], 24)
        self.assertEqual(row["right_binding_start_1based"], 134)
        self.assertEqual(row["right_binding_end_1based"], 154)

    def test_no_pair_keeps_explanation(self):
        def empty_engine(**kwargs):
            return {
                "PRIMER_PAIR_NUM_RETURNED": 0,
                "PRIMER_PAIR_EXPLAIN": "considered 10, unacceptable product size 10, ok 0",
            }

        row = design_record(
            TemplateRecord("x", "ACGT" * 50), DesignConfig(), empty_engine
        )[0]
        self.assertEqual(row["status"], "no_pair")
        self.assertIn("unacceptable product size", row["message"])

    def test_csv_input_and_target_validation(self):
        with tempfile.TemporaryDirectory() as tmp:
            path = Path(tmp) / "input.csv"
            path.write_text(
                "id,sequence,target_start,target_length\nA,acgtacgtacgt,3,4\n",
                encoding="utf-8",
            )
            record = list(read_records(path))[0]
        self.assertEqual(record, TemplateRecord("A", "ACGTACGTACGT", 3, 4))

    def test_invalid_base_is_rejected(self):
        with tempfile.TemporaryDirectory() as tmp:
            path = Path(tmp) / "input.fasta"
            path.write_text(">A\nACGTZ\n", encoding="utf-8")
            rows = list(design_batch(read_records(path), DesignConfig(), fake_engine))
        self.assertEqual(rows[0]["status"], "error")
        self.assertIn("unsupported base", rows[0]["message"])

    def test_duplicate_id_is_reported(self):
        def never_called(**kwargs):
            return {"PRIMER_PAIR_NUM_RETURNED": 0}

        records = [
            TemplateRecord("same", "ACGT" * 50),
            TemplateRecord("same", "ACGT" * 50),
        ]
        rows = list(design_batch(records, DesignConfig(), never_called))
        self.assertEqual(rows[1]["status"], "error")
        self.assertEqual(rows[1]["message"], "duplicate sample id")


if __name__ == "__main__":
    unittest.main()

