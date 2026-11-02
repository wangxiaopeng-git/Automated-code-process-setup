import sys
import unittest
from pathlib import Path

sys.path.insert(0, str(Path(__file__).resolve().parents[1] / "src"))

from batch_primer_design.core import InputError
from batch_primer_design.web_gui import _input_from_payload


class WebGuiInputTests(unittest.TestCase):
    def test_plain_pasted_sequence_becomes_fasta(self):
        filename, content = _input_from_payload(
            {
                "input_mode": "paste",
                "sample_id": "Gene A / test",
                "content": "acgt acgt\nACGT",
            }
        )
        self.assertEqual(filename, "Gene_A_test.fasta")
        self.assertEqual(content, ">Gene_A_test\nACGTACGTACGT\n")

    def test_pasted_fasta_is_preserved(self):
        filename, content = _input_from_payload(
            {"input_mode": "paste", "content": ">one\nACGT\n>two\nTGCA"}
        )
        self.assertEqual(filename, "pasted_sequences.fasta")
        self.assertEqual(content, ">one\nACGT\n>two\nTGCA\n")

    def test_file_mode_requires_filename(self):
        with self.assertRaisesRegex(InputError, "请选择有效"):
            _input_from_payload({"input_mode": "file", "content": "ACGT"})


if __name__ == "__main__":
    unittest.main()
