import tempfile
import unittest
from pathlib import Path
from batch_primer_design.core import DesignConfig, InputError, read_records, design_record, TemplateRecord
from batch_primer_design.cli import main

class RegressionTests(unittest.TestCase):
    def test_empty_fasta_title_is_input_error(self):
        with tempfile.TemporaryDirectory() as d:
            p=Path(d)/'bad.fa'; p.write_text('>\nACGT\n')
            with self.assertRaises(InputError): list(read_records(p))

    def test_empty_table_not_success(self):
        with tempfile.TemporaryDirectory() as d:
            p=Path(d)/'empty.csv';p.write_text('id,sequence\n')
            self.assertEqual(main([str(p),'-o',str(Path(d)/'out.csv')]),2)

    def test_invalid_chemistry_rejected(self):
        for kw in ({'salt_monovalent_mm':-1},{'primer_dna_nm':float('nan')},{'dntp_mm':float('inf')},{'primer_dna_nm':0}):
            with self.subTest(kw=kw),self.assertRaises(InputError): DesignConfig(**kw).validate()

    def test_engine_error_is_not_no_pair(self):
        row=design_record(TemplateRecord('a','ACGT'*50),DesignConfig(),lambda **kw:{'PRIMER_ERROR':'bad chemistry'})[0]
        self.assertEqual(row['status'],'error')
        self.assertIn('bad chemistry',row['message'])
