import hashlib
import io
import json
import tempfile
import unittest
import zipfile
from pathlib import Path
from unittest.mock import patch
import requests
from batch_primer_design.core import InputError


class NcbiTests(unittest.TestCase):
    def archive(self, bad=False):
        content = io.BytesIO()
        with zipfile.ZipFile(content, "w") as z:
            if bad:
                z.writestr("../escape", "bad")
            z.writestr(
                "ncbi_dataset/data/assembly_data_report.jsonl",
                json.dumps(
                    {"accession": "GCF_000005845.2", "organism": {"taxId": 511145}}
                )
                + "\n",
            )
            z.writestr(
                "ncbi_dataset/data/GCF_000005845.2/genomic.fna", ">c\nACGTACGT\n"
            )
            z.writestr(
                "ncbi_dataset/data/GCF_000005845.2/genomic.gff",
                "c\tx\tgene\t1\t4\t.\t+\t.\tID=a\n",
            )
        return content.getvalue()

    def response(self, content):
        r = requests.Response()
        r.status_code = 200
        r._content = content
        r.raw = io.BytesIO(content)
        return r

    def test_search_retains_accession_and_handles_query(self):
        from batch_primer_design.ncbi import search_ncbi

        data = {
            "reports": [
                {
                    "accession": "GCF_000005845.2",
                    "organism": {"tax_id": 511145},
                    "assembly_info": {"assembly_name": "ASM584v2"},
                }
            ]
        }
        with patch(
            "requests.Session.get",
            return_value=self.response(json.dumps(data).encode()),
        ):
            r = search_ncbi("Escherichia coli")
            self.assertEqual(r["reports"][0]["accession"], "GCF_000005845.2")

    def test_download_imports_reference(self):
        from batch_primer_design.ncbi import download_ncbi

        with (
            tempfile.TemporaryDirectory() as d,
            patch("requests.Session.get", return_value=self.response(self.archive())),
        ):
            manifest = download_ncbi("GCF_000005845.2", Path(d) / "ref", taxid="511145")
            self.assertEqual(manifest["assembly"], "GCF_000005845.2")
            self.assertIn("genome.fa", manifest["files"])
            self.assertTrue((Path(d) / "ref" / "download.json").exists())

    def test_unsafe_archive_never_imported(self):
        from batch_primer_design.ncbi import download_ncbi

        with (
            tempfile.TemporaryDirectory() as d,
            patch(
                "requests.Session.get", return_value=self.response(self.archive(True))
            ),
        ):
            with self.assertRaises(InputError):
                download_ncbi("GCF_000005845.2", Path(d) / "ref", taxid="511145")
            self.assertFalse((Path(d) / "ref").exists())
            self.assertFalse((Path(d) / "escape").exists())

    def test_accession_must_have_version(self):
        from batch_primer_design.ncbi import download_ncbi

        with tempfile.TemporaryDirectory() as d, self.assertRaises(InputError):
            download_ncbi("GCF_000005845", Path(d) / "ref", taxid="511145")

    def test_taxonomy_mismatch_rejected(self):
        from batch_primer_design.ncbi import download_ncbi

        with (
            tempfile.TemporaryDirectory() as d,
            patch("requests.Session.get", return_value=self.response(self.archive())),
        ):
            with self.assertRaisesRegex(InputError, "tax"):
                download_ncbi("GCF_000005845.2", Path(d) / "ref", taxid="9606")
            self.assertFalse((Path(d) / "ref").exists())
