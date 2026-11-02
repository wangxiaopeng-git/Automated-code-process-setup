import io
import json
import random
import tempfile
import unittest
from pathlib import Path
from unittest.mock import patch
from Bio.Seq import Seq
from Bio.SeqRecord import SeqRecord
from Bio.SeqFeature import SeqFeature, FeatureLocation
from Bio import SeqIO
from batch_primer_design.core import InputError

ASSEMBLY = "GCF_000005845.2"
REPORT = {
    "reports": [
        {
            "accession": ASSEMBLY,
            "organism": {"tax_id": 511145, "organism_name": "E. coli"},
        }
    ]
}
SEQUENCES = {
    "reports": [
        {
            "assembly_accession": ASSEMBLY,
            "refseq_accession": "NC_000913.3",
            "genbank_accession": "U00096.3",
            "chr_name": "1",
            "length": 10000,
            "sequence_name": "Chromosome",
        }
    ]
}


def genbank(start, end):
    rng = random.Random(33)
    seq = "".join(rng.choices("ACGT", k=end - start))
    r = SeqRecord(Seq(seq), id=f"NC_000913.3", annotations={"molecule_type": "DNA"})
    r.features = [
        SeqFeature(
            FeatureLocation(2, 8, strand=-1),
            type="gene",
            qualifiers={"gene": ["example"]},
        )
    ]
    out = io.StringIO()
    SeqIO.write(r, out, "genbank")
    return out.getvalue()


class OnlineReferenceTests(unittest.TestCase):
    def client(self, root):
        from batch_primer_design.online_reference import OnlineReference

        return OnlineReference(ASSEMBLY, root)

    def response(self, url, params):
        if url.endswith("/dataset_report"):
            return json.dumps(REPORT)
        if url.endswith("/sequence_reports"):
            return json.dumps(SEQUENCES)
        if url.endswith("efetch.fcgi"):
            return genbank(int(params["seq_start"]) - 1, int(params["seq_stop"]))
        if url.endswith("esummary.fcgi"):
            return json.dumps(
                {
                    "result": {
                        "uids": ["1"],
                        "1": {
                            "genomicinfo": [
                                {
                                    "chraccver": "NC_000913.3",
                                    "chrstart": 400,
                                    "chrstop": 409,
                                }
                            ]
                        },
                    }
                }
            )
        raise AssertionError(url)

    def test_region_is_bounded_and_annotated_without_full_genome(self):
        from batch_primer_design.online_reference import NcbiTransport

        with tempfile.TemporaryDirectory() as d, patch.object(
            NcbiTransport, "fetch_text", side_effect=self.response
        ):
            client = self.client(Path(d))
            target = client.region("a", "1", 400, 410, "-", flank=20)
            self.assertEqual(len(target.record.sequence), 50)
            self.assertEqual((target.start, target.end, target.strand), (380, 430, "-"))
            self.assertEqual(target.record.target_start, 21)
            self.assertEqual(target.genomic_interval(1, 5), (425, 430))
            feature = client.regions[0]["features"][0]
            self.assertEqual(feature["parts"][0]["start_0based"], 382)
            self.assertEqual(feature["parts"][0]["end_exclusive"], 388)

    def test_ambiguous_alias_not_guessed(self):
        from batch_primer_design.online_reference import NcbiTransport

        data = {
            "reports": SEQUENCES["reports"]
            + [
                dict(
                    SEQUENCES["reports"][0],
                    refseq_accession="NC_999999.1",
                    genbank_accession="X99999.1",
                )
            ]
        }

        def fetch(url, params):
            return (
                json.dumps(data)
                if url.endswith("/sequence_reports")
                else self.response(url, params)
            )

        with tempfile.TemporaryDirectory() as d, patch.object(
            NcbiTransport, "fetch_text", side_effect=fetch
        ):
            c = self.client(Path(d))
            with self.assertRaisesRegex(InputError, "ambiguous"):
                c.resolve_contig("1")

    def test_gene_coordinates_must_belong_to_requested_assembly(self):
        from batch_primer_design.online_reference import NcbiTransport

        def fetch(url, params):
            if url.endswith("esummary.fcgi"):
                return json.dumps(
                    {
                        "result": {
                            "1": {
                                "genomicinfo": [
                                    {
                                        "chraccver": "NC_000913.2",
                                        "chrstart": 400,
                                        "chrstop": 409,
                                    }
                                ]
                            }
                        }
                    }
                )
            return self.response(url, params)

        with tempfile.TemporaryDirectory() as d, patch.object(
            NcbiTransport, "fetch_text", side_effect=fetch
        ):
            with self.assertRaisesRegex(InputError, "assembly"):
                self.client(Path(d)).gene("1")

    def test_metadata_rejects_other_assembly(self):
        from batch_primer_design.online_reference import NcbiTransport

        def fetch(url, params):
            if url.endswith("/dataset_report"):
                return json.dumps(
                    {
                        "reports": [
                            dict(REPORT["reports"][0], accession="GCF_999999999.1")
                        ]
                    }
                )
            return self.response(url, params)

        with tempfile.TemporaryDirectory() as d, patch.object(
            NcbiTransport, "fetch_text", side_effect=fetch
        ):
            with self.assertRaises(InputError):
                self.client(Path(d))

    def test_whole_chromosome_request_exceeds_lightweight_cap(self):
        from batch_primer_design.online_reference import NcbiTransport

        with tempfile.TemporaryDirectory() as d, patch.object(
            NcbiTransport, "fetch_text", side_effect=self.response
        ):
            c = self.client(Path(d))
            c.max_region_bp = 100
            with self.assertRaises(InputError):
                c.region("a", "1", 100, 1000)

    def test_gene_id_resolves_and_rejects_symbols(self):
        from batch_primer_design.online_reference import NcbiTransport

        with tempfile.TemporaryDirectory() as d, patch.object(
            NcbiTransport, "fetch_text", side_effect=self.response
        ):
            c = self.client(Path(d))
            t = c.gene("1", flank=20)
            self.assertEqual((t.start, t.end), (380, 430))
            with self.assertRaisesRegex(InputError, "numeric"):
                c.gene("BRCA1")

    def test_sequence_report_pagination_resolves_later_contigs(self):
        from batch_primer_design.online_reference import NcbiTransport

        def fetch(url, params):
            if url.endswith("/sequence_reports"):
                if params.get("page_token") == "second":
                    return json.dumps(
                        {
                            "reports": [
                                dict(
                                    SEQUENCES["reports"][0],
                                    refseq_accession="NC_999999.1",
                                    genbank_accession="X99999.1",
                                    chr_name="2",
                                )
                            ]
                        }
                    )
                return json.dumps(
                    dict(SEQUENCES, total_count=2, next_page_token="second")
                )
            return self.response(url, params)

        with tempfile.TemporaryDirectory() as d, patch.object(
            NcbiTransport, "fetch_text", side_effect=fetch
        ):
            c = self.client(Path(d))
            self.assertEqual(c.resolve_contig("2"), "NC_999999.1")

    def test_truncated_sequence_report_is_rejected(self):
        from batch_primer_design.online_reference import NcbiTransport

        def fetch(url, params):
            return (
                json.dumps(dict(SEQUENCES, total_count=2))
                if url.endswith("/sequence_reports")
                else self.response(url, params)
            )

        with tempfile.TemporaryDirectory() as d, patch.object(
            NcbiTransport, "fetch_text", side_effect=fetch
        ):
            with self.assertRaisesRegex(InputError, "incomplete"):
                self.client(Path(d))

    def test_returned_region_version_is_checked(self):
        from batch_primer_design.online_reference import NcbiTransport

        def fetch(url, params):
            return (
                self.response(url, params).replace("NC_000913.3", "NC_000913.2")
                if url.endswith("efetch.fcgi")
                else self.response(url, params)
            )

        with tempfile.TemporaryDirectory() as d, patch.object(
            NcbiTransport, "fetch_text", side_effect=fetch
        ):
            c = self.client(Path(d))
            with self.assertRaisesRegex(InputError, "accession/version"):
                c.region("a", "1", 400, 410, flank=0)

    def test_retry_attempts_obey_throttling(self):
        from batch_primer_design.online_reference import NcbiTransport
        import requests

        class Clock:
            now = 1000.0

            def monotonic(self):
                return self.now

            def sleep(self, delay):
                self.now += delay

        clock = Clock()
        attempts = []

        def get(*args, **kwargs):
            attempts.append(clock.now)
            response = requests.Response()
            response.status_code = 503 if len(attempts) == 1 else 200
            response._content = b"hello"
            response.raw = io.BytesIO(b"hello")
            return response

        with tempfile.TemporaryDirectory() as d, patch(
            "batch_primer_design.online_reference.time.monotonic",
            side_effect=clock.monotonic,
        ), patch(
            "batch_primer_design.online_reference.time.sleep", side_effect=clock.sleep
        ):
            c = NcbiTransport(d)
            with patch.object(c.client, "get", side_effect=get):
                self.assertEqual(
                    c.fetch_text("https://eutils.ncbi.nlm.nih.gov/x", {}), "hello"
                )
            self.assertEqual(len(attempts), 2)
            self.assertGreaterEqual(attempts[1] - attempts[0], 0.4)
            c.close()

    def test_inspect_metadata_contains_download_free_contig_names(self):
        from batch_primer_design.online_reference import NcbiTransport

        with tempfile.TemporaryDirectory() as d, patch.object(
            NcbiTransport, "fetch_text", side_effect=self.response
        ):
            info = self.client(Path(d)).metadata()
            self.assertEqual(info["sequences"]["NC_000913.3"]["chr_name"], "1")
            self.assertEqual(info["sequences"]["NC_000913.3"]["length"], 10000)
            self.assertFalse(info["full_genome_downloaded"])

    def test_response_cache_reuses_and_detects_tampering(self):
        from batch_primer_design.online_reference import NcbiTransport
        import requests

        response = requests.Response()
        response.status_code = 200
        response._content = b"hello"
        response.raw = io.BytesIO(b"hello")
        with tempfile.TemporaryDirectory() as d:
            c = NcbiTransport(d)
            with patch.object(c.client, "get", return_value=response):
                self.assertEqual(
                    c.fetch_text("https://eutils.ncbi.nlm.nih.gov/cache", {}), "hello"
                )
            with patch.object(
                c.client,
                "get",
                side_effect=AssertionError("cache should not contact service"),
            ):
                self.assertEqual(
                    c.fetch_text("https://eutils.ncbi.nlm.nih.gov/cache", {}), "hello"
                )
                path = next(Path(d).glob("*.json"))
                payload = json.loads(path.read_text())
                payload["text"] = "changed"
                path.write_text(json.dumps(payload))
                with self.assertRaisesRegex(InputError, "checksum"):
                    c.fetch_text("https://eutils.ncbi.nlm.nih.gov/cache", {})
            c.close()
