import json
import tempfile
import unittest
from pathlib import Path
from unittest.mock import patch
from types import SimpleNamespace
from batch_primer_design.core import InputError, TemplateRecord
from batch_primer_design.reference import Target
from batch_primer_design.specificity import Hit

XML = """<?xml version="1.0"?><BlastOutput><BlastOutput_iterations><Iteration><Iteration_query-def>q1</Iteration_query-def><Iteration_query-len>20</Iteration_query-len><Iteration_hits><Hit><Hit_id>ref|NC_000913.3|</Hit_id><Hit_accession>NC_000913</Hit_accession><Hit_hsps><Hsp><Hsp_query-from>2</Hsp_query-from><Hsp_query-to>20</Hsp_query-to><Hsp_hit-from>12</Hsp_hit-from><Hsp_hit-to>30</Hsp_hit-to><Hsp_gaps>0</Hsp_gaps></Hsp></Hit_hsps></Hit></Iteration_hits></Iteration></BlastOutput_iterations></BlastOutput>"""


class FakeReference:
    assembly = "GCF_000005845.2"
    taxid = 511145
    accession_aliases = {"NC_000913.3": "NC_000913.3"}
    sequences = {"NC_000913.3": {"length": 1000}}

    def fetch_record(self, contig, start, end):
        return contig, SimpleNamespace(seq="T" + "A" * 19)


class RemoteBlastTests(unittest.TestCase):
    def test_partial_hsp_is_extended_and_terminal_mismatch_counted(self):
        from batch_primer_design.remote_blast import decode_hits

        hits, limits = decode_hits(
            XML, {"q1": "A" * 20}, FakeReference(), max_mismatches=2, hit_limit=50
        )
        self.assertEqual(hits["q1"], [Hit("q1", "NC_000913.3", 10, 30, "+", (1,))])
        self.assertEqual(limits["missing_queries"], [])

    def test_hit_limit_is_incomplete_even_if_one_assembly_hit(self):
        from batch_primer_design.remote_blast import decode_hits

        _, limits = decode_hits(
            XML, {"q1": "A" * 20}, FakeReference(), max_mismatches=2, hit_limit=1
        )
        self.assertTrue(limits["hit_limit_reached"])

    def test_missing_query_cannot_be_treated_as_no_hits(self):
        from batch_primer_design.remote_blast import decode_hits

        _, limits = decode_hits(
            XML,
            {"q1": "A" * 20, "q2": "C" * 20},
            FakeReference(),
            max_mismatches=2,
            hit_limit=50,
        )
        self.assertEqual(limits["missing_queries"], ["q2"])

    def test_other_sequence_version_is_not_mapped_to_requested_version(self):
        from batch_primer_design.remote_blast import decode_hits

        xml = XML.replace("NC_000913.3", "NC_000913.2")
        hits, limits = decode_hits(
            xml, {"q1": "A" * 20}, FakeReference(), max_mismatches=2, hit_limit=50
        )
        self.assertEqual(hits["q1"], [])
        self.assertEqual(limits["excluded_sequence_hits"], 1)

    def test_failed_status_stops_without_resubmitting(self):
        from batch_primer_design.remote_blast import RemoteBlast

        with tempfile.TemporaryDirectory() as d:
            client = RemoteBlast(d, poll_timeout=120)
            responses = ["RID = SAMPLE123\nRTOE = 10\n", "Status=FAILED\n"]
            with patch.object(client, "request", side_effect=responses):
                with self.assertRaisesRegex(InputError, "FAILED"):
                    client.run({"q1": "A" * 20}, 511145)
            client.close()

    def test_ready_returns_xml_and_records_rid(self):
        from batch_primer_design.remote_blast import RemoteBlast

        with tempfile.TemporaryDirectory() as d:
            client = RemoteBlast(d, poll_timeout=120)
            with patch.object(
                client,
                "request",
                side_effect=[
                    "RID = SAMPLE123\nRTOE = 10\n",
                    "Status=READY\nThereAreHits=yes\n",
                    XML,
                ],
            ):
                xml, job = client.run({"q1": "A" * 20}, 511145)
            self.assertEqual(xml, XML)
            self.assertEqual(job["rid"], "SAMPLE123")
            client.close()

    def test_assembly_check_remains_limited_even_when_expected_is_unique(self):
        from batch_primer_design.remote_blast import OnlineChecker

        fake = FakeReference()
        with tempfile.TemporaryDirectory() as d:
            checker = OnlineChecker(fake, d)
            row = {
                "record_id": 1,
                "pair_rank": 1,
                "status": "ok",
                "left_sequence_5to3": "A" * 20,
                "right_sequence_5to3": "C" * 20,
                "left_binding_start_1based": 1,
                "right_binding_end_1based": 100,
            }
            target = Target(
                TemplateRecord("a", "A" * 150),
                "NC_000913.3",
                0,
                150,
                "+",
                fake.assembly,
            )
            hits = {
                "r1p1L": [Hit("r1p1L", "NC_000913.3", 0, 20, "+")],
                "r1p1R": [Hit("r1p1R", "NC_000913.3", 80, 100, "-")],
            }
            limits = {
                "hit_limit_reached": False,
                "missing_queries": [],
                "excluded_sequence_hits": 0,
                "unversioned_hits": 0,
            }
            with patch.object(
                checker.blast, "run", return_value=(XML, {"rid": "X"})
            ), patch(
                "batch_primer_design.remote_blast.decode_hits",
                return_value=(hits, limits),
            ):
                result = checker.check(
                    [row], [target], max_product=200, max_mismatches=2
                )[(1, 1)]
            self.assertEqual(
                result["specificity_status"], "single_product_remote_limited"
            )
            self.assertFalse(
                result["specificity_limits"]["whole_assembly_coverage_verified"]
            )
            checker.close()

    def test_duplicate_candidate_sequences_are_submitted_once(self):
        from batch_primer_design.remote_blast import OnlineChecker

        fake = FakeReference()
        with tempfile.TemporaryDirectory() as d:
            checker = OnlineChecker(fake, d)
            row = {
                "record_id": 1,
                "pair_rank": 1,
                "status": "ok",
                "left_sequence_5to3": "A" * 20,
                "right_sequence_5to3": "C" * 20,
                "left_binding_start_1based": 1,
                "right_binding_end_1based": 100,
            }
            target = Target(
                TemplateRecord("a", "A" * 150),
                "NC_000913.3",
                0,
                150,
                "+",
                fake.assembly,
            )
            queries = []

            def run(q, taxid):
                queries.extend(q)
                return XML, {"rid": "X"}

            limits = {
                "hit_limit_reached": False,
                "missing_queries": [],
                "excluded_sequence_hits": 0,
                "unversioned_hits": 0,
            }

            def decode(xml, q, *args, **kwargs):
                return {key: [] for key in q}, limits

            with patch.object(checker.blast, "run", side_effect=run), patch(
                "batch_primer_design.remote_blast.decode_hits", side_effect=decode
            ):
                checker.check(
                    [row, dict(row, pair_rank=2)],
                    [target],
                    max_product=200,
                    max_mismatches=2,
                )
            self.assertEqual(len(queries), 2)
            checker.close()

    def test_negative_subject_hsp_extension(self):
        from batch_primer_design.remote_blast import decode_hits

        xml = XML.replace(
            "<Hsp_hit-from>12</Hsp_hit-from>", "<Hsp_hit-from>29</Hsp_hit-from>"
        ).replace("<Hsp_hit-to>30</Hsp_hit-to>", "<Hsp_hit-to>11</Hsp_hit-to>")
        ref = FakeReference()
        ref.fetch_record = lambda *args: ("NC_000913.3", SimpleNamespace(seq="T" * 20))
        hits, _ = decode_hits(
            xml, {"q1": "A" * 20}, ref, max_mismatches=2, hit_limit=50
        )
        self.assertEqual(hits["q1"], [Hit("q1", "NC_000913.3", 10, 30, "-", ())])

    def test_redundant_nt_hit_is_checked_on_all_assembly_contigs(self):
        from batch_primer_design.remote_blast import decode_hits

        ref = FakeReference()
        ref.accession_aliases = dict(ref.accession_aliases, NC_999999_1="unused")
        ref.accession_aliases["NC_999999.1"] = "NC_999999.1"
        ref.sequences = dict(ref.sequences, **{"NC_999999.1": {"length": 1000}})
        xml = XML.replace(
            "<Hit_accession>",
            "<Hit_def>also ref|NC_999999.1| identical</Hit_def><Hit_accession>",
        )
        hits, _ = decode_hits(
            xml, {"q1": "A" * 20}, ref, max_mismatches=2, hit_limit=50
        )
        self.assertEqual({h.contig for h in hits["q1"]}, {"NC_000913.3", "NC_999999.1"})

    def test_poll_deadline_and_resumed_actual_dispatch_interval(self):
        from batch_primer_design.remote_blast import RemoteBlast

        class Clock:
            now = 1000.0

            def monotonic(self):
                return self.now

            def sleep(self, delay):
                self.now += delay

        clock = Clock()
        polls = []

        def response(text):
            return SimpleNamespace(
                raise_for_status=lambda: None,
                iter_content=lambda n: iter([text.encode()]),
                close=lambda: None,
            )

        def request(method, url, **kwargs):
            params = kwargs.get("params", kwargs.get("data", {}))
            if method == "POST":
                return response("RID = SAMPLE123\n")
            if params.get("FORMAT_OBJECT"):
                polls.append(clock.now)
                return response(
                    "Status=WAITING\n" if len(polls) == 1 else "Status=READY\n"
                )
            return response(XML)

        with tempfile.TemporaryDirectory() as d, patch(
            "batch_primer_design.remote_blast.time.monotonic",
            side_effect=clock.monotonic,
        ), patch(
            "batch_primer_design.remote_blast.time.time", side_effect=clock.monotonic
        ), patch(
            "batch_primer_design.remote_blast.time.sleep", side_effect=clock.sleep
        ):
            c = RemoteBlast(d, poll_timeout=60)
            with patch.object(
                c.client, "request", side_effect=request
            ), self.assertRaisesRegex(InputError, "timeout"):
                c.run({"q1": "A" * 20}, 511145)
            self.assertTrue(all(t < 1060 for t in polls))
            cached = json.loads(next(Path(d).glob("*.blast.json")).read_text())
            self.assertEqual(cached["last_poll_at"], polls[-1])
            c.close()
            clock.now = 1040
            c = RemoteBlast(d, poll_timeout=120)
            with patch.object(c.client, "request", side_effect=request):
                c.run({"q1": "A" * 20}, 511145)
            self.assertGreaterEqual(polls[-1] - polls[-2], 60)
            c.close()
