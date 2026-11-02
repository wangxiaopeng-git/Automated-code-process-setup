import contextlib
import io
import json
import random
import tempfile
import unittest
from pathlib import Path
from unittest.mock import patch
from Bio import SeqIO
from Bio.Seq import Seq
from Bio.SeqRecord import SeqRecord
from batch_primer_design.advanced_cli import main
from batch_primer_design.online_reference import NcbiTransport


class OnlineWorkflowTests(unittest.TestCase):
    def setUp(self):
        self.tmp = tempfile.TemporaryDirectory()
        self.addCleanup(self.tmp.cleanup)
        self.root = Path(self.tmp.name)
        rng = random.Random(93)
        self.genome = "".join(rng.choices("ACGT", k=10000))
        self.assembly = "GCF_000005845.2"
        self.bed = self.root / "targets.bed"
        self.bed.write_text("NC_000913.3\t4000\t4010\tgene\t0\t-\n")
        self.ids = self.root / "genes.txt"
        self.ids.write_text("1\n")
        self.queries = {}

    def ncbi(self, url, params):
        if url.endswith("/dataset_report"):
            return json.dumps(
                {
                    "reports": [
                        {"accession": self.assembly, "organism": {"tax_id": 511145}}
                    ]
                }
            )
        if url.endswith("/sequence_reports"):
            return json.dumps(
                {
                    "reports": [
                        {
                            "assembly_accession": self.assembly,
                            "refseq_accession": "NC_000913.3",
                            "length": len(self.genome),
                        }
                    ]
                }
            )
        if url.endswith("esummary.fcgi"):
            return json.dumps(
                {
                    "result": {
                        "1": {
                            "genomicinfo": [
                                {
                                    "chraccver": "NC_000913.3",
                                    "chrstart": 4000,
                                    "chrstop": 4009,
                                }
                            ]
                        }
                    }
                }
            )
        if url.endswith("efetch.fcgi"):
            a, b = int(params["seq_start"]) - 1, int(params["seq_stop"])
            record = SeqRecord(
                Seq(self.genome[a:b]),
                id="NC_000913.3",
                annotations={"molecule_type": "DNA"},
            )
            text = io.StringIO()
            SeqIO.write(record, text, "genbank")
            return text.getvalue()
        raise AssertionError(url)

    def blast(self, method, params, **kwargs):
        if method == "POST":
            self.queries = {
                r.id: str(r.seq)
                for r in SeqIO.parse(io.StringIO(params["QUERY"]), "fasta")
            }
            return "RID = JOB1\n"
        if params.get("FORMAT_OBJECT"):
            return "Status=READY\nThereAreHits=yes\n"
        text = "<BlastOutput><BlastOutput_iterations>"
        for key, seq in self.queries.items():
            text += f"<Iteration><Iteration_query-def>{key}</Iteration_query-def><Iteration_query-len>{len(seq)}</Iteration_query-len><Iteration_hits>"
            found = []
            for strand, needle in [
                ("+", seq),
                ("-", str(Seq(seq).reverse_complement())),
            ]:
                pos = self.genome.find(needle)
                if pos >= 0:
                    found.append((pos, strand))
            if found:
                text += "<Hit><Hit_id>ref|NC_000913.3|</Hit_id><Hit_hsps>"
                for pos, strand in found:
                    a, b = (
                        (pos + 1, pos + len(seq))
                        if strand == "+"
                        else (pos + len(seq), pos + 1)
                    )
                    text += f"<Hsp><Hsp_query-from>1</Hsp_query-from><Hsp_query-to>{len(seq)}</Hsp_query-to><Hsp_hit-from>{a}</Hsp_hit-from><Hsp_hit-to>{b}</Hsp_hit-to><Hsp_gaps>0</Hsp_gaps></Hsp>"
                text += "</Hit_hsps></Hit>"
            text += "</Iteration_hits></Iteration>"
        return text + "</BlastOutput_iterations></BlastOutput>"

    def run_cli(self, args):
        with contextlib.redirect_stdout(io.StringIO()), contextlib.redirect_stderr(
            io.StringIO()
        ):
            return main(args)

    def test_region_fetch_and_batched_remote_screening_outputs(self):
        from batch_primer_design.remote_blast import RemoteBlast

        out = self.root / "out"
        with patch.object(
            NcbiTransport, "fetch_text", side_effect=self.ncbi
        ), patch.object(RemoteBlast, "request", side_effect=self.blast):
            code = self.run_cli(
                [
                    "design",
                    str(self.bed),
                    "--assembly",
                    self.assembly,
                    "--num-return",
                    "2",
                    "-o",
                    str(out),
                ]
            )
        self.assertEqual(code, 0)
        report = json.loads((out / "results.json").read_text())
        self.assertEqual(report["summary"]["candidate_pairs"], 2)
        self.assertEqual(
            len(self.queries), 3
        )  # repeated primer shared by two pairs is submitted once
        self.assertTrue(
            all(
                r["specificity_status"] == "single_product_remote_limited"
                for r in report["results"]
            )
        )
        self.assertTrue((out / "online_reference.json").exists())
        self.assertTrue((out / "online_templates.fasta").exists())
        self.assertTrue((out / "region_features.tsv").exists())
        template = next(SeqIO.parse(out / "online_templates.fasta", "fasta"))
        self.assertEqual(
            str(template.seq), str(Seq(self.genome[3700:4310]).reverse_complement())
        )
        self.assertFalse(report["reference"]["full_genome_downloaded"])
        self.assertLess(len(report["reference"]["regions"][0]["features"]), 10)

    def test_fetch_only_gene_id_does_not_submit_blast(self):
        out = self.root / "out"
        with patch.object(NcbiTransport, "fetch_text", side_effect=self.ncbi):
            code = self.run_cli(
                [
                    "design",
                    str(self.ids),
                    "--assembly",
                    self.assembly,
                    "--specificity",
                    "none",
                    "-o",
                    str(out),
                ]
            )
        self.assertEqual(code, 0)
        self.assertEqual(
            json.loads((out / "results.json").read_text())["results"][0][
                "specificity_status"
            ],
            "not_run",
        )

    def test_local_reference_and_assembly_conflict(self):
        with self.assertRaises(SystemExit):
            self.run_cli(
                [
                    "design",
                    str(self.bed),
                    "--assembly",
                    self.assembly,
                    "--reference",
                    "g.fa",
                ]
            )

    def test_no_remote_product_is_unconfirmed_not_unchecked(self):
        from batch_primer_design.remote_blast import OnlineChecker

        def check(rows, targets, **kwargs):
            return {
                (r["record_id"], r["pair_rank"]): {
                    "specificity_status": "no_product_found_remote_limited",
                    "specificity_method": "remote",
                }
                for r in rows
                if r["status"] == "ok"
            }

        out = self.root / "out"
        with patch.object(
            NcbiTransport, "fetch_text", side_effect=self.ncbi
        ), patch.object(OnlineChecker, "check", side_effect=check):
            self.assertEqual(
                self.run_cli(
                    [
                        "design",
                        str(self.bed),
                        "--assembly",
                        self.assembly,
                        "--num-return",
                        "1",
                        "-o",
                        str(out),
                    ]
                ),
                0,
            )
        row = json.loads((out / "results.json").read_text())["results"][0]
        self.assertEqual(row["screening_priority"], "unconfirmed_target")
