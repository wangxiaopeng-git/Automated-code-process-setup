"""NCBI BLAST URLAPI batch screening. A limited search, never an assembly-wide proof."""

from __future__ import annotations
import hashlib
import json
import re
import time
import xml.etree.ElementTree as ET
from dataclasses import replace
from pathlib import Path
from Bio.Seq import Seq
import requests
from .core import InputError
from .specificity import Hit, pair_hits

URL = "https://blast.ncbi.nlm.nih.gov/Blast.cgi"


class RemoteBlast:
    def __init__(
        self, cache_dir, *, email=None, hit_limit=100, poll_timeout=900, max_queries=40
    ):
        if (
            not 1 <= hit_limit <= 500
            or not 1 <= max_queries <= 100
            or not 60 <= poll_timeout <= 3600
        ):
            raise InputError(
                "remote hit limit 1..500, max queries 1..100, timeout 60..3600 seconds"
            )
        self.cache_dir = Path(cache_dir)
        self.email = email
        self.hit_limit = hit_limit
        self.poll_timeout = poll_timeout
        self.max_queries = max_queries
        self.client = requests.Session()
        self.client.headers["User-Agent"] = "batch-primer-designer/0.3.0"
        self.last_request = 0.0
        self.last_poll = 0.0
        self.deadline = None
        self.poll_checkpoint = None

    def close(self):
        self.client.close()

    def request(self, method, params, *, poll=False):
        # No implicit retry: a failed POST may already have created a remote job.
        wait = max(
            0.0,
            10 - (time.monotonic() - self.last_request),
            60 - (time.monotonic() - self.last_poll) if poll else 0.0,
        )
        if self.deadline is not None and time.monotonic() + wait >= self.deadline:
            raise InputError(
                f'online BLAST polling timeout for RID {params.get("RID", "pending")}; rerun with the same cache to resume'
            )
        if wait > 0:
            time.sleep(wait)
        self.last_request = time.monotonic()
        if poll:
            self.last_poll = self.last_request
            if self.poll_checkpoint:
                self.poll_checkpoint(time.time())
        params = dict(params, TOOL="batch_primer_designer")
        if self.email:
            params["EMAIL"] = self.email
        try:
            response = self.client.request(
                method,
                URL,
                **({"data": params} if method == "POST" else {"params": params}),
                stream=True,
                timeout=(
                    (
                        min(15, max(0.1, self.deadline - time.monotonic())),
                        min(120, max(0.1, self.deadline - time.monotonic())),
                    )
                    if self.deadline is not None
                    else (15, 120)
                ),
            )
            try:
                response.raise_for_status()
                chunks = []
                size = 0
                for chunk in response.iter_content(65536):
                    if self.deadline is not None and time.monotonic() >= self.deadline:
                        raise InputError(
                            "online BLAST polling timeout while reading response; resume with the same cache"
                        )
                    size += len(chunk)
                    if size > 50_000_000:
                        raise InputError("remote BLAST response exceeds 50 MB limit")
                    chunks.append(chunk)
                return b"".join(chunks).decode("utf-8")
            finally:
                response.close()
        except (requests.RequestException, UnicodeError) as exc:
            raise InputError(
                f"online BLAST request failed; do not assume screening completed: {exc}"
            ) from exc

    def run(self, queries, taxid):
        self.deadline = None
        self.poll_checkpoint = None
        if not queries or len(queries) > self.max_queries:
            raise InputError(
                f"online BLAST supports at most {self.max_queries} primer queries per run; reduce targets/num-return or use local reference"
            )
        if any(
            not 12 <= len(s) <= 60 or re.search("[^ACGTN]", s) for s in queries.values()
        ):
            raise InputError("online BLAST requires DNA primers of 12..60 nt")
        fasta = "".join(f">{key}\n{seq}\n" for key, seq in queries.items())
        params = {
            "CMD": "Put",
            "PROGRAM": "blastn",
            "DATABASE": "nt",
            "QUERY": fasta,
            "ENTREZ_QUERY": f"txid{int(taxid)}[ORGN]",
            "EXPECT": 1000,
            "WORD_SIZE": 7,
            "FILTER": "F",
            "HITLIST_SIZE": self.hit_limit,
            "NUCL_REWARD": 1,
            "NUCL_PENALTY": -3,
        }
        key = hashlib.sha256(json.dumps(params, sort_keys=True).encode()).hexdigest()
        self.cache_dir.mkdir(parents=True, exist_ok=True)
        path = self.cache_dir / (key + ".blast.json")
        xml_path = self.cache_dir / (key + ".xml")
        job = json.loads(path.read_text()) if path.exists() else {}
        if job and time.time() - job.get("submitted_at", 0) < 86400:
            if job.get("failed"):
                raise InputError(
                    "previous online BLAST job failed; inspect cache/job record and retry with a new cache directory"
                )
            if xml_path.exists() and job.get("xml_sha256"):
                xml = xml_path.read_text()
                if hashlib.sha256(xml.encode()).hexdigest() != job["xml_sha256"]:
                    raise InputError("cached online BLAST XML checksum mismatch")
                return xml, job
        else:
            job = {}
        if not job:
            result = self.request("POST", params)
            match = re.search(r"\bRID\s*=\s*([A-Za-z0-9_-]+)", result)
            if not match:
                raise InputError(
                    "NCBI did not return a BLAST RID; request was not silently resubmitted"
                )
            job = {
                "rid": match.group(1),
                "submitted_at": time.time(),
                "database": "nt",
                "taxid": str(taxid),
                "hit_limit": self.hit_limit,
                "query_sha256": hashlib.sha256(fasta.encode()).hexdigest(),
            }
            path.write_text(json.dumps(job, indent=2))
        print(
            f"在线 BLAST RID={job['rid']}；遵守服务限流，结果可能需要数分钟。",
            flush=True,
        )
        if job.get("last_poll_at"):
            self.last_poll = time.monotonic() - max(
                0, time.time() - job["last_poll_at"]
            )
        self.deadline = time.monotonic() + self.poll_timeout

        def checkpoint(dispatched_at):
            job["last_poll_at"] = dispatched_at
            path.write_text(json.dumps(job, indent=2))

        self.poll_checkpoint = checkpoint
        while time.monotonic() < self.deadline:
            status = self.request(
                "GET",
                {"CMD": "Get", "RID": job["rid"], "FORMAT_OBJECT": "SearchInfo"},
                poll=True,
            )
            if "Status=WAITING" in status:
                continue
            if "Status=FAILED" in status or "Status=UNKNOWN" in status:
                job["failed"] = True
                path.write_text(json.dumps(job, indent=2))
                raise InputError(
                    f'online BLAST job {job["rid"]}: FAILED/UNKNOWN; no completed specificity result'
                )
            if "Status=READY" not in status:
                raise InputError("unrecognized NCBI BLAST status response")
            xml = self.request(
                "GET",
                {
                    "CMD": "Get",
                    "RID": job["rid"],
                    "FORMAT_TYPE": "XML",
                    "HITLIST_SIZE": self.hit_limit,
                },
            )
            job["xml_sha256"] = hashlib.sha256(xml.encode()).hexdigest()
            job["completed_at"] = time.time()
            xml_path.write_text(xml)
            path.write_text(json.dumps(job, indent=2))
            return xml, job
        raise InputError(
            f'online BLAST polling timeout for RID {job["rid"]}; rerun with the same inputs/cache to resume pending job'
        )


def decode_hits(xml, queries, reference, *, max_mismatches, hit_limit):
    if not 0 <= max_mismatches <= 5:
        raise InputError("max-mismatches must be 0..5")
    if "<!ENTITY" in xml:
        raise InputError("unexpected XML entity declaration in BLAST response")
    try:
        root = ET.fromstring(xml)
    except ET.ParseError as exc:
        raise InputError("online BLAST returned invalid XML") from exc
    if root.tag != "BlastOutput":
        raise InputError("expected NCBI BLAST XML output")
    results = {key: [] for key in queries}
    seen = set()
    sites = 0
    limits = {
        "hit_limit_reached": False,
        "missing_queries": [],
        "excluded_sequence_hits": 0,
        "unversioned_hits": 0,
    }
    for iteration in root.findall("./BlastOutput_iterations/Iteration"):
        title = (iteration.findtext("Iteration_query-def") or "").split()
        if not title:
            raise InputError("online BLAST XML missing query definition")
        key = title[0]
        if key not in queries:
            raise InputError("online BLAST returned an unexpected query ID")
        seen.add(key)
        primer = queries[key]
        if int(iteration.findtext("Iteration_query-len", "0")) != len(primer):
            raise InputError("online BLAST query length mismatch")
        hits = iteration.findall("./Iteration_hits/Hit")
        if len(hits) >= hit_limit:
            limits["hit_limit_reached"] = True
        for hit in hits:
            identity = " ".join(
                hit.findtext(k, "") for k in ("Hit_accession", "Hit_id", "Hit_def")
            )
            accessions = re.findall(r"[A-Z]{1,6}_?\d+\.\d+", identity)
            known = {
                reference.accession_aliases[a]
                for a in accessions
                if a in reference.accession_aliases
            }
            if not known:
                limits["excluded_sequence_hits"] += 1
                if not accessions:
                    limits["unversioned_hits"] += 1
                continue
            for contig in sorted(known):
                for hsp in hit.findall("./Hit_hsps/Hsp"):
                    if int(hsp.findtext("Hsp_gaps", "0")):
                        continue
                    qstart, qend = int(hsp.findtext("Hsp_query-from")), int(
                        hsp.findtext("Hsp_query-to")
                    )
                    a, b = int(hsp.findtext("Hsp_hit-from")), int(
                        hsp.findtext("Hsp_hit-to")
                    )
                    if not (
                        1 <= qstart <= len(primer)
                        and 1 <= qend <= len(primer)
                        and a >= 1
                        and b >= 1
                    ):
                        raise InputError(
                            "online BLAST contains invalid alignment coordinates"
                        )
                    if abs(a - b) != abs(qstart - qend):
                        raise InputError(
                            "online BLAST ungapped alignment span mismatch"
                        )
                    if qstart > qend:
                        qstart, qend = qend, qstart
                        a, b = b, a
                    strand = "+" if a <= b else "-"
                    if strand == "+":
                        start = a - qstart
                        end = start + len(primer)
                    else:
                        end = a + qstart - 1
                        start = end - len(primer)
                    if start < 0 or end > reference.sequences[contig]["length"]:
                        continue
                    sites += 1
                    if sites > 500:
                        raise InputError(
                            "online hit extension exceeds 500 sites; screening incomplete, use local reference"
                        )
                    _, record = reference.fetch_record(contig, start, end)
                    segment = record.seq
                    if strand == "-":
                        segment = Seq(str(segment)).reverse_complement()
                    mismatch = tuple(
                        i + 1
                        for i, (x, y) in enumerate(zip(primer, str(segment).upper()))
                        if x != y
                    )
                    if len(mismatch) <= max_mismatches:
                        results[key].append(
                            Hit(key, contig, start, end, strand, mismatch)
                        )
    for key in results:
        results[key] = list(set(results[key]))
    limits["missing_queries"] = sorted(set(queries) - seen)
    return results, limits


class OnlineChecker:
    def __init__(self, reference, cache_dir, **kwargs):
        self.reference = reference
        self.blast = RemoteBlast(cache_dir, **kwargs)

    def close(self):
        self.blast.close()

    def check(self, rows, targets, *, max_product, max_mismatches):
        queries = {}
        names = {}
        by_sequence = {}
        for row in rows:
            if row["status"] != "ok":
                continue
            key = (row["record_id"], row["pair_rank"])
            left = f"r{key[0]}p{key[1]}L"
            right = f"r{key[0]}p{key[1]}R"
            for name, sequence in [
                (left, row["left_sequence_5to3"]),
                (right, row["right_sequence_5to3"]),
            ]:
                if sequence not in by_sequence:
                    by_sequence[sequence] = name
                    queries[name] = sequence
            names[key] = (
                by_sequence[row["left_sequence_5to3"]],
                by_sequence[row["right_sequence_5to3"]],
            )
        if not queries:
            return {}
        xml, job = self.blast.run(queries, self.reference.taxid)
        hits, limits = decode_hits(
            xml,
            queries,
            self.reference,
            max_mismatches=max_mismatches,
            hit_limit=self.blast.hit_limit,
        )
        results = {}
        for row in rows:
            if row["status"] != "ok":
                continue
            key = (row["record_id"], row["pair_rank"])
            left, right = names[key]
            pair = [replace(h, primer="L") for h in hits[left]] + [
                replace(h, primer="R") for h in hits[right]
            ]
            products = pair_hits(
                pair, len(queries[left]) + len(queries[right]), max_product
            )
            target = targets[row["record_id"] - 1]
            expected = None
            if target.contig:
                a, b = target.genomic_interval(
                    row["left_binding_start_1based"], row["right_binding_end_1based"]
                )
                expected = (target.contig, a, b, target.strand)
            for product in products:
                product["is_expected"] = (
                    (
                        product["contig"],
                        product["start"],
                        product["end"],
                        product["strand"],
                    )
                    == expected
                    if expected
                    else None
                )
            incomplete = (
                limits["hit_limit_reached"]
                or limits["missing_queries"]
                or limits["unversioned_hits"]
            )
            if incomplete:
                status = "incomplete_remote"
            elif expected and not any(p["is_expected"] for p in products):
                status = "expected_not_found_remote_limited"
            elif len(products) > 1:
                status = "multiple_products_remote"
            elif len(products) == 1:
                status = (
                    "single_product_remote_limited"
                    if expected
                    else "single_product_unassigned_remote_limited"
                )
            else:
                status = "no_product_found_remote_limited"
            results[key] = {
                "specificity_status": status,
                "specificity_method": "remote",
                "specificity_reference": self.reference.assembly,
                "predicted_products": products,
                "specificity_hit_count": len(pair),
                "remote_blast_job": job,
                "specificity_limits": {
                    **limits,
                    "whole_assembly_coverage_verified": False,
                    "database": "NCBI nt (dynamic)",
                    "organism_filter_taxid": str(self.reference.taxid),
                    "assembly_filter": "returned versioned sequence IDs intersected with requested assembly report",
                    "max_product_bp": max_product,
                    "max_mismatches": max_mismatches,
                    "word_size": 7,
                    "scope": "heuristic online organism search; no proof all assembly contigs are indexed; gaps and same-primer products excluded; partial HSPs extended by EFetch",
                },
            }
        return results
