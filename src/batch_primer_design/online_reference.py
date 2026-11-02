"""Assembly-pinned metadata and bounded NCBI E-utilities regions, not full genomes."""

from __future__ import annotations
import hashlib
import io
import json
import re
import tempfile
import time
from datetime import datetime, timezone
from pathlib import Path
from Bio import SeqIO
from Bio.Seq import Seq
import requests
from .core import InputError, TemplateRecord, validate_record
from .ncbi import BASE, session
from .reference import Target

EUTILS = "https://eutils.ncbi.nlm.nih.gov/entrez/eutils/"
VERSIONED = re.compile(r"^[A-Z][A-Z0-9_]*\d+\.\d+$")


class NcbiTransport:
    """Public API requests at <=2.5/s, bounded responses and 7-day checksum caches."""

    def __init__(self, cache_dir, email=None):
        self.cache_dir = Path(cache_dir)
        self.email = email
        self.client = requests.Session()
        self.last_request = 0.0
        self.client.headers["User-Agent"] = "batch-primer-designer/0.3.0"
        self.transfers = []

    def close(self):
        self.client.close()

    def fetch_text(self, url, params):
        key = hashlib.sha256(
            json.dumps([url, params], sort_keys=True).encode()
        ).hexdigest()
        path = self.cache_dir / (key + ".json")
        if path.exists():
            payload = json.loads(path.read_text())
            digest = hashlib.sha256(payload["text"].encode()).hexdigest()
            if digest != payload["sha256"]:
                raise InputError(
                    f"online cache checksum mismatch: {path}; remove this cache entry and retry"
                )
            if time.time() - payload["saved_at"] < 7 * 86400:
                self.transfers.append(
                    {
                        "url": url,
                        "sha256": digest,
                        "cached": True,
                        "retrieved_at": payload["retrieved_at"],
                    }
                )
                return payload["text"]
        request_params = (
            dict(params, tool="batch_primer_designer")
            if "eutils" in url
            else dict(params)
        )
        if self.email:
            request_params["email"] = self.email
        for attempt in range(4):
            wait = 0.4 - (time.monotonic() - self.last_request)
            if wait > 0:
                time.sleep(wait)
            self.last_request = time.monotonic()
            try:
                response = self.client.get(
                    url, params=request_params, stream=True, timeout=(15, 120)
                )
                try:
                    if (
                        response.status_code in {429, 500, 502, 503, 504}
                        and attempt < 3
                    ):
                        retry_after = response.headers.get("Retry-After", "0")
                        if not retry_after.isdigit():
                            raise InputError(
                                "NCBI asks to retry later; retry after the server indicated time"
                            )
                        delay = max(2**attempt, int(retry_after))
                        if delay > 60:
                            raise InputError(
                                f"NCBI asks to retry after {delay} seconds; stop and retry later"
                            )
                        response.close()
                        time.sleep(delay)
                        continue
                    response.raise_for_status()
                    chunks = []
                    size = 0
                    for chunk in response.iter_content(65536):
                        size += len(chunk)
                        if size > 20_000_000:
                            raise InputError(
                                "NCBI response exceeds 20 MB lightweight request limit"
                            )
                        chunks.append(chunk)
                    text = b"".join(chunks).decode("utf-8")
                    break
                finally:
                    response.close()
            except requests.RequestException as exc:
                if exc.response is not None and exc.response.status_code not in {
                    429,
                    500,
                    502,
                    503,
                    504,
                }:
                    raise InputError(f"NCBI request failed: {exc}") from exc
                if attempt >= 3:
                    raise InputError(f"NCBI request failed: {exc}") from exc
                time.sleep(2**attempt)
            except UnicodeError as exc:
                raise InputError("NCBI response is not valid UTF-8") from exc
        now = datetime.now(timezone.utc).isoformat()
        digest = hashlib.sha256(text.encode()).hexdigest()
        payload = {
            "text": text,
            "sha256": digest,
            "retrieved_at": now,
            "saved_at": time.time(),
        }
        self.cache_dir.mkdir(parents=True, exist_ok=True)
        with tempfile.NamedTemporaryFile(
            mode="w", dir=self.cache_dir, delete=False, encoding="utf-8"
        ) as tmp:
            json.dump(payload, tmp)
            temporary = Path(tmp.name)
        temporary.replace(path)
        self.transfers.append(
            {"url": url, "sha256": digest, "cached": False, "retrieved_at": now}
        )
        return text


class OnlineReference:
    def __init__(self, assembly, cache_dir, *, email=None, max_region_bp=200000):
        if not re.fullmatch(r"GC[AF]_\d+\.\d+", assembly):
            raise InputError(
                "online mode requires a versioned GCF_/GCA_ assembly accession"
            )
        if not 1 <= max_region_bp <= 1_000_000:
            raise InputError("max-region-bp must be 1..1000000")
        self.assembly = assembly
        self.max_region_bp = max_region_bp
        self.regions = []
        self.transport = NcbiTransport(cache_dir, email=email)
        try:
            reports = self.json(
                BASE + "/genome/accession/" + assembly + "/dataset_report", {}
            ).get("reports", [])
            matches = [r for r in reports if r.get("accession") == assembly]
            if len(matches) != 1:
                raise InputError(
                    "NCBI metadata does not match requested assembly version"
                )
            self.report = matches[0]
            organism = self.report.get("organism", {})
            self.taxid = organism.get("tax_id", organism.get("taxId"))
            if not str(self.taxid).isdigit():
                raise InputError("NCBI assembly metadata has no valid TaxID")
            self.sequences = {}
            self.aliases = {}
            self.accession_aliases = {}
            tokens = set()
            token = ""
            expected_count = None
            while True:
                page = self.json(
                    BASE + "/genome/accession/" + assembly + "/sequence_reports",
                    {"page_size": 1000, **({"page_token": token} if token else {})},
                )
                if expected_count is None and "total_count" in page:
                    expected_count = int(page["total_count"])
                for r in page.get("reports", []):
                    if r.get("assembly_accession") not in {None, assembly}:
                        raise InputError(
                            "sequence report contains another assembly version"
                        )
                    accessions = [
                        r[k]
                        for k in (
                            ("refseq_accession", "genbank_accession")
                            if assembly.startswith("GCF_")
                            else ("genbank_accession", "refseq_accession")
                        )
                        if r.get(k)
                    ]
                    if not accessions or not all(
                        VERSIONED.fullmatch(a) for a in accessions
                    ):
                        raise InputError(
                            "sequence report lacks versioned nucleotide accession"
                        )
                    canonical = accessions[0]
                    length = int(r.get("length", 0))
                    if length <= 0:
                        raise InputError("sequence report has invalid contig length")
                    if canonical in self.sequences:
                        raise InputError(
                            "duplicate contig in paginated sequence report"
                        )
                    self.sequences[canonical] = {
                        "length": length,
                        "aliases": accessions,
                        "chr_name": r.get("chr_name"),
                        "sequence_name": r.get("sequence_name"),
                        "ucsc_style_name": r.get("ucsc_style_name"),
                    }
                    if len(self.sequences) > 200000:
                        raise InputError(
                            "sequence report exceeds 200000 contig metadata limit"
                        )
                    for alias in accessions:
                        self.accession_aliases[alias] = canonical
                    for alias in accessions + [
                        str(r[k])
                        for k in ("chr_name", "ucsc_style_name", "sequence_name")
                        if r.get(k)
                    ]:
                        self.aliases.setdefault(alias, set()).add(canonical)
                token = page.get("next_page_token", "")
                if not token:
                    break
                if token in tokens or len(self.sequences) > 200000:
                    raise InputError(
                        "NCBI sequence-report pagination is incomplete or exceeds 200000 contigs"
                    )
                tokens.add(token)
            if expected_count is not None and expected_count != len(self.sequences):
                raise InputError("incomplete NCBI sequence report pagination")
            if not self.sequences:
                raise InputError(
                    "assembly contains no downloadable versioned sequence records"
                )
        except Exception:
            self.transport.close()
            raise

    def __enter__(self):
        return self

    def __exit__(self, *args):
        self.transport.close()

    def json(self, url, params):
        try:
            return json.loads(self.transport.fetch_text(url, params))
        except (ValueError, TypeError) as exc:
            raise InputError("NCBI returned invalid JSON") from exc

    def resolve_contig(self, name):
        candidates = self.aliases.get(str(name), set())
        if not candidates:
            raise InputError(f"contig is not in assembly {self.assembly}: {name}")
        if len(candidates) != 1:
            raise InputError(f"ambiguous contig alias: {name}; use versioned accession")
        return next(iter(candidates))

    def fetch_record(self, contig, start, end):
        contig = self.resolve_contig(contig)
        if not 0 <= start < end <= self.sequences[contig]["length"]:
            raise InputError("online interval outside assembly sequence")
        if end - start > self.max_region_bp:
            raise InputError(
                "target plus flanks exceeds max-region-bp; use a shorter BED interval or local reference"
            )
        text = self.transport.fetch_text(
            EUTILS + "efetch.fcgi",
            {
                "db": "nuccore",
                "id": contig,
                "seq_start": start + 1,
                "seq_stop": end,
                "strand": 1,
                "rettype": "gbwithparts",
                "retmode": "text",
            },
        )
        try:
            record = SeqIO.read(io.StringIO(text), "genbank")
        except ValueError as exc:
            raise InputError(
                "NCBI did not return one valid region GenBank record"
            ) from exc
        if len(record) != end - start:
            raise InputError("NCBI region length differs from requested interval")
        # EFetch may identify a subsequence as ACCESSION:start-end; require its base accession/version.
        returned = record.id.split(":", 1)[0]
        if returned not in self.sequences[contig]["aliases"]:
            raise InputError(
                "NCBI region accession/version differs from requested assembly sequence"
            )
        return contig, record

    def region(self, sid, contig, start, end, strand="+", flank=300):
        contig = self.resolve_contig(contig)
        if strand not in {"+", "-"} or flank < 0:
            raise InputError("strand must be +/- and flank nonnegative")
        if not 0 <= start < end <= self.sequences[contig]["length"]:
            raise InputError("target outside assembly sequence")
        left = max(0, start - flank)
        right = min(self.sequences[contig]["length"], end + flank)
        _, record = self.fetch_record(contig, left, right)
        seq = record.seq
        local_start = start - left + 1
        if strand == "-":
            seq = seq.reverse_complement()
            local_start = right - end + 1
        sequence = re.sub("[^ACGTN]", "N", str(seq).upper())
        features = []
        for feature in record.features:
            if feature.location is None:
                continue
            features.append(
                {
                    "type": feature.type,
                    "parts": [
                        {
                            "start_0based": left + int(p.start),
                            "end_exclusive": left + int(p.end),
                            "strand": p.strand,
                        }
                        for p in feature.location.parts
                    ],
                    "qualifiers": feature.qualifiers,
                    "location_in_fetched_record": str(feature.location),
                }
            )
        self.regions.append(
            {
                "sample_id": sid,
                "contig": contig,
                "start_0based": left,
                "end_exclusive": right,
                "target_start_0based": start,
                "target_end_exclusive": end,
                "strand": strand,
                "features": features,
                "template_sequence_5to3": sequence,
                "annotation_status": (
                    "regional_features_available"
                    if any(f["type"] != "source" for f in features)
                    else "no_regional_gene_features"
                ),
                "annotation_scope": "features returned by NCBI nuccore for this region; may be partial, absent or differ from a current Datasets annotation release",
            }
        )
        return Target(
            validate_record(TemplateRecord(sid, sequence, local_start, end - start)),
            contig,
            left,
            right,
            strand,
            self.assembly,
        )

    def gene(self, gene_id, flank=300):
        if not str(gene_id).isdigit():
            raise InputError(
                "online .txt input requires numeric NCBI Gene IDs; use BED for symbols or old assemblies"
            )
        result = (
            self.json(
                EUTILS + "esummary.fcgi",
                {"db": "gene", "id": str(gene_id), "retmode": "json"},
            )
            .get("result", {})
            .get(str(gene_id), {})
        )
        matches = []
        for item in result.get("genomicinfo", []):
            accession = item.get("chraccver")
            if accession in self.accession_aliases:
                a, b = int(item["chrstart"]), int(item["chrstop"])
                matches.append(
                    (
                        self.accession_aliases[accession],
                        min(a, b),
                        max(a, b) + 1,
                        "+" if a <= b else "-",
                    )
                )
        matches = list(dict.fromkeys(matches))
        if len(matches) != 1:
            raise InputError(
                f"Gene ID {gene_id} has no unique coordinates on requested assembly {self.assembly}; supply assembly-specific BED instead"
            )
        contig, a, b, strand = matches[0]
        return self.region(str(gene_id), contig, a, b, strand, flank)

    def read_bed(self, path, flank=300):
        from .inputs import open_text

        targets = []
        with open_text(path) as handle:
            for n, line in enumerate(handle, 1):
                if not line.strip() or line.startswith(("#", "track", "browser")):
                    continue
                f = line.rstrip().split("\t")
                if len(f) < 3:
                    raise InputError(f"BED line {n} requires at least 3 columns")
                try:
                    a, b = int(f[1]), int(f[2])
                except ValueError as exc:
                    raise InputError(f"BED line {n} has invalid coordinates") from exc
                targets.append(
                    self.region(
                        f[3] if len(f) > 3 else f"{f[0]}:{a}-{b}",
                        f[0],
                        a,
                        b,
                        f[5] if len(f) > 5 and f[5] != "." else "+",
                        flank,
                    )
                )
        if not targets:
            raise InputError("BED contains no targets")
        return targets

    def metadata(self):
        return {
            "provider": "NCBI",
            "assembly": self.assembly,
            "taxid": str(self.taxid),
            "assembly_report": self.report,
            "sequence_count": len(self.sequences),
            "sequences": self.sequences,
            "sequence_report_sha256": hashlib.sha256(
                json.dumps(self.sequences, sort_keys=True).encode()
            ).hexdigest(),
            "regions": self.regions,
            "transfers": self.transport.transfers,
            "max_region_bp": self.max_region_bp,
            "full_genome_downloaded": False,
        }
