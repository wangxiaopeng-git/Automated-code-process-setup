"""Bounded in-silico PCR; statuses describe the search, not experimental specificity."""

from __future__ import annotations
import bisect
import subprocess
import tempfile
from dataclasses import asdict, dataclass
from pathlib import Path
from Bio.Seq import Seq
from Bio import SeqIO
from .core import InputError

MAX_HITS = 5000
MAX_PRODUCTS = 10000


@dataclass(frozen=True)
class Hit:
    primer: str
    contig: str
    start: int
    end: int
    strand: str
    mismatches: tuple = ()  # query positions, 1-based from primer 5′ end


def exact_hits(sequences, left, right):
    hits = []
    for contig, record in sequences.items():
        sequence = str(record.seq if hasattr(record, "seq") else record).upper()
        for key, primer in [("L", left), ("R", right)]:
            for strand, query in [
                ("+", primer),
                ("-", str(Seq(primer).reverse_complement())),
            ]:
                start = sequence.find(query)
                while start >= 0:
                    hits.append(Hit(key, contig, start, start + len(query), strand))
                    if len(hits) >= MAX_HITS:
                        raise InputError(
                            "specificity hit limit reached; search incomplete"
                        )
                    start = sequence.find(query, start + 1)
    return hits


def pair_hits(hits, min_size, max_size):
    groups = {}
    for h in hits:
        groups.setdefault((h.contig, h.primer, h.strand), []).append(h)
    products = {}
    for (contig, key, strand), forward in groups.items():
        if strand != "+":
            continue
        other = "R" if key == "L" else "L"
        reverse = sorted(groups.get((contig, other, "-"), []), key=lambda h: h.end)
        ends = [h.end for h in reverse]
        for a in forward:
            for b in reverse[
                bisect.bisect_left(ends, a.start + min_size) : bisect.bisect_right(
                    ends, a.start + max_size
                )
            ]:
                if a.end > b.start:
                    continue
                identity = (contig, a.start, b.end, "+" if key == "L" else "-")
                products[identity] = {
                    "contig": contig,
                    "start": a.start,
                    "end": b.end,
                    "strand": identity[3],
                    "size_bp": b.end - a.start,
                    "forward_hit": asdict(a),
                    "reverse_hit": asdict(b),
                }
                if len(products) >= MAX_PRODUCTS:
                    raise InputError("amplicon limit reached; search incomplete")
    return [products[k] for k in sorted(products)]


def blast_hits(genome, left, right, *, max_mismatches=2):
    if not 0 <= max_mismatches <= 5:
        raise InputError("BLAST max_mismatches must be 0..5")
    # -subject requires no database writes and keeps the reference bundle immutable.
    with tempfile.TemporaryDirectory(prefix="primer-blast-") as temp:
        query = Path(temp) / "primers.fa"
        query.write_text(f">L\n{left}\n>R\n{right}\n")
        cmd = [
            "blastn",
            "-task",
            "blastn-short",
            "-query",
            str(query),
            "-subject",
            str(Path(genome).resolve()),
            "-dust",
            "no",
            "-soft_masking",
            "false",
            "-evalue",
            "1000",
            "-word_size",
            "4",
            "-max_target_seqs",
            "1000000",
            "-max_hsps",
            str(MAX_HITS),
            "-outfmt",
            "6 qseqid sseqid qlen length gapopen qstart qend sstart send qseq sseq",
        ]
        try:
            run = subprocess.run(
                cmd, capture_output=True, text=True, check=True, timeout=180
            )
        except FileNotFoundError as exc:
            raise InputError(
                "BLAST+ blastn is required; install ncbi-blast+ or choose exact"
            ) from exc
        except subprocess.CalledProcessError as exc:
            raise InputError(f"BLAST failed: {exc.stderr[-2000:]}") from exc
        except subprocess.TimeoutExpired as exc:
            raise InputError("BLAST timed out; specificity not completed") from exc
    hits = []
    lines = run.stdout.splitlines()
    # BLAST HSPs can trim terminal mismatches. Extend each ungapped seed to the
    # complete query span and independently re-count all substitutions.
    if len(lines) >= MAX_HITS:
        raise InputError("BLAST HSP limit reached; search incomplete")
    sequences = SeqIO.index(str(genome), "fasta")
    try:
        for line in lines:
            key, contig, qlen, length, gaps, qstart, qend, sstart, send, qseq, sseq = (
                line.split("\t")
            )
            if int(gaps):
                continue
            primer = left if key == "L" else right
            qlen = int(qlen)
            a, b = int(sstart), int(send)
            strand = "+" if a <= b else "-"
            if strand == "+":
                start = a - int(qstart)
                end = start + qlen
            else:
                end = a + int(qstart) - 1
                start = end - qlen
            record = sequences[contig]
            if start < 0 or end > len(record):
                continue
            segment = record.seq[start:end]
            if strand == "-":
                segment = segment.reverse_complement()
            mismatches = tuple(
                i + 1
                for i, (x, y) in enumerate(zip(primer.upper(), str(segment).upper()))
                if x != y
            )
            if len(mismatches) > max_mismatches:
                continue
            hits.append(Hit(key, contig, start, end, strand, mismatches))
    finally:
        sequences.close()
    return list(set(hits))


def check_pair(row, target, bundle, method, *, max_product=2000, max_mismatches=2):
    if method == "exact":
        if sum(bundle.manifest["contigs"].values()) > 5_000_000:
            raise InputError(
                "exact mode limited to 5 Mb references; use BLAST for larger genomes"
            )
        hits = exact_hits(
            bundle.sequences, row["left_sequence_5to3"], row["right_sequence_5to3"]
        )
    elif method == "blast":
        hits = blast_hits(
            bundle.path / "genome.fa",
            row["left_sequence_5to3"],
            row["right_sequence_5to3"],
            max_mismatches=max_mismatches,
        )
    else:
        raise InputError(f"unknown specificity method: {method}")
    products = pair_hits(
        hits,
        len(row["left_sequence_5to3"]) + len(row["right_sequence_5to3"]),
        max_product,
    )
    expected = None
    if target.contig:
        start, end = target.genomic_interval(
            row["left_binding_start_1based"], row["right_binding_end_1based"]
        )
        expected = (target.contig, start, end, target.strand)
    for p in products:
        p["is_expected"] = (
            (p["contig"], p["start"], p["end"], p["strand"]) == expected
            if expected
            else None
        )
    if expected and not any(p["is_expected"] for p in products):
        status = "expected_not_found"
    elif len(products) > 1:
        status = "multiple_products"
    elif len(products) == 1:
        status = (
            f"single_product_{method}"
            if expected
            else f"single_product_unassigned_{method}"
        )
    else:
        status = "no_product_found"
    return {
        "specificity_status": status,
        "specificity_method": method,
        "predicted_products": products,
        "specificity_hit_count": len(hits),
        "specificity_reference": bundle.manifest["assembly"],
        "specificity_limits": {
            "max_product_bp": max_product,
            "max_mismatches": 0 if method == "exact" else max_mismatches,
            "alignment": (
                "exact full-length"
                if method == "exact"
                else "BLAST heuristic seeds extended to full primer length, ungapped, word_size=4"
            ),
            "scope": "supplied contigs; hetero-primer products only; no experimental guarantee",
        },
    }
