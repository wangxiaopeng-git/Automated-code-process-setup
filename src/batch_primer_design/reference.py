"""Local reference bundles with validated annotations and explicit strand mapping."""

from __future__ import annotations
import hashlib
import json
import re
import shutil
import sqlite3
import tempfile
from dataclasses import dataclass
from datetime import datetime, timezone
from pathlib import Path
from urllib.parse import unquote
from Bio import SeqIO
from Bio.Seq import Seq
from .core import InputError, TemplateRecord, validate_record
from .inputs import open_text


def sha256(path):
    h = hashlib.sha256()
    with open(path, "rb") as f:
        for block in iter(lambda: f.read(1024 * 1024), b""):
            h.update(block)
    return h.hexdigest()


@dataclass
class Target:
    record: TemplateRecord
    contig: str = ""
    start: int = 0
    end: int = 0
    strand: str = "+"
    reference: str = ""

    def genomic_interval(self, start_1based, end_1based):
        if not self.contig:
            raise InputError("template has no genomic mapping")
        if not (1 <= start_1based <= end_1based <= len(self.record.sequence)):
            raise InputError("primer coordinate outside template")
        if self.strand == "+":
            return self.start + start_1based - 1, self.start + end_1based
        return self.end - end_1based, self.end - start_1based + 1


def parse_attributes(text):
    if "=" in text:
        return {
            unquote(k): unquote(v)
            for item in text.split(";")
            if "=" in item
            for k, v in [item.split("=", 1)]
        }
    return dict(re.findall(r'(\S+)\s+"([^"]*)"', text))


def build_annotation(path, db, lengths):
    db.execute(
        "CREATE TABLE features (id TEXT, label TEXT, kind TEXT, contig TEXT, start INT, end INT, strand TEXT)"
    )
    grouped = {}
    count = 0
    with open_text(path) as h:
        for n, line in enumerate(h, 1):
            if line.startswith("##FASTA"):
                break
            if not line.strip() or line.startswith("#"):
                continue
            fields = line.rstrip("\n").split("\t")
            if len(fields) != 9:
                raise InputError(f"annotation line {n}: expected 9 columns")
            contig, _, kind, a, b, _, strand, _, attrs = fields
            try:
                start, end = int(a) - 1, int(b)
            except ValueError as exc:
                raise InputError(f"annotation line {n}: invalid coordinates") from exc
            if contig not in lengths or not (0 <= start < end <= lengths[contig]):
                raise InputError(
                    f"annotation line {n}: contig mismatch or coordinates outside genome"
                )
            if strand not in {"+", "-", ".", "?"}:
                raise InputError(f"annotation line {n}: invalid strand")
            attrs = parse_attributes(attrs)
            count += 1
            # GTF often contains only exons/CDS: aggregate genomic spans, never splice silently.
            if "gene_id" in attrs:
                for key, k in [("gene_id", "gene"), ("transcript_id", "transcript")]:
                    if key not in attrs:
                        continue
                    identity = (attrs[key], k)
                    value = [
                        attrs[key],
                        (
                            attrs.get("gene_name", "")
                            if k == "gene"
                            else attrs.get("transcript_name", "")
                        ),
                        k,
                        contig,
                        start,
                        end,
                        strand,
                    ]
                    if identity in grouped:
                        old = grouped[identity]
                        if old[3] != contig or old[6] != strand:
                            raise InputError(f"conflicting annotation ID: {attrs[key]}")
                        old[4] = min(old[4], start)
                        old[5] = max(old[5], end)
                        if not old[1]:
                            old[1] = value[1]
                    else:
                        grouped[identity] = value
            elif kind.lower() in {
                "gene",
                "mrna",
                "transcript",
                "ncrna",
                "rrna",
                "trna",
            }:
                fid = attrs.get("ID")
                if not fid:
                    continue
                identity = (fid, kind)
                if any(k[0] == fid for k in grouped):
                    raise InputError(f"duplicate annotation ID: {fid}")
                grouped[identity] = [
                    fid,
                    attrs.get("Name", ""),
                    kind,
                    contig,
                    start,
                    end,
                    strand,
                ]
    if not count:
        raise InputError("annotation contains no features")
    db.executemany("INSERT INTO features VALUES (?,?,?,?,?,?,?)", grouped.values())
    db.execute("CREATE INDEX by_id ON features(id)")
    db.execute("CREATE INDEX by_label ON features(label)")
    db.commit()


def import_reference(
    genome,
    annotation,
    output,
    *,
    taxid,
    assembly,
    source="local",
    annotation_release="unspecified",
):
    output = Path(output).absolute()
    if output.exists():
        raise InputError(f"bundle already exists: {output}")
    if taxid is not None and (not str(taxid).isdigit() or int(taxid) <= 0):
        raise InputError("taxid must be a positive NCBI taxonomy ID")
    if not assembly.strip():
        raise InputError("assembly identifier is required")
    output.parent.mkdir(parents=True, exist_ok=True)
    with tempfile.TemporaryDirectory(prefix=".reference-", dir=output.parent) as tmp:
        stage = Path(tmp)
        fasta = stage / "genome.fa"
        lengths = {}
        with open_text(genome) as src, fasta.open("w") as dest:
            for rec in SeqIO.parse(src, "fasta"):
                if rec.id in lengths:
                    raise InputError(f"duplicate contig: {rec.id}")
                seq = str(rec.seq).upper()
                if not seq or re.search("[^ACGTRYSWKMBDHVN]", seq):
                    raise InputError(f"invalid reference bases: {rec.id}")
                lengths[rec.id] = len(seq)
                rec.seq = Seq(seq)
                SeqIO.write(rec, dest, "fasta")
        if not lengths:
            raise InputError("genome has no sequences")
        if annotation:
            with (
                open_text(annotation) as src,
                (stage / "annotation.txt").open("w") as dst,
            ):
                shutil.copyfileobj(src, dst)
            with sqlite3.connect(stage / "annotation.sqlite") as db:
                build_annotation(stage / "annotation.txt", db, lengths)
        manifest = {
            "schema_version": 1,
            "taxid": str(taxid) if taxid is not None else None,
            "assembly": assembly,
            "annotation_release": annotation_release,
            "source": source,
            "created_at": datetime.now(timezone.utc).isoformat(),
            "contigs": lengths,
            "files": {p.name: sha256(p) for p in stage.iterdir() if p.is_file()},
        }
        (stage / "manifest.json").write_text(
            json.dumps(manifest, indent=2), encoding="utf-8"
        )
        # Rename only after full validation. No incomplete bundle can appear at output.
        stage.rename(output)
    return manifest


class ReferenceBundle:
    def __init__(self, path):
        self.path = Path(path)
        self.manifest = json.loads((self.path / "manifest.json").read_text())
        self.sequences = SeqIO.index(str(self.path / "genome.fa"), "fasta")
        self.db = None
        if (self.path / "annotation.sqlite").exists():
            self.db = sqlite3.connect(
                f'file:{(self.path/"annotation.sqlite").resolve()}?mode=ro', uri=True
            )

    def __enter__(self):
        return self

    def __exit__(self, *args):
        self.sequences.close()
        if self.db:
            self.db.close()

    def verify(self):
        for name, digest in self.manifest["files"].items():
            if Path(name).name != name or sha256(self.path / name) != digest:
                raise InputError(f"reference checksum mismatch: {name}")

    def target(self, sid, contig, start, end, strand="+", flank=300):
        if strand not in {"+", "-"}:
            raise InputError("target strand must be + or -")
        if flank < 0:
            raise InputError("flank must be nonnegative")
        if contig not in self.sequences:
            raise InputError(f"unknown contig: {contig}")
        size = self.manifest["contigs"][contig]
        if not 0 <= start < end <= size:
            raise InputError("target interval outside reference")
        left = max(0, start - flank)
        right = min(size, end + flank)
        seq = self.sequences[contig].seq[left:right]
        target_start = start - left + 1
        if strand == "-":
            seq = seq.reverse_complement()
            target_start = right - end + 1
        # Ambiguous reference bases become N in design templates, without claiming degeneracy.
        sequence = re.sub("[^ACGTN]", "N", str(seq).upper())
        rec = validate_record(TemplateRecord(sid, sequence, target_start, end - start))
        return Target(rec, contig, left, right, strand, self.manifest["assembly"])

    def resolve_gene(self, identifier, flank=300):
        if not self.db:
            raise InputError("reference has no annotation")
        rows = self.db.execute(
            "SELECT * FROM features WHERE id=?", (identifier,)
        ).fetchall()
        if not rows:
            rows = self.db.execute(
                "SELECT * FROM features WHERE label=?", (identifier,)
            ).fetchall()
        if not rows:
            raise InputError(f"gene/transcript not found: {identifier}")
        if len(rows) != 1:
            raise InputError(
                f"ambiguous gene/transcript: {identifier}; use exact unique ID"
            )
        fid, label, kind, contig, start, end, strand = rows[0]
        if strand not in {"+", "-"}:
            raise InputError(f"unknown strand for {identifier}")
        return self.target(identifier, contig, start, end, strand, flank)

    def read_bed(self, path, flank=300):
        result = []
        with open_text(path) as h:
            for n, line in enumerate(h, 1):
                if not line.strip() or line.startswith(("#", "track", "browser")):
                    continue
                f = line.rstrip().split("\t")
                if len(f) < 3:
                    raise InputError(f"BED line {n}: requires at least 3 columns")
                try:
                    a, b = int(f[1]), int(f[2])
                except ValueError as exc:
                    raise InputError(f"BED line {n}: invalid coordinates") from exc
                strand = f[5] if len(f) > 5 and f[5] != "." else "+"
                result.append(
                    self.target(
                        f[3] if len(f) > 3 else f"{f[0]}:{a}-{b}",
                        f[0],
                        a,
                        b,
                        strand,
                        flank,
                    )
                )
        if not result:
            raise InputError("BED contains no targets")
        return result
