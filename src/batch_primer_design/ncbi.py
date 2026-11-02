"""NCBI Datasets v2 client. Online transport is separate from local bundle use."""

from __future__ import annotations
import hashlib
import json
import re
import stat
import tempfile
import zipfile
from pathlib import Path, PurePosixPath
from urllib.parse import quote
import requests
from requests.adapters import HTTPAdapter
from urllib3.util.retry import Retry
from .core import InputError
from .reference import import_reference, sha256

BASE = "https://api.ncbi.nlm.nih.gov/datasets/v2"


def session():
    client = requests.Session()
    retry = Retry(
        total=3,
        backoff_factor=1,
        status_forcelist=[429, 500, 502, 503, 504],
        allowed_methods=["GET"],
        respect_retry_after_header=True,
    )
    client.mount("https://", HTTPAdapter(max_retries=retry))
    client.headers["User-Agent"] = "batch-primer-designer/0.2.0"
    return client


def search_ncbi(taxon, *, limit=20):
    if not str(taxon).strip() or not 1 <= limit <= 1000:
        raise InputError("taxon required; limit must be 1..1000")
    try:
        with session() as client:
            response = client.get(
                BASE
                + "/genome/taxon/"
                + quote(str(taxon), safe="")
                + "/dataset_report",
                params={"page_size": limit},
                timeout=(15, 120),
            )
            response.raise_for_status()
            return response.json()
    except (requests.RequestException, ValueError) as exc:
        raise InputError(f"NCBI search failed: {exc}") from exc


def safe_extract(archive, destination, max_bytes):
    with zipfile.ZipFile(archive) as z:
        total = 0
        for member in z.infolist():
            p = PurePosixPath(member.filename)
            mode = member.external_attr >> 16
            total += member.file_size
            if (
                p.is_absolute()
                or ".." in p.parts
                or "\\" in member.filename
                or ":" in member.filename
                or stat.S_ISLNK(mode)
            ):
                raise InputError("unsafe path in NCBI archive")
            if total > max_bytes:
                raise InputError("NCBI archive exceeds extraction limit")
        z.extractall(destination)
    verified = False
    covered = set()
    checksum = Path(destination) / "md5sum.txt"
    if checksum.exists():
        for line in checksum.read_text().splitlines():
            if not line.strip():
                continue
            try:
                digest, name = line.split(maxsplit=1)
            except ValueError as exc:
                raise InputError("malformed NCBI checksum manifest") from exc
            name = name.lstrip("*")
            p = PurePosixPath(name)
            if p.is_absolute() or ".." in p.parts or "\\" in name:
                raise InputError("unsafe checksum path")
            file = Path(destination) / name
            h = (
                hashlib.md5()
            )  # NCBI uses MD5 for transfer verification, not authentication.
            with file.open("rb") as f:
                for block in iter(lambda: f.read(1024 * 1024), b""):
                    h.update(block)
            if h.hexdigest().lower() != digest.lower():
                raise InputError(f"NCBI checksum mismatch: {name}")
            covered.add(file.resolve())
        required = {
            p.resolve()
            for p in Path(destination).rglob("*")
            if p.is_file() and p != checksum
        }
        verified = bool(covered) and bool(required) and required.issubset(covered)
    return verified


def download_ncbi(
    accession,
    output,
    *,
    taxid,
    max_bytes=10_000_000_000,
    annotation_release="NCBI download",
):
    if not re.fullmatch(r"GC[AF]_\d+\.\d+", accession):
        raise InputError("use a versioned GCF_/GCA_ assembly accession")
    output = Path(output).absolute()
    if output.exists():
        raise InputError(
            "reference output already exists; reuse it offline or choose another directory"
        )
    output.parent.mkdir(parents=True, exist_ok=True)
    url = BASE + "/genome/accession/" + accession + "/download"
    with tempfile.TemporaryDirectory(prefix=".ncbi-", dir=output.parent) as temp:
        temp = Path(temp)
        archive = temp / "download.zip"
        unpacked = temp / "unpacked"
        try:
            with session() as client:
                response = client.get(
                    url,
                    params={
                        "include_annotation_type": [
                            "GENOME_FASTA",
                            "GENOME_GFF",
                            "RNA_FASTA",
                        ],
                        "filename": "dataset.zip",
                    },
                    stream=True,
                    timeout=(15, 180),
                )
                response.raise_for_status()
                size = 0
                try:
                    with archive.open("wb") as h:
                        for chunk in response.iter_content(1024 * 1024):
                            size += len(chunk)
                            if size > max_bytes:
                                raise InputError("NCBI download exceeds size limit")
                            h.write(chunk)
                finally:
                    response.close()
            provider_verified = safe_extract(archive, unpacked, max_bytes * 4)
        except (requests.RequestException, zipfile.BadZipFile, OSError) as exc:
            raise InputError(f"NCBI download failed: {exc}") from exc
        report_path = unpacked / "ncbi_dataset" / "data" / "assembly_data_report.jsonl"
        if not report_path.is_file():
            raise InputError(
                "NCBI archive missing assembly metadata; cannot verify taxonomy"
            )
        metadata = [
            json.loads(line)
            for line in report_path.read_text().splitlines()
            if line.strip()
        ]
        matched = [r for r in metadata if r.get("accession") == accession]
        if len(matched) != 1:
            raise InputError(
                "NCBI assembly metadata does not match requested accession"
            )
        organism = matched[0].get("organism", {})
        actual_taxid = organism.get("taxId", organism.get("tax_id"))
        if str(actual_taxid) != str(taxid):
            raise InputError(
                f"NCBI taxid mismatch: assembly taxid is {actual_taxid}, requested {taxid}"
            )
        data = unpacked / "ncbi_dataset" / "data" / accession
        genomes = list(data.glob("*.fna"))
        genomes = [p for p in genomes if p.name.endswith("genomic.fna")]
        if len(genomes) != 1:
            raise InputError(
                "NCBI archive must contain one genomic FASTA for requested accession"
            )
        annotations = list(data.glob("*.gff"))
        if len(annotations) > 1:
            raise InputError("ambiguous NCBI annotation files")
        bundle = temp / "bundle"
        manifest = import_reference(
            genomes[0],
            annotations[0] if annotations else None,
            bundle,
            taxid=taxid,
            assembly=accession,
            source=url,
            annotation_release=annotation_release,
        )
        # Keep transcript FASTA and metadata if supplied, but do not claim spliced design support.
        for file in data.iterdir():
            if (
                file.is_file()
                and file.name.endswith((".fna", ".jsonl"))
                and file not in genomes
            ):
                import shutil

                shutil.copy2(file, bundle / file.name)
                manifest["files"][file.name] = sha256(bundle / file.name)
        info = {
            "url": url,
            "accession": accession,
            "archive_sha256": sha256(archive),
            "ncbi_md5_manifest_verified": provider_verified,
            "taxid_source": "verified against downloaded NCBI assembly metadata",
            "annotation_available": bool(annotations),
        }
        (bundle / "download.json").write_text(json.dumps(info, indent=2))
        manifest["files"]["download.json"] = sha256(bundle / "download.json")
        (bundle / "manifest.json").write_text(json.dumps(manifest, indent=2))
        bundle.rename(output)
    return manifest
