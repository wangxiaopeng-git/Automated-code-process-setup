"""Prepare raw reference files transparently; reuse validated, content-keyed caches."""

from __future__ import annotations
import hashlib
import json
from pathlib import Path
from .core import InputError
from .reference import import_reference, sha256


def prepare_reference(reference, annotation=None, *, cache_dir):
    reference = Path(reference).resolve()
    annotation = Path(annotation).resolve() if annotation else None
    if reference.is_dir():
        if annotation:
            raise InputError(
                "Existing reference bundle already contains its annotation; use --reference genome FASTA with --gtf/--gff"
            )
        if not (reference / "manifest.json").is_file():
            raise InputError(
                "--reference directory must be an imported reference bundle; otherwise provide the genome FASTA file"
            )
        return reference
    if not reference.is_file():
        raise InputError(f"Reference genome file does not exist: {reference}")
    if annotation and not annotation.is_file():
        raise InputError(f"Annotation file does not exist: {annotation}")
    inputs = {
        "genome_sha256": sha256(reference),
        "annotation_sha256": sha256(annotation) if annotation else None,
    }
    # Include cache schema so future normalization changes cannot reuse incompatible files.
    key = hashlib.sha256(
        json.dumps({"cache_schema": 1, **inputs}, sort_keys=True).encode()
    ).hexdigest()
    destination = Path(cache_dir) / key
    if not destination.exists():
        manifest = import_reference(
            reference,
            annotation,
            destination,
            taxid=None,
            assembly="local-sha256:" + inputs["genome_sha256"],
            source="local raw files; taxonomy unspecified",
        )
        manifest["original_inputs"] = inputs
        (destination / "manifest.json").write_text(
            json.dumps(manifest, indent=2), encoding="utf-8"
        )
        print(f"已自动准备本地参考缓存：{destination}")
    else:
        print(f"复用本地参考缓存：{destination}")
    return destination
