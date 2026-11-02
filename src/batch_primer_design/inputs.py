"""Input adapters. Formats never change coordinate conventions implicitly."""

from __future__ import annotations
import csv
import gzip
from pathlib import Path
from Bio import SeqIO
from .core import (
    InputError,
    RecordIssue,
    TemplateRecord,
    validate_record,
    _optional_int,
)


def open_text(path):
    return (
        gzip.open(path, "rt", encoding="utf-8-sig")
        if str(path).lower().endswith(".gz")
        else open(path, encoding="utf-8-sig", newline="")
    )


def file_kind(path):
    name = Path(path).name.lower()
    return Path(name[:-3] if name.endswith(".gz") else name).suffix


def table_records(rows):
    rows = iter(rows)
    try:
        header = [str(x).strip() if x is not None else "" for x in next(rows)]
    except StopIteration:
        raise InputError("input table is empty")
    if len(set(header)) != len(header) or not {"id", "sequence"}.issubset(header):
        raise InputError("table requires unique id and sequence columns")
    for n, values in enumerate(rows, 2):
        if not any(v is not None and str(v).strip() for v in values):
            continue
        row = dict(zip(header, values))
        sid = str(row.get("id") or "").strip()
        try:
            yield validate_record(
                TemplateRecord(
                    sid,
                    str(row.get("sequence") or ""),
                    _optional_int(
                        (
                            str(row["target_start"])
                            if row.get("target_start") is not None
                            else None
                        ),
                        "target_start",
                        sid,
                    ),
                    _optional_int(
                        (
                            str(row["target_length"])
                            if row.get("target_length") is not None
                            else None
                        ),
                        "target_length",
                        sid,
                    ),
                )
            )
        except InputError as exc:
            yield RecordIssue(sid or f"line_{n}", str(exc))


def load_templates(path):
    path = Path(path)
    kind = file_kind(path)
    if kind in {".fa", ".fna", ".fasta"}:
        # Keep stricter FASTA validation, including empty headers and leading sequence.
        result = []
        sid = None
        parts = []

        def append():
            try:
                result.append(validate_record(TemplateRecord(sid, "".join(parts))))
            except InputError as exc:
                result.append(RecordIssue(sid, str(exc)))

        with open_text(path) as handle:
            for n, line in enumerate(handle, 1):
                line = line.strip()
                if not line:
                    continue
                if line.startswith(">"):
                    if sid is not None:
                        append()
                    if not line[1:].strip():
                        raise InputError(f"line {n}: empty FASTA header")
                    sid = line[1:].split()[0]
                    parts = []
                elif sid is None:
                    raise InputError(f"line {n}: sequence before FASTA header")
                else:
                    parts.append(line)
        if sid is not None:
            append()
    elif kind in {".csv", ".tsv"}:
        with open_text(path) as h:
            result = list(
                table_records(csv.reader(h, delimiter="\t" if kind == ".tsv" else ","))
            )
    elif kind == ".xlsx":
        from openpyxl import load_workbook

        book = load_workbook(path, read_only=True, data_only=False)
        try:
            rows = list(book.active.iter_rows(values_only=True))
            if any(
                isinstance(v, str) and v.startswith("=") for row in rows for v in row
            ):
                raise InputError("spreadsheet formulas are not accepted")
            result = list(table_records(rows))
        finally:
            book.close()
    elif kind in {".gb", ".gbk", ".genbank"}:
        result = []
        with open_text(path) as h:
            for r in SeqIO.parse(h, "genbank"):
                try:
                    result.append(validate_record(TemplateRecord(r.id, str(r.seq))))
                except InputError as exc:
                    result.append(RecordIssue(r.id, str(exc)))
    else:
        raise InputError(
            "supported templates: FASTA(.gz), CSV/TSV(.gz), XLSX, GenBank(.gz)"
        )
    if not result:
        raise InputError("input contains no records")
    return result
