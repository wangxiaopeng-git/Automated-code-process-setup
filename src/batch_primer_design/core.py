from __future__ import annotations

import csv
import math
import re
from dataclasses import dataclass
from pathlib import Path
from typing import Any, Callable, Dict, Iterable, Iterator, List, Optional, Sequence, Union


DNA_RE = re.compile(r"^[ACGTN]+$")


class InputError(ValueError):
    """Raised when an input record is invalid."""


@dataclass(frozen=True)
class TemplateRecord:
    sample_id: str
    sequence: str
    target_start: Optional[int] = None  # user-facing 1-based coordinate
    target_length: Optional[int] = None


@dataclass(frozen=True)
class RecordIssue:
    sample_id: str
    message: str


@dataclass(frozen=True)
class DesignConfig:
    num_return: int = 5
    primer_min_size: int = 18
    primer_opt_size: int = 20
    primer_max_size: int = 25
    primer_min_tm: float = 57.0
    primer_opt_tm: float = 60.0
    primer_max_tm: float = 63.0
    primer_min_gc: float = 35.0
    primer_max_gc: float = 65.0
    product_min_size: int = 100
    product_max_size: int = 300
    max_poly_x: int = 4
    max_ns: int = 0
    salt_monovalent_mm: float = 50.0
    salt_divalent_mm: float = 1.5
    dntp_mm: float = 0.6
    primer_dna_nm: float = 50.0
    pick_anyway: bool = False

    def validate(self) -> None:
        for name in ("primer_min_tm", "primer_opt_tm", "primer_max_tm", "primer_min_gc", "primer_max_gc", "salt_monovalent_mm", "salt_divalent_mm", "dntp_mm", "primer_dna_nm"):
            if not math.isfinite(getattr(self, name)):
                raise InputError(f"{name} must be finite")
        for name in ("salt_monovalent_mm", "salt_divalent_mm", "dntp_mm"):
            if getattr(self, name) < 0:
                raise InputError(f"{name} must be nonnegative")
        if self.primer_dna_nm <= 0:
            raise InputError("primer_dna_nm must be positive")
        if self.num_return > 100:
            raise InputError("num_return must be <= 100")
        if self.num_return < 1:
            raise InputError("num_return must be >= 1")
        if not (1 <= self.primer_min_size <= self.primer_opt_size <= self.primer_max_size):
            raise InputError("primer sizes must satisfy 1 <= min <= opt <= max")
        if not (self.primer_min_tm <= self.primer_opt_tm <= self.primer_max_tm):
            raise InputError("primer Tm values must satisfy min <= opt <= max")
        if not (0 <= self.primer_min_gc <= self.primer_max_gc <= 100):
            raise InputError("primer GC values must satisfy 0 <= min <= max <= 100")
        if not (1 <= self.product_min_size <= self.product_max_size):
            raise InputError("product sizes must satisfy 1 <= min <= max")
        if self.max_poly_x < 1 or self.max_ns < 0:
            raise InputError("max_poly_x must be >= 1 and max_ns must be >= 0")

    def primer3_args(self) -> Dict[str, Any]:
        self.validate()
        return {
            "PRIMER_TASK": "generic",
            "PRIMER_PICK_LEFT_PRIMER": 1,
            "PRIMER_PICK_RIGHT_PRIMER": 1,
            "PRIMER_PICK_INTERNAL_OLIGO": 0,
            "PRIMER_NUM_RETURN": self.num_return,
            "PRIMER_MIN_SIZE": self.primer_min_size,
            "PRIMER_OPT_SIZE": self.primer_opt_size,
            "PRIMER_MAX_SIZE": self.primer_max_size,
            "PRIMER_MIN_TM": self.primer_min_tm,
            "PRIMER_OPT_TM": self.primer_opt_tm,
            "PRIMER_MAX_TM": self.primer_max_tm,
            "PRIMER_MIN_GC": self.primer_min_gc,
            "PRIMER_MAX_GC": self.primer_max_gc,
            "PRIMER_PRODUCT_SIZE_RANGE": [
                [self.product_min_size, self.product_max_size]
            ],
            "PRIMER_MAX_POLY_X": self.max_poly_x,
            "PRIMER_MAX_NS_ACCEPTED": self.max_ns,
            "PRIMER_SALT_MONOVALENT": self.salt_monovalent_mm,
            "PRIMER_SALT_DIVALENT": self.salt_divalent_mm,
            "PRIMER_DNTP_CONC": self.dntp_mm,
            "PRIMER_DNA_CONC": self.primer_dna_nm,
            "PRIMER_EXPLAIN_FLAG": 1,
            "PRIMER_PICK_ANYWAY": int(self.pick_anyway),
        }


RESULT_FIELDS = [
    "sample_id",
    "pair_rank",
    "left_sequence_5to3",
    "right_sequence_5to3",
    "left_tm_c",
    "right_tm_c",
    "left_gc_percent",
    "right_gc_percent",
    "left_binding_start_1based",
    "left_binding_end_1based",
    "right_binding_start_1based",
    "right_binding_end_1based",
    "product_size_bp",
    "pair_penalty",
    "status",
    "message",
]


def normalize_sequence(sequence: str) -> str:
    return re.sub(r"\s+", "", sequence).upper()


def validate_record(record: TemplateRecord) -> TemplateRecord:
    sample_id = record.sample_id.strip()
    sequence = normalize_sequence(record.sequence)
    if not sample_id:
        raise InputError("sample id is empty")
    if not sequence:
        raise InputError(f"{sample_id}: sequence is empty")
    if not DNA_RE.fullmatch(sequence):
        invalid = "".join(sorted(set(sequence) - set("ACGTN")))
        raise InputError(
            f"{sample_id}: sequence contains unsupported base(s): {invalid}; "
            "allowed bases are A/C/G/T/N"
        )

    one_missing = (record.target_start is None) != (record.target_length is None)
    if one_missing:
        raise InputError(
            f"{sample_id}: target_start and target_length must be supplied together"
        )
    if record.target_start is not None and record.target_length is not None:
        if record.target_start < 1 or record.target_length < 1:
            raise InputError(f"{sample_id}: target coordinates must be positive")
        if record.target_start + record.target_length - 1 > len(sequence):
            raise InputError(f"{sample_id}: target region extends beyond the template")

    return TemplateRecord(
        sample_id=sample_id,
        sequence=sequence,
        target_start=record.target_start,
        target_length=record.target_length,
    )


def _optional_int(value: Optional[str], field: str, sample_id: str) -> Optional[int]:
    if value is None or not value.strip():
        return None
    try:
        return int(value)
    except ValueError as exc:
        raise InputError(f"{sample_id}: {field} must be an integer") from exc


def read_delimited(path: Path) -> Iterator[Union[TemplateRecord, RecordIssue]]:
    delimiter = "\t" if path.suffix.lower() == ".tsv" else ","
    with path.open("r", encoding="utf-8-sig", newline="") as handle:
        reader = csv.DictReader(handle, delimiter=delimiter)
        if reader.fieldnames is None:
            raise InputError("input table has no header")
        required = {"id", "sequence"}
        missing = required - set(reader.fieldnames)
        if missing:
            raise InputError("input table is missing column(s): " + ", ".join(sorted(missing)))
        for line_number, row in enumerate(reader, start=2):
            sample_id = (row.get("id") or "").strip() or f"line_{line_number}"
            try:
                record = TemplateRecord(
                    sample_id=sample_id,
                    sequence=row.get("sequence") or "",
                    target_start=_optional_int(
                        row.get("target_start"), "target_start", sample_id
                    ),
                    target_length=_optional_int(
                        row.get("target_length"), "target_length", sample_id
                    ),
                )
                yield validate_record(record)
            except InputError as exc:
                yield RecordIssue(sample_id, str(exc))


def read_fasta(path: Path) -> Iterator[Union[TemplateRecord, RecordIssue]]:
    header: Optional[str] = None
    sequence_lines: List[str] = []
    with path.open("r", encoding="utf-8-sig") as handle:
        for line_number, raw_line in enumerate(handle, start=1):
            line = raw_line.strip()
            if not line:
                continue
            if line.startswith(">"):
                if header is not None:
                    try:
                        yield validate_record(
                            TemplateRecord(header, "".join(sequence_lines))
                        )
                    except InputError as exc:
                        yield RecordIssue(header, str(exc))
                title = line[1:].strip()
                if not title:
                    raise InputError(f"line {line_number}: FASTA header is empty")
                header = title.split(maxsplit=1)[0]
                sequence_lines = []
            else:
                if header is None:
                    raise InputError(f"line {line_number}: sequence appears before FASTA header")
                sequence_lines.append(line)
    if header is None:
        raise InputError("FASTA input contains no records")
    try:
        yield validate_record(TemplateRecord(header, "".join(sequence_lines)))
    except InputError as exc:
        yield RecordIssue(header, str(exc))


def read_records(path: Path) -> Iterator[Union[TemplateRecord, RecordIssue]]:
    suffix = path.suffix.lower()
    if suffix in {".fa", ".fasta", ".fna"}:
        yield from read_fasta(path)
    elif suffix in {".csv", ".tsv"}:
        yield from read_delimited(path)
    else:
        raise InputError("input must be .fasta/.fa/.fna, .csv, or .tsv")


def _blank_result(sample_id: str, status: str, message: str) -> Dict[str, Any]:
    row = {field: "" for field in RESULT_FIELDS}
    row.update({"sample_id": sample_id, "status": status, "message": message})
    return row


def _message_from_result(result: Dict[str, Any]) -> str:
    parts = []
    for key in ("PRIMER_PAIR_EXPLAIN", "PRIMER_LEFT_EXPLAIN", "PRIMER_RIGHT_EXPLAIN"):
        if result.get(key):
            parts.append(f"{key}={result[key]}")
    return "; ".join(parts) or "Primer3 returned no primer pair"


def design_record(
    record: TemplateRecord,
    config: DesignConfig,
    engine: Callable[..., Dict[str, Any]],
) -> List[Dict[str, Any]]:
    seq_args: Dict[str, Any] = {
        "SEQUENCE_ID": record.sample_id,
        "SEQUENCE_TEMPLATE": record.sequence,
    }
    if record.target_start is not None and record.target_length is not None:
        seq_args["SEQUENCE_TARGET"] = [record.target_start - 1, record.target_length]

    try:
        result = engine(seq_args=seq_args, global_args=config.primer3_args())
    except Exception as exc:  # Primer3 exposes several exception types across versions.
        return [_blank_result(record.sample_id, "error", str(exc))]

    if result.get("PRIMER_ERROR"):
        return [_blank_result(record.sample_id, "error", str(result["PRIMER_ERROR"]))]
    pair_count = int(result.get("PRIMER_PAIR_NUM_RETURNED", 0) or 0)
    if pair_count == 0:
        return [_blank_result(record.sample_id, "no_pair", _message_from_result(result))]

    rows: List[Dict[str, Any]] = []
    for index in range(pair_count):
        left_pos, left_len = result[f"PRIMER_LEFT_{index}"]
        right_pos, right_len = result[f"PRIMER_RIGHT_{index}"]
        rows.append(
            {
                "sample_id": record.sample_id,
                "pair_rank": index + 1,
                "left_sequence_5to3": result[f"PRIMER_LEFT_{index}_SEQUENCE"],
                "right_sequence_5to3": result[f"PRIMER_RIGHT_{index}_SEQUENCE"],
                "left_tm_c": round(float(result[f"PRIMER_LEFT_{index}_TM"]), 3),
                "right_tm_c": round(float(result[f"PRIMER_RIGHT_{index}_TM"]), 3),
                "left_gc_percent": round(
                    float(result[f"PRIMER_LEFT_{index}_GC_PERCENT"]), 3
                ),
                "right_gc_percent": round(
                    float(result[f"PRIMER_RIGHT_{index}_GC_PERCENT"]), 3
                ),
                "left_binding_start_1based": left_pos + 1,
                "left_binding_end_1based": left_pos + left_len,
                # Primer3 reports the rightmost template base (primer 5′ end).
                "right_binding_start_1based": right_pos - right_len + 2,
                "right_binding_end_1based": right_pos + 1,
                "product_size_bp": result[f"PRIMER_PAIR_{index}_PRODUCT_SIZE"],
                "pair_penalty": round(
                    float(result.get(f"PRIMER_PAIR_{index}_PENALTY", 0.0)), 6
                ),
                "status": "ok",
                "message": str(result.get("PRIMER_WARNING", "")),
            }
        )
    return rows


def design_batch(
    records: Iterable[Union[TemplateRecord, RecordIssue]],
    config: DesignConfig,
    engine: Callable[..., Dict[str, Any]],
) -> Iterator[Dict[str, Any]]:
    seen = set()
    for record in records:
        if isinstance(record, RecordIssue):
            yield _blank_result(record.sample_id, "error", record.message)
            continue
        if record.sample_id in seen:
            yield _blank_result(record.sample_id, "error", "duplicate sample id")
            continue
        seen.add(record.sample_id)
        yield from design_record(record, config, engine)


def write_results(rows: Sequence[Dict[str, Any]], output_path: Path) -> None:
    output_path.parent.mkdir(parents=True, exist_ok=True)
    with output_path.open("w", encoding="utf-8-sig", newline="") as handle:
        writer = csv.DictWriter(handle, fieldnames=RESULT_FIELDS)
        writer.writeheader()
        writer.writerows(rows)

