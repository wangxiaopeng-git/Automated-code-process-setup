"""Advanced workflow: separate design, mapping, thermodynamics and specificity."""

from __future__ import annotations
import json
import platform
import tempfile
from dataclasses import asdict
from datetime import datetime, timezone
from importlib.metadata import version
from pathlib import Path
from primer3 import bindings
from .core import InputError, RecordIssue, _blank_result, design_record, validate_record
from .exporters import export_report
from .reference import sha256
from .specificity import check_pair


def thermo(row, config):
    kwargs = {
        "mv_conc": config.salt_monovalent_mm,
        "dv_conc": config.salt_divalent_mm,
        "dntp_conc": config.dntp_mm,
        "dna_conc": config.primer_dna_nm,
    }
    left, right = row["left_sequence_5to3"], row["right_sequence_5to3"]
    values = {}
    for name, fn, args in [
        ("left_hairpin", bindings.calc_hairpin, (left,)),
        ("right_hairpin", bindings.calc_hairpin, (right,)),
        ("left_homodimer", bindings.calc_homodimer, (left,)),
        ("right_homodimer", bindings.calc_homodimer, (right,)),
        ("heterodimer", bindings.calc_heterodimer, (left, right)),
    ]:
        result = fn(*args, **kwargs)
        values[name] = {
            "structure_found": result.structure_found,
            "tm_c": result.tm,
            "dg_cal_per_mol": result.dg,
            "calculation_temperature_c": 37.0,
        }
    return values


def run_design(
    targets,
    config,
    output,
    *,
    bundle=None,
    specificity="none",
    max_product=2000,
    max_mismatches=2,
    provenance=None,
):
    config.validate()
    targets = list(targets)
    output = Path(output).absolute()
    if not targets:
        raise InputError("input contains no targets")
    if output.exists():
        raise InputError("output exists; choose a new run directory")
    if specificity != "none" and bundle is None:
        raise InputError("specificity requires a reference bundle")
    if specificity not in {"none", "exact", "blast"}:
        raise InputError("specificity must be none/exact/blast")
    if max_product < config.product_max_size:
        raise InputError("specificity max_product must cover design product maximum")
    rows = []
    seen = set()
    success = 0
    checked_errors = 0
    for n, target in enumerate(targets, 1):
        record = target if isinstance(target, RecordIssue) else target.record
        if isinstance(record, RecordIssue):
            current = [_blank_result(record.sample_id, "error", record.message)]
        elif record.sample_id in seen:
            current = [_blank_result(record.sample_id, "error", "duplicate sample id")]
        else:
            seen.add(record.sample_id)
            try:
                record = validate_record(record)
            except InputError as exc:
                current = [_blank_result(record.sample_id, "error", str(exc))]
            else:
                current = design_record(record, config, bindings.design_primers)
        success += int(any(r["status"] == "ok" for r in current))
        for row in current:
            row.update(
                record_id=n, specificity_status="not_run", specificity_method="none"
            )
            if row["status"] == "ok":
                a, b = row["left_binding_start_1based"], row["right_binding_end_1based"]
                row["amplicon_sequence_5to3"] = record.sequence[a - 1 : b]
                row["tm_difference_c"] = round(
                    abs(row["left_tm_c"] - row["right_tm_c"]), 3
                )
                row["left_length_nt"] = len(row["left_sequence_5to3"])
                row["right_length_nt"] = len(row["right_sequence_5to3"])
                row["thermodynamics"] = thermo(row, config)
                if target.contig:
                    row.update(
                        contig=target.contig,
                        assembly=target.reference,
                        template_strand=target.strand,
                    )
                    for side in ("left", "right"):
                        start, end = target.genomic_interval(
                            row[side + "_binding_start_1based"],
                            row[side + "_binding_end_1based"],
                        )
                        strand = (
                            target.strand
                            if side == "left"
                            else ("-" if target.strand == "+" else "+")
                        )
                        row.update(
                            {
                                side + "_genomic_start_0based": start,
                                side + "_genomic_end_exclusive": end,
                                side + "_genomic_strand": strand,
                            }
                        )
                if specificity != "none":
                    try:
                        row.update(
                            check_pair(
                                row,
                                target,
                                bundle,
                                specificity,
                                max_product=max_product,
                                max_mismatches=max_mismatches,
                            )
                        )
                    except InputError as exc:
                        row.update(
                            specificity_status="error",
                            specificity_method=specificity,
                            specificity_error=str(exc),
                        )
                        checked_errors += 1
            rows.append(row)
    report = {
        "schema_version": 1,
        "created_at": datetime.now(timezone.utc).isoformat(),
        "software": {
            "python": platform.python_version(),
            "primer3-py": version("primer3-py"),
            "batch-primer-designer": "0.2.0",
        },
        "parameters": asdict(config),
        "provenance": provenance or {},
        "reference": bundle.manifest if bundle else None,
        "summary": {
            "input_records": len(targets),
            "successful_records": success,
            "failed_records": len(targets) - success,
            "candidate_pairs": sum(r["status"] == "ok" for r in rows),
            "specificity_errors": checked_errors,
        },
        "results": rows,
    }
    output.parent.mkdir(parents=True, exist_ok=True)
    with tempfile.TemporaryDirectory(prefix=".design-", dir=output.parent) as tmp:
        stage = Path(tmp)
        export_report(report, stage)
        manifest = {
            "schema_version": 1,
            "files": {p.name: sha256(p) for p in stage.iterdir()},
        }
        (stage / "manifest.json").write_text(json.dumps(manifest, indent=2))
        stage.rename(output)
    return report
