"""Human reports and machine-readable outputs with explicit coordinate systems."""

import csv
import html
import json
from pathlib import Path
from openpyxl import Workbook


def spreadsheet_safe(value):
    if isinstance(value, str) and value.startswith(("=", "+", "-", "@", "\t", "\r")):
        return "'" + value
    return value


def export_report(report, output):
    output = Path(output)
    rows = report["results"]
    scalar = [
        {
            k: (json.dumps(v, ensure_ascii=False) if isinstance(v, (dict, list)) else v)
            for k, v in r.items()
        }
        for r in rows
    ]
    fields = list(dict.fromkeys(k for r in scalar for k in r))
    for filename, delimiter in [("results.csv", ","), ("results.tsv", "\t")]:
        with (output / filename).open("w", encoding="utf-8-sig", newline="") as h:
            w = csv.DictWriter(h, fieldnames=fields, delimiter=delimiter)
            w.writeheader()
            for r in scalar:
                w.writerow({k: spreadsheet_safe(v) for k, v in r.items()})
    book = Workbook()
    sheet = book.active
    sheet.title = "Candidates"
    sheet.append(fields)
    for r in scalar:
        sheet.append([spreadsheet_safe(r.get(k, "")) for k in fields])
    sheet.freeze_panes = "A2"
    sheet.auto_filter.ref = sheet.dimensions
    book.save(output / "results.xlsx")
    (output / "results.json").write_text(
        json.dumps(report, ensure_ascii=False, indent=2, allow_nan=False),
        encoding="utf-8",
    )
    if (report.get("reference") or {}).get("provider") == "NCBI":
        (output / "online_reference.json").write_text(
            json.dumps(report["reference"], ensure_ascii=False, indent=2),
            encoding="utf-8",
        )
        regions = report["reference"].get("regions", [])
        if regions:
            (output / "online_templates.fasta").write_text(
                "".join(
                    f">online_target{i} sample={str(r['sample_id']).replace(chr(10),' ')} {r['contig']}:{r['start_0based']}-{r['end_exclusive']} strand={r['strand']}\n{r['template_sequence_5to3']}\n"
                    for i, r in enumerate(regions, 1)
                ),
                encoding="utf-8",
            )
            with (output / "region_features.tsv").open(
                "w", encoding="utf-8-sig", newline=""
            ) as handle:
                writer = csv.writer(handle, delimiter="\t")
                writer.writerow(
                    [
                        "sample_id",
                        "contig",
                        "feature_type",
                        "part_index",
                        "start_1based",
                        "end_1based",
                        "strand",
                        "qualifiers_json",
                        "location_in_fetched_record",
                    ]
                )
                for r in regions:
                    for f in r["features"]:
                        for i, part in enumerate(f["parts"], 1):
                            writer.writerow(
                                [
                                    spreadsheet_safe(str(r["sample_id"])),
                                    r["contig"],
                                    f["type"],
                                    i,
                                    part["start_0based"] + 1,
                                    part["end_exclusive"],
                                    part["strand"],
                                    json.dumps(f["qualifiers"], ensure_ascii=False),
                                    f["location_in_fetched_record"],
                                ]
                            )
    primers = []
    amplicons = []
    bed = []
    for row in rows:
        if row["status"] != "ok":
            continue
        key = f"r{row['record_id']}_pair{row['pair_rank']}"
        for side in ("left", "right"):
            primers.append(f">{key}_{side}\n{row[side+'_sequence_5to3']}\n")
            if row.get("contig"):
                bed.append(
                    "\t".join(
                        map(
                            str,
                            [
                                row["contig"],
                                row[side + "_genomic_start_0based"],
                                row[side + "_genomic_end_exclusive"],
                                key + "_" + side,
                                0,
                                row[side + "_genomic_strand"],
                            ],
                        )
                    )
                    + "\n"
                )
        amplicons.append(
            f">{key} sample={str(row['sample_id']).replace(chr(10),' ')}\n{row['amplicon_sequence_5to3']}\n"
        )
    (output / "primers.fasta").write_text("".join(primers))
    (output / "amplicons.fasta").write_text("".join(amplicons))
    (output / "primers.bed").write_text("".join(bed))
    headers = [
        "record_id",
        "sample_id",
        "pair_rank",
        "status",
        "specificity_status",
        "screening_priority",
        "left_sequence_5to3",
        "right_sequence_5to3",
        "product_size_bp",
        "message",
    ]
    table = (
        "<tr>" + "".join("<th>" + html.escape(k) + "</th>" for k in headers) + "</tr>"
    )
    for row in rows[:1000]:
        table += (
            "<tr>"
            + "".join(
                "<td>" + html.escape(str(row.get(k, ""))) + "</td>" for k in headers
            )
            + "</tr>"
        )
    (output / "report.html").write_text(
        '<!doctype html><meta charset="utf-8"><title>Primer report</title><style>body{font:14px sans-serif;margin:2em}td,th{border:1px solid #ccc;padding:6px}table{border-collapse:collapse}</style><h1>引物设计报告</h1><p>候选生成不等于实验验证。not_run 表示未检查特异性；单一预测产物仅适用于记录的检索条件。最多显示 1000 行，完整数据见 JSON/表格。</p><pre>'
        + html.escape(json.dumps(report["summary"], ensure_ascii=False, indent=2))
        + "</pre><table>"
        + table
        + "</table>",
        encoding="utf-8",
    )
