"""Local-first multi-species CLI. Run `primer-workflow --help`."""

from __future__ import annotations
import argparse
import json
import sys
from contextlib import ExitStack
from dataclasses import fields
from pathlib import Path
from .core import DesignConfig, InputError, RecordIssue
from .inputs import file_kind, load_templates, open_text
from .reference import ReferenceBundle, Target, import_reference, sha256
from .workflow import run_design


def parser():
    p = argparse.ArgumentParser(
        prog="primer-workflow",
        description="Multi-format genomic PCR primer design and reference management",
    )
    commands = p.add_subparsers(dest="command", required=True)
    design = commands.add_parser(
        "design", help="design candidates and export a complete run directory"
    )
    design.add_argument("input", type=Path)
    design.add_argument("-o", "--output", type=Path, required=True)
    design.add_argument("--reference", type=Path)
    design.add_argument("--flank", type=int, default=300)
    design.add_argument(
        "--config",
        type=Path,
        help="JSON object with DesignConfig fields, concentrations in mM/nM",
    )
    design.add_argument(
        "--specificity", choices=["none", "exact", "blast"], default="none"
    )
    design.add_argument("--specificity-max-product", type=int, default=2000)
    design.add_argument("--max-mismatches", type=int, default=2)
    design.add_argument("--num-return", type=int)
    design.add_argument("--product-min-size", type=int)
    design.add_argument("--product-max-size", type=int)
    reference = commands.add_parser("reference")
    actions = reference.add_subparsers(dest="action", required=True)
    imp = actions.add_parser(
        "import", help="create a validated reference bundle from local files"
    )
    imp.add_argument("--genome", type=Path, required=True)
    imp.add_argument("--annotation", type=Path)
    imp.add_argument("--taxid", required=True)
    imp.add_argument("--assembly", required=True)
    imp.add_argument("--annotation-release", default="unspecified")
    imp.add_argument("-o", "--output", type=Path, required=True)
    verify = actions.add_parser("verify")
    verify.add_argument("bundle", type=Path)
    search = actions.add_parser(
        "search",
        help="search NCBI Datasets; choose a versioned accession from its reports",
    )
    search.add_argument("--taxon", required=True)
    search.add_argument("--limit", type=int, default=20)
    download = actions.add_parser(
        "download", help="download versioned NCBI assembly and available annotation"
    )
    download.add_argument("--accession", required=True)
    download.add_argument("--taxid", required=True)
    download.add_argument("-o", "--output", type=Path, required=True)
    download.add_argument("--max-download-gb", type=float, default=10)
    return p


def main(argv=None):
    args = parser().parse_args(argv)
    try:
        if args.command == "reference":
            if args.action == "import":
                result = import_reference(
                    args.genome,
                    args.annotation,
                    args.output,
                    taxid=args.taxid,
                    assembly=args.assembly,
                    annotation_release=args.annotation_release,
                )
            elif args.action == "verify":
                with ReferenceBundle(args.bundle) as bundle:
                    bundle.verify()
                    result = {
                        "verified": str(args.bundle),
                        "assembly": bundle.manifest["assembly"],
                    }
            elif args.action == "search":
                from .ncbi import search_ncbi

                result = search_ncbi(args.taxon, limit=args.limit)
            else:
                from .ncbi import download_ncbi

                if not 0 < args.max_download_gb <= 100:
                    raise InputError("max-download-gb must be >0 and <=100")
                result = download_ncbi(
                    args.accession,
                    args.output,
                    taxid=args.taxid,
                    max_bytes=int(args.max_download_gb * 1e9),
                )
            print(json.dumps(result, indent=2, ensure_ascii=False))
            return 0
        settings = json.loads(args.config.read_text()) if args.config else {}
        if not isinstance(settings, dict):
            raise InputError("config must be a JSON object")
        if set(settings) - {f.name for f in fields(DesignConfig)}:
            raise InputError("unknown DesignConfig fields")
        for key in ("num_return", "product_min_size", "product_max_size"):
            if getattr(args, key) is not None:
                settings[key] = getattr(args, key)
        config = DesignConfig(**settings)
        config.validate()
        with ExitStack() as stack:
            bundle = (
                stack.enter_context(ReferenceBundle(args.reference))
                if args.reference
                else None
            )
            if bundle:
                bundle.verify()
            kind = file_kind(args.input)
            if kind in {".bed", ".txt"}:
                if bundle is None:
                    raise InputError(
                        "BED and gene/transcript ID lists require --reference"
                    )
                if kind == ".bed":
                    targets = bundle.read_bed(args.input, args.flank)
                else:
                    targets = []
                    with open_text(args.input) as h:
                        for line in h:
                            identifier = line.strip()
                            if not identifier or identifier.startswith("#"):
                                continue
                            try:
                                targets.append(
                                    bundle.resolve_gene(identifier, args.flank)
                                )
                            except InputError as exc:
                                targets.append(RecordIssue(identifier, str(exc)))
            else:
                targets = [
                    r if isinstance(r, RecordIssue) else Target(r)
                    for r in load_templates(args.input)
                ]
            provenance = {
                "input_name": args.input.name,
                "input_sha256": sha256(args.input),
                "coordinate_convention": "template 1-based closed; reference/BED 0-based half-open",
                "flank_bp": args.flank,
                "mode": "genomic PCR; gene/transcript IDs target their full genomic span",
            }
            report = run_design(
                targets,
                config,
                args.output,
                bundle=bundle,
                specificity=args.specificity,
                max_product=args.specificity_max_product,
                max_mismatches=args.max_mismatches,
                provenance=provenance,
            )
        print(json.dumps(report["summary"], ensure_ascii=False))
        print(f"output: {args.output}")
        return (
            1
            if report["summary"]["failed_records"]
            or report["summary"]["specificity_errors"]
            else 0
        )
    except (InputError, OSError, ValueError, TypeError, KeyError) as exc:
        print(f"error: {exc}", file=sys.stderr)
        return 2


if __name__ == "__main__":
    raise SystemExit(main())
