"""Local-first multi-species CLI. Run `primer-workflow --help`."""

from __future__ import annotations
import argparse
import json
import shutil
import sys
from contextlib import ExitStack
from dataclasses import fields
from pathlib import Path
from .core import DesignConfig, InputError, RecordIssue
from .inputs import file_kind, load_templates, open_text
from .reference import ReferenceBundle, Target, import_reference, sha256
from .workflow import run_design
from . import __version__


def parser():
    p = argparse.ArgumentParser(
        prog="primer-workflow",
        description="Multi-format genomic PCR primer design and reference management",
    )
    p.add_argument("--version", action="version", version="%(prog)s " + __version__)
    commands = p.add_subparsers(dest="command", required=True)
    design = commands.add_parser(
        "design", help="design candidates and export a complete run directory"
    )
    design.add_argument("input", type=Path)
    design.add_argument(
        "-o",
        "--output",
        type=Path,
        help="Output directory; default: results/<input project folder>",
    )
    source = design.add_mutually_exclusive_group()
    source.add_argument(
        "--reference",
        type=Path,
        help="Genome FASTA(.gz) file, or an existing reference bundle",
    )
    source.add_argument(
        "--assembly",
        help="Experimental online NCBI mode: versioned GCF_/GCA_ assembly, no full genome download",
    )
    design.add_argument(
        "--ncbi-email", help="Contact email supplied to NCBI requests (optional)"
    )
    design.add_argument(
        "--max-region-bp",
        type=int,
        default=200000,
        help="Maximum online target plus flanks, 1..1000000 bp",
    )
    design.add_argument(
        "--remote-hit-limit",
        type=int,
        default=100,
        help="NCBI BLAST hits per primer, 1..500; saturation is incomplete",
    )
    design.add_argument(
        "--remote-max-queries",
        type=int,
        default=40,
        help="Maximum online primer queries per run, 1..100",
    )
    design.add_argument(
        "--remote-timeout",
        type=int,
        default=900,
        help="Online BLAST polling timeout, 60..3600 seconds",
    )
    annotation = design.add_mutually_exclusive_group()
    annotation.add_argument(
        "--gtf", type=Path, help="Matching GTF annotation (optional for sequence input)"
    )
    annotation.add_argument("--gff", type=Path, help="Matching GFF/GFF3 annotation")
    design.add_argument(
        "--cache-dir",
        type=Path,
        help="Reference cache; default: .primer-cache inside the input project folder",
    )
    design.add_argument("--flank", type=int, default=300)
    design.add_argument(
        "--config",
        type=Path,
        help="JSON object with DesignConfig fields, concentrations in mM/nM",
    )
    design.add_argument(
        "--specificity",
        choices=["none", "exact", "blast", "remote"],
        default=None,
        help="Default: remote with --assembly; blast with --reference; otherwise none",
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
    inspect = actions.add_parser(
        "inspect",
        help="Retrieve assembly and sequence metadata only, without whole genome download",
    )
    inspect.add_argument("--assembly", required=True)
    inspect.add_argument("--cache-dir", type=Path, default=Path("data/.ncbi-cache"))
    inspect.add_argument("--ncbi-email")
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
            elif args.action == "inspect":
                from .online_reference import OnlineReference

                with OnlineReference(
                    args.assembly, args.cache_dir, email=args.ncbi_email
                ) as online:
                    result = online.metadata()
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
        if (args.gtf or args.gff) and not args.reference:
            raise InputError(
                "--gtf/--gff requires local --reference genome FASTA; online --assembly retrieves region features automatically"
            )
        args.specificity = args.specificity or (
            "remote" if args.assembly else "blast" if args.reference else "none"
        )
        if args.specificity != "none" and not (args.reference or args.assembly):
            raise InputError("Specificity checking requires --reference genome FASTA")
        if args.assembly and args.specificity not in {"none", "remote"}:
            raise InputError(
                "--assembly supports --specificity remote/none; exact/blast require a local --reference"
            )
        if args.specificity == "remote" and not args.assembly:
            raise InputError("--specificity remote requires --assembly")
        if args.specificity == "blast" and not shutil.which("blastn"):
            raise InputError(
                "BLAST+ is required for reference specificity checking. Install ncbi-blast+ (Linux), brew install blast (macOS), or the official Windows BLAST+ package; run blastn -version to verify. Use --specificity exact only for references <=5 Mb; --specificity none explicitly skips checking."
            )
        args.output = args.output or Path("results") / (
            args.input.resolve().parent.name or args.input.stem
        )
        if not args.input.is_file():
            raise InputError(f"input file does not exist: {args.input}")
        if args.output.exists():
            raise InputError("output exists; choose a new run directory")
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
            bundle = None
            online = None
            checker = None
            if args.assembly:
                from .online_reference import OnlineReference

                online = stack.enter_context(
                    OnlineReference(
                        args.assembly,
                        args.cache_dir or args.input.resolve().parent / ".ncbi-cache",
                        email=args.ncbi_email,
                        max_region_bp=args.max_region_bp,
                    )
                )
                if args.specificity == "remote":
                    from .remote_blast import OnlineChecker

                    checker = OnlineChecker(
                        online,
                        args.cache_dir or args.input.resolve().parent / ".ncbi-cache",
                        email=args.ncbi_email,
                        hit_limit=args.remote_hit_limit,
                        poll_timeout=args.remote_timeout,
                        max_queries=args.remote_max_queries,
                    )
                    stack.callback(checker.close)
            if args.reference:
                from .auto_reference import prepare_reference

                reference_path = prepare_reference(
                    args.reference,
                    args.gtf or args.gff,
                    cache_dir=args.cache_dir
                    or args.input.resolve().parent / ".primer-cache",
                )
                bundle = stack.enter_context(ReferenceBundle(reference_path))
            if bundle:
                bundle.verify()
            kind = file_kind(args.input)
            if kind in {".bed", ".txt"}:
                if bundle is None and online is None:
                    raise InputError(
                        "BED and gene/transcript ID lists require --reference or --assembly"
                    )
                if kind == ".bed":
                    targets = (online or bundle).read_bed(args.input, args.flank)
                else:
                    targets = []
                    with open_text(args.input) as h:
                        for line in h:
                            identifier = line.strip()
                            if not identifier or identifier.startswith("#"):
                                continue
                            try:
                                targets.append(
                                    online.gene(identifier, args.flank)
                                    if online
                                    else bundle.resolve_gene(identifier, args.flank)
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
                "online_assembly": args.assembly,
                "reference_input": str(args.reference) if args.reference else None,
                "annotation_input": (
                    str(args.gtf or args.gff) if (args.gtf or args.gff) else None
                ),
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
                remote_checker=checker,
                reference_metadata=online.metadata() if online else None,
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
