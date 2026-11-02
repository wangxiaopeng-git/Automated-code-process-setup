from __future__ import annotations

import argparse
import sys
from pathlib import Path
from typing import Optional, Sequence

from .core import DesignConfig, InputError, design_batch, read_records, write_results


def build_parser() -> argparse.ArgumentParser:
    parser = argparse.ArgumentParser(
        prog="batch-primer-design",
        description="Use Primer3 to design PCR primer pairs for FASTA/CSV/TSV templates.",
    )
    parser.add_argument("input", type=Path, help="input .fasta/.fa/.fna, .csv, or .tsv")
    parser.add_argument("-o", "--output", type=Path, default=Path("primer_results.csv"))
    parser.add_argument("--num-return", type=int, default=5)
    parser.add_argument("--primer-min-size", type=int, default=18)
    parser.add_argument("--primer-opt-size", type=int, default=20)
    parser.add_argument("--primer-max-size", type=int, default=25)
    parser.add_argument("--primer-min-tm", type=float, default=57.0)
    parser.add_argument("--primer-opt-tm", type=float, default=60.0)
    parser.add_argument("--primer-max-tm", type=float, default=63.0)
    parser.add_argument("--primer-min-gc", type=float, default=35.0)
    parser.add_argument("--primer-max-gc", type=float, default=65.0)
    parser.add_argument("--product-min-size", type=int, default=100)
    parser.add_argument("--product-max-size", type=int, default=300)
    parser.add_argument("--max-poly-x", type=int, default=4)
    parser.add_argument("--max-ns", type=int, default=0)
    parser.add_argument("--salt-monovalent-mm", type=float, default=50.0)
    parser.add_argument("--salt-divalent-mm", type=float, default=1.5)
    parser.add_argument("--dntp-mm", type=float, default=0.6)
    parser.add_argument("--primer-dna-nm", type=float, default=50.0)
    parser.add_argument(
        "--pick-anyway",
        action="store_true",
        help="allow Primer3 to return forced primers that violate constraints",
    )
    return parser


def main(argv: Optional[Sequence[str]] = None) -> int:
    args = build_parser().parse_args(argv)
    config = DesignConfig(
        num_return=args.num_return,
        primer_min_size=args.primer_min_size,
        primer_opt_size=args.primer_opt_size,
        primer_max_size=args.primer_max_size,
        primer_min_tm=args.primer_min_tm,
        primer_opt_tm=args.primer_opt_tm,
        primer_max_tm=args.primer_max_tm,
        primer_min_gc=args.primer_min_gc,
        primer_max_gc=args.primer_max_gc,
        product_min_size=args.product_min_size,
        product_max_size=args.product_max_size,
        max_poly_x=args.max_poly_x,
        max_ns=args.max_ns,
        salt_monovalent_mm=args.salt_monovalent_mm,
        salt_divalent_mm=args.salt_divalent_mm,
        dntp_mm=args.dntp_mm,
        primer_dna_nm=args.primer_dna_nm,
        pick_anyway=args.pick_anyway,
    )
    try:
        config.validate()
        if not args.input.is_file():
            raise InputError(f"input file does not exist: {args.input}")
        try:
            from primer3.bindings import design_primers
        except ImportError as exc:
            raise InputError(
                "primer3-py is not installed; run: python -m pip install -e ."
            ) from exc

        records = list(read_records(args.input))
        if not records:
            raise InputError("input contains no records")
        rows = list(design_batch(records, config, design_primers))
        write_results(rows, args.output)
    except (InputError, OSError) as exc:
        print(f"error: {exc}", file=sys.stderr)
        return 2

    ok_samples = len({row["record_id"] for row in rows if row["status"] == "ok"})
    failed_samples = len({row["record_id"] for row in rows if row["status"] != "ok"})
    print(
        f"done: {ok_samples} sample(s) with primers, "
        f"{failed_samples} failed; output: {args.output}"
    )
    return 0 if failed_samples == 0 else 1


if __name__ == "__main__":
    raise SystemExit(main())
