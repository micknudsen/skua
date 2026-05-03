"""Command-line interface for skua."""

import argparse
from pathlib import Path

import pysam

from .core import verify_snv_vcf_to_json


def build_parser() -> argparse.ArgumentParser:
    """Build and return the top-level argument parser."""
    parser = argparse.ArgumentParser(prog="skua")
    subparsers = parser.add_subparsers(dest="command", required=True)

    verify_parser = subparsers.add_parser("verify", help="Verify SNV evidence from VCF and BAM/CRAM")
    verify_parser.add_argument("--vcf", required=True, help="Input VCF path")
    verify_parser.add_argument("--alignment", required=True, help="Input BAM/CRAM path")
    verify_parser.add_argument("--reference", help="Reference FASTA path (required for CRAM)")
    verify_parser.add_argument("--output", help="Optional output JSON path")
    verify_parser.add_argument("--min-baseq", type=int, default=20, help="Minimum base quality")
    verify_parser.add_argument("--min-mapq", type=int, default=20, help="Minimum mapping quality")

    return parser


def main(argv: list[str] | None = None) -> int:
    """Run the skua CLI."""
    parser = build_parser()
    args = parser.parse_args(argv)

    if args.command == "verify":
        alignment_path = Path(args.alignment)
        if alignment_path.suffix.lower() == ".cram" and args.reference is None:
            parser.error("--reference is required for CRAM input")

        alignment_kwargs: dict[str, str] = {}
        if args.reference is not None:
            alignment_kwargs["reference_filename"] = args.reference

        with pysam.AlignmentFile(args.alignment, "rb", **alignment_kwargs) as alignment_file:
            payload = verify_snv_vcf_to_json(
                alignment_file,
                Path(args.vcf),
                output_path=args.output,
                min_baseq=args.min_baseq,
                min_mapq=args.min_mapq,
            )

        if args.output is None:
            print(payload)
        return 0

    parser.error(f"Unknown command: {args.command}")
    return 2
