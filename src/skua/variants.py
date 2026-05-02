"""Variant parsing and normalization helpers."""

from dataclasses import dataclass


@dataclass(frozen=True)
class Variant:
    """Minimal SNV variant model using 0-based reference position."""

    contig: str
    ref_pos0: int
    ref: str
    alt: str

    @classmethod
    def from_vcf_fields(cls, *, contig: str, pos1: int, ref: str, alt: str) -> "Variant":
        """Build a Variant from basic VCF fields.

        This is intentionally left unimplemented for TDD.
        """
        raise NotImplementedError("Variant.from_vcf_fields is not implemented yet")


def parse_vcf_snv_line(line: str) -> Variant | None:
    """Parse one VCF line and return a SNV Variant when applicable.

    This is intentionally left unimplemented for TDD.
    """
    raise NotImplementedError("parse_vcf_snv_line is not implemented yet")
