"""Read-level evidence classification primitives for variant verification."""

from dataclasses import dataclass
from enum import Enum
from typing import Any


class AlleleSupport(str, Enum):
    """High-level read support classification at a variant locus."""

    ALT = "alt"
    NON_ALT = "non_alt"
    UNUSABLE = "unusable"


class UnusableReason(str, Enum):
    """Reason for excluding a read from evidence counting."""

    LOW_MAPQ = "low_mapq"
    LOW_BASEQ = "low_baseq"
    NO_BASE_AT_SITE = "no_base_at_site"
    INVALID_BASE = "invalid_base"


@dataclass(frozen=True)
class ReadAlleleCall:
    """Result of classifying one read at one SNV locus."""

    support: AlleleSupport
    is_reverse: bool
    reason: UnusableReason | None = None
    observed_base: str | None = None
    base_quality: int | None = None


def classify_snv_read(
    read: Any,
    *,
    ref_pos0: int,
    ref_base: str,
    alt_base: str,
    min_baseq: int = 20,
    min_mapq: int = 20,
) -> ReadAlleleCall:
    """Classify one read as ALT, NON_ALT, or UNUSABLE for a single SNV.

    This is intentionally left unimplemented for TDD; tests define the behavior.
    """
    raise NotImplementedError("read-level SNV classification is not implemented yet")
