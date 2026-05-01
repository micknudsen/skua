from skua.evidence import (
    AlleleSupport,
    ReadAlleleCall,
    UnusableReason,
    aggregate_read_calls,
)


def test_aggregate_counts_strand_aware_alt_and_non_alt() -> None:
    calls = [
        ReadAlleleCall(support=AlleleSupport.ALT, is_reverse=False, observed_base="T", base_quality=35),
        ReadAlleleCall(support=AlleleSupport.ALT, is_reverse=False, observed_base="T", base_quality=30),
        ReadAlleleCall(support=AlleleSupport.ALT, is_reverse=True, observed_base="T", base_quality=32),
        ReadAlleleCall(support=AlleleSupport.NON_ALT, is_reverse=False, observed_base="A", base_quality=33),
        ReadAlleleCall(support=AlleleSupport.NON_ALT, is_reverse=True, observed_base="A", base_quality=31),
        ReadAlleleCall(support=AlleleSupport.NON_ALT, is_reverse=True, observed_base="A", base_quality=29),
    ]

    counts = aggregate_read_calls(calls)

    assert counts.alt_forward == 2
    assert counts.alt_reverse == 1
    assert counts.non_alt_forward == 1
    assert counts.non_alt_reverse == 2
    assert counts.usable == 6
    assert counts.unusable == 0


def test_aggregate_excludes_unusable_from_evidence_counts() -> None:
    calls = [
        ReadAlleleCall(support=AlleleSupport.ALT, is_reverse=False, observed_base="T", base_quality=35),
        ReadAlleleCall(
            support=AlleleSupport.UNUSABLE,
            is_reverse=False,
            reason=UnusableReason.LOW_BASEQ,
            observed_base="T",
            base_quality=10,
        ),
        ReadAlleleCall(
            support=AlleleSupport.UNUSABLE,
            is_reverse=True,
            reason=UnusableReason.LOW_MAPQ,
        ),
    ]

    counts = aggregate_read_calls(calls)

    assert counts.alt_forward == 1
    assert counts.alt_reverse == 0
    assert counts.non_alt_forward == 0
    assert counts.non_alt_reverse == 0
    assert counts.usable == 1
    assert counts.unusable == 2


def test_aggregate_tracks_unusable_reasons() -> None:
    calls = [
        ReadAlleleCall(
            support=AlleleSupport.UNUSABLE,
            is_reverse=False,
            reason=UnusableReason.LOW_BASEQ,
            observed_base="T",
            base_quality=10,
        ),
        ReadAlleleCall(
            support=AlleleSupport.UNUSABLE,
            is_reverse=True,
            reason=UnusableReason.LOW_BASEQ,
            observed_base="A",
            base_quality=12,
        ),
        ReadAlleleCall(
            support=AlleleSupport.UNUSABLE,
            is_reverse=False,
            reason=UnusableReason.INVALID_BASE,
            observed_base="N",
            base_quality=35,
        ),
    ]

    counts = aggregate_read_calls(calls)

    assert counts.unusable == 3
    assert counts.unusable_by_reason[UnusableReason.LOW_BASEQ] == 2
    assert counts.unusable_by_reason[UnusableReason.INVALID_BASE] == 1


def test_aggregate_empty_calls() -> None:
    counts = aggregate_read_calls([])

    assert counts.alt_forward == 0
    assert counts.alt_reverse == 0
    assert counts.non_alt_forward == 0
    assert counts.non_alt_reverse == 0
    assert counts.usable == 0
    assert counts.unusable == 0
    assert counts.unusable_by_reason == {}
