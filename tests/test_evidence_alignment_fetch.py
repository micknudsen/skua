from dataclasses import dataclass

from skua.evidence import UnusableReason, collect_snv_evidence_from_alignment


@dataclass
class FakeRead:
    mapping_quality: int
    is_reverse: bool
    query_sequence: str
    query_qualities: list[int]
    aligned_pairs: list[tuple[int | None, int | None]]


class FakeAlignmentFile:
    def __init__(self, reads: list[FakeRead]) -> None:
        self._reads = reads
        self.fetch_calls: list[tuple[str, int, int]] = []

    def fetch(self, contig: str, start: int, stop: int):
        self.fetch_calls.append((contig, start, stop))
        return iter(self._reads)



def build_linear_pairs(read_len: int, ref_start: int) -> list[tuple[int, int]]:
    return [(qpos, ref_start + qpos) for qpos in range(read_len)]



def test_collect_snv_evidence_from_alignment_fetches_one_locus_window() -> None:
    reads = [
        FakeRead(
            mapping_quality=60,
            is_reverse=False,
            query_sequence="AAAAATAAAA",
            query_qualities=[35] * 10,
            aligned_pairs=build_linear_pairs(10, 100),
        )
    ]
    alignment_file = FakeAlignmentFile(reads)

    counts = collect_snv_evidence_from_alignment(
        alignment_file,
        contig="chr7",
        ref_pos0=105,
        ref_base="A",
        alt_base="T",
    )

    assert alignment_file.fetch_calls == [("chr7", 105, 106)]
    assert counts.alt_forward == 1
    assert counts.alt_reverse == 0
    assert counts.non_alt_forward == 0
    assert counts.non_alt_reverse == 0



def test_collect_snv_evidence_from_alignment_propagates_mixed_counts() -> None:
    reads = [
        FakeRead(
            mapping_quality=60,
            is_reverse=False,
            query_sequence="AAAAATAAAA",
            query_qualities=[35] * 10,
            aligned_pairs=build_linear_pairs(10, 100),
        ),
        FakeRead(
            mapping_quality=60,
            is_reverse=True,
            query_sequence="AAAAAAAAAA",
            query_qualities=[35] * 10,
            aligned_pairs=build_linear_pairs(10, 100),
        ),
        FakeRead(
            mapping_quality=5,
            is_reverse=False,
            query_sequence="AAAAATAAAA",
            query_qualities=[35] * 10,
            aligned_pairs=build_linear_pairs(10, 100),
        ),
    ]
    alignment_file = FakeAlignmentFile(reads)

    counts = collect_snv_evidence_from_alignment(
        alignment_file,
        contig="chr1",
        ref_pos0=105,
        ref_base="A",
        alt_base="T",
        min_baseq=20,
        min_mapq=20,
    )

    assert counts.alt_forward == 1
    assert counts.alt_reverse == 0
    assert counts.non_alt_forward == 0
    assert counts.non_alt_reverse == 1
    assert counts.usable == 2
    assert counts.unusable == 1
    assert counts.unusable_by_reason[UnusableReason.LOW_MAPQ] == 1
