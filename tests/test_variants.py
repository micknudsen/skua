import pytest

from skua.variants import Variant, parse_vcf_snv_line


def test_variant_from_vcf_fields_converts_pos1_to_pos0() -> None:
    variant = Variant.from_vcf_fields(contig="chr7", pos1=106, ref="A", alt="T")

    assert variant.contig == "chr7"
    assert variant.ref_pos0 == 105
    assert variant.ref == "A"
    assert variant.alt == "T"


def test_variant_from_vcf_fields_rejects_non_snv_ref() -> None:
    with pytest.raises(ValueError, match="SNV"):
        Variant.from_vcf_fields(contig="chr1", pos1=10, ref="AT", alt="A")


def test_variant_from_vcf_fields_rejects_non_snv_alt() -> None:
    with pytest.raises(ValueError, match="SNV"):
        Variant.from_vcf_fields(contig="chr1", pos1=10, ref="A", alt="AT")


def test_parse_vcf_snv_line_parses_data_line() -> None:
    line = "chr1\t106\t.\tA\tT\t.\tPASS\t." 

    variant = parse_vcf_snv_line(line)

    assert variant is not None
    assert variant.contig == "chr1"
    assert variant.ref_pos0 == 105
    assert variant.ref == "A"
    assert variant.alt == "T"


def test_parse_vcf_snv_line_skips_header_lines() -> None:
    assert parse_vcf_snv_line("##fileformat=VCFv4.2") is None
    assert parse_vcf_snv_line("#CHROM\tPOS\tID\tREF\tALT\tQUAL\tFILTER\tINFO") is None


def test_parse_vcf_snv_line_skips_non_snv_records() -> None:
    assert parse_vcf_snv_line("chr1\t106\t.\tA\tAT\t.\tPASS\t.") is None
    assert parse_vcf_snv_line("chr1\t106\t.\tAT\tA\t.\tPASS\t.") is None


def test_parse_vcf_snv_line_skips_multiallelic_records() -> None:
    assert parse_vcf_snv_line("chr1\t106\t.\tA\tT,C\t.\tPASS\t.") is None
