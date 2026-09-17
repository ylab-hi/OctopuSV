from octopusv.utils.svcf_coordinate_parser import parse_svcf_co


def test_parse_simple_co():
    assert parse_svcf_co("chr1_100-chr2_200") == ("chr1", 100, "chr2", 200)


def test_parse_contigs_with_underscores():
    assert parse_svcf_co(
        "chrY_KI270740v1_random_100-NC_007605_200"
    ) == ("chrY_KI270740v1_random", 100, "NC_007605", 200)


def test_parse_contigs_with_hyphens_when_unambiguous():
    assert parse_svcf_co(
        "chr1-alt_100-chr2-alt_200"
    ) == ("chr1-alt", 100, "chr2-alt", 200)


def test_ambiguous_co_is_not_guessed():
    # Two hyphens can each produce a syntactically valid chrom_pos split.
    assert parse_svcf_co(
        "chrY_KI270740v1_random_123-NC_007605-alt_456"
    ) is None


def test_missing_or_malformed_co_returns_none():
    assert parse_svcf_co(".") is None
    assert parse_svcf_co(None) is None
    assert parse_svcf_co("chr1-x-chr2_200") is None


def test_zero_coordinate_is_structurally_parseable():
    assert parse_svcf_co("chr1_0-chr2_200") == ("chr1", 0, "chr2", 200)
