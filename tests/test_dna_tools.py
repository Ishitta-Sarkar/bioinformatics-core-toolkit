"""Unit tests for DNATools, including the checked examples from 24 Nov 2026."""

import warnings
from pathlib import Path

import pytest

from src.dna_tools import DNATools, LengthMismatchWarning


def test_gc_and_translation():
    s = DNATools("ATGTTTAAAGGGCCCTAG")
    assert s.gc_content() == 44.44
    assert s.translate_to_protein() == "MFKGP*"
    assert s.reverse_complement() == "CTAGGGCCCTTTAAACAT"


def test_mutations():
    muts = DNATools.detect_mutations("ACGT", "AGGT")
    assert muts == [{"position": 2, "reference": "C", "sample": "G"}]


def test_cleaning_and_validation():
    dna = DNATools(" atg\ntt n ")
    assert dna.sequence == "ATGTTN"
    assert dna.length() == 6
    assert dna.is_valid_dna() is False
    assert DNATools("atgc").is_valid_dna() is True
    assert DNATools("").is_valid_dna() is True
    assert DNATools("   \n").sequence == ""
    assert DNATools("u").is_valid_dna() is False


def test_sequence_must_be_a_string():
    with pytest.raises(TypeError):
        DNATools(None)


def test_gc_content_edges_and_rounding():
    assert DNATools("").gc_content() == 0
    assert DNATools("G").gc_content() == 100.0
    assert DNATools("AT").gc_content() == 0.0
    assert DNATools("GC").gc_content() == 100.0
    # 1/6 * 100 = 16.666... -> 16.67; 1/3 * 100 = 33.333... -> 33.33.
    assert DNATools("GAAAAA").gc_content() == 16.67
    assert DNATools("GAA").gc_content() == 33.33
    # N is not G or C, but it still counts in the length.
    assert DNATools("ATGN").gc_content() == 25.0


def test_invalid_dna_returns_a_message_instead_of_a_sequence():
    dna = DNATools("ATGN")
    assert dna.transcribe() == "Invalid DNA sequence"
    assert dna.reverse_complement() == "Invalid DNA sequence"
    assert dna.translate_to_protein() == "Invalid DNA sequence"


def test_transcription_and_reverse_complement_are_invertible_for_valid_dna():
    dna = DNATools("atg ttt aaa ggg ccc tag")
    assert dna.transcribe() == "AUGUUUAAAGGGCCCUAG"
    assert dna.reverse_complement() == "CTAGGGCCCTTTAAACAT"
    assert DNATools(dna.reverse_complement()).reverse_complement() == dna.sequence
    assert DNATools("").transcribe() == ""
    assert DNATools("").reverse_complement() == ""
    assert DNATools("").translate_to_protein() == ""


def test_translation_drops_a_trailing_incomplete_codon_and_keeps_stops():
    assert DNATools("ATGTT").translate_to_protein() == "M"
    assert DNATools("TAATAGTGA").translate_to_protein() == "***"
    assert DNATools("ATGGTTTAA").translate_to_protein() == "MV*"
    assert len(DNATools.codon_table) == 64
    assert DNATools.codon_table["ATG"] == "M"
    assert DNATools.codon_table["TGG"] == "W"
    assert DNATools.codon_table["TAA"] == "*"
    assert DNATools.codon_table["TAG"] == "*"
    assert DNATools.codon_table["TGA"] == "*"


def test_detect_mutations_reports_every_substitution_in_order():
    mutations = DNATools.detect_mutations("ACGT", "TGCA")
    assert mutations == [
        {"position": 1, "reference": "A", "sample": "T"},
        {"position": 2, "reference": "C", "sample": "G"},
        {"position": 3, "reference": "G", "sample": "C"},
        {"position": 4, "reference": "T", "sample": "A"},
    ]
    assert DNATools.detect_mutations("ACGT", "ACGT") == []
    assert DNATools.detect_mutations("acgt", "AgGt") == [
        {"position": 2, "reference": "C", "sample": "G"},
    ]


def test_equal_lengths_do_not_warn():
    with warnings.catch_warnings():
        warnings.simplefilter("error", LengthMismatchWarning)
        assert DNATools.detect_mutations("", "") == []
        DNATools.detect_mutations("ACGT", "AGGT")


def test_length_mismatch_warns_and_keeps_only_the_shared_prefix():
    with pytest.warns(
        LengthMismatchWarning, match="does not call insertions or deletions"
    ):
        mutations = DNATools.detect_mutations("ACGT", "AGGAAA")

    assert mutations == [
        {"position": 2, "reference": "C", "sample": "G"},
        {"position": 4, "reference": "T", "sample": "A"},
    ]

    with pytest.warns(LengthMismatchWarning, match=r"Reference length \(4\)"):
        assert DNATools.detect_mutations("ACGT", "") == []

    with pytest.warns(LengthMismatchWarning, match=r"sample length \(6\)"):
        assert DNATools.detect_mutations("ACGT", "ACGTAA") == []


def test_whitespace_is_not_stripped_inside_detect_mutations():
    with pytest.warns(LengthMismatchWarning):
        mutations = DNATools.detect_mutations("AC GT", "ACGT")
    assert mutations == [
        {"position": 3, "reference": " ", "sample": "G"},
        {"position": 4, "reference": "G", "sample": "T"},
    ]


def test_equal_length_can_hide_a_compensating_indel():
    """Deleting one base and inserting another keeps the length the same.

    Reference ``ATGCAT`` with the T at position 2 removed and a T added at
    the end is ``AGCATT``. The lengths match, so no warning is emitted, and
    the ungapped comparison reports substitutions. Those positions are not
    an alignment.
    """
    reference = "ATGCAT"
    sample = reference[0] + reference[2:] + "T"
    assert sample == "AGCATT"
    with warnings.catch_warnings():
        warnings.simplefilter("error", LengthMismatchWarning)
        mutations = DNATools.detect_mutations(reference, sample)
    assert [item["position"] for item in mutations] == [2, 3, 4, 5]


def test_detect_mutations_rejects_non_strings():
    with pytest.raises(TypeError):
        DNATools.detect_mutations(None, "ACGT")
    with pytest.raises(TypeError):
        DNATools.detect_mutations("ACGT", b"ACGT")


def test_parse_fasta_multiline_duplicate_headers_and_leading_sequence(tmp_path: Path):
    path = tmp_path / "records.fasta"
    path.write_text(
        "ATGC\n>alpha example\natg\nc\n;not-a-comment\n>alpha example\nTT\n>beta\nGG\n",
        encoding="utf-8",
    )
    parsed = DNATools.parse_fasta(path)
    assert parsed == {
        "alpha example": "TT",
        "beta": "GG",
    }

    comment_path = tmp_path / "comment.fasta"
    comment_path.write_text(">only\natgc\n;not-a-comment\n", encoding="utf-8")
    assert DNATools.parse_fasta(comment_path) == {"only": "atgc;not-a-comment"}


def test_parse_fasta_missing_file():
    with pytest.raises(FileNotFoundError):
        DNATools.parse_fasta("no-such-sequence.fasta")
