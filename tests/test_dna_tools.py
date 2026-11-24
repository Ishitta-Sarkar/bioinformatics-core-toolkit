"""Unit tests for DNATools, including the checked example values."""

from __future__ import annotations

import warnings
from pathlib import Path

import pytest

from src.dna_tools import DNATools, LengthDifferenceWarning

# NP_000509.1, the /translation qualifier of NM_000518.5 CDS 51..494.
# The terminal TAA is not included in that qualifier.
NP_000509_1 = (
    "MVHLTPEEKSAVTALWGKVNVDEVGGEALGRLLVVYPWTQRFFESFGDLSTPDAVMGNPKVKAHGKK"
    "VLGAFSDGLAHLDNLKGTFATLSELHCDKLHVDPENFRLLGNVLVCVLAHHFGKEFTPPVQAAYQKVV"
    "AGVANALAHKYH"
)
EXAMPLE = Path(__file__).resolve().parents[1] / "examples" / "sample_sequence.fasta"


def test_gc_and_translation():
    s = DNATools("ATGTTTAAAGGGCCCTAG")
    assert s.gc_content() == 44.44
    assert s.translate_to_protein() == "MFKGP*"
    assert s.reverse_complement() == "CTAGGGCCCTTTAAACAT"


def test_mutations():
    muts = DNATools.detect_mutations("ACGT", "AGGT")
    assert muts == [{"position": 2, "reference": "C", "sample": "G"}]


def test_gc_content_edges():
    assert DNATools("").gc_content() == 0
    assert DNATools("   \n").gc_content() == 0
    assert DNATools("GGGG").gc_content() == 100
    assert DNATools("ATAT").gc_content() == 0
    assert DNATools("GCGC").gc_content() == 100
    assert DNATools("ACGT").gc_content() == 50
    # N counts in the length and is not GC: 1 GC base of 2.
    assert DNATools("GN").gc_content() == 50


def test_case_and_whitespace_are_removed():
    dna = DNATools(" atg ttt aaa ggg ccc tag\n")
    assert dna.sequence == "ATGTTTAAAGGGCCCTAG"
    assert dna.is_valid_dna()
    assert dna.length() == 18
    assert dna.gc_content() == 44.44


def test_is_valid_dna():
    assert DNATools("acgt").is_valid_dna()
    assert DNATools("").is_valid_dna()
    assert not DNATools("ACGN").is_valid_dna()
    assert not DNATools("ACGU").is_valid_dna()


def test_invalid_dna_returns_message_and_does_not_raise():
    dna = DNATools("ATGNA")
    assert dna.transcribe() == "Invalid DNA sequence"
    assert dna.reverse_complement() == "Invalid DNA sequence"
    assert dna.translate_to_protein() == "Invalid DNA sequence"
    assert dna.length() == 5


def test_transcribe_and_reverse_complement_round_trip():
    dna = DNATools("ATGTTTAAAGGGCCCTAG")
    assert dna.transcribe() == "AUGUUUAAAGGGCCCUAG"
    assert DNATools(dna.reverse_complement()).reverse_complement() == dna.sequence


def test_translation_drops_a_trailing_incomplete_codon():
    assert DNATools("ATGTT").translate_to_protein() == "M"
    assert DNATools("ATG").translate_to_protein() == "M"
    assert DNATools("AT").translate_to_protein() == ""
    assert DNATools("TAATAGTGAATG").translate_to_protein() == "***M"


def test_codon_table_has_the_standard_sixty_four_codons():
    assert len(DNATools.codon_table) == 64
    assert set(DNATools.codon_table) == {
        a + b + c for a in "ACGT" for b in "ACGT" for c in "ACGT"
    }
    assert DNATools.codon_table["ATG"] == "M"
    assert DNATools.codon_table["TAA"] == "*"
    assert DNATools.codon_table["TAG"] == "*"
    assert DNATools.codon_table["TGA"] == "*"


def test_non_string_sequence_raises_type_error():
    with pytest.raises(TypeError):
        DNATools(None)  # type: ignore[arg-type]
    with pytest.raises(TypeError):
        DNATools.detect_mutations("ACGT", b"ACGT")  # type: ignore[arg-type]


def test_equal_length_scan_emits_no_warning():
    with warnings.catch_warnings():
        warnings.simplefilter("error")
        mutations = DNATools.detect_mutations("ac gt", "AGGT")
    assert mutations == [{"position": 2, "reference": "C", "sample": "G"}]


def test_identical_sequences_have_no_mutations():
    with warnings.catch_warnings():
        warnings.simplefilter("error")
        assert DNATools.detect_mutations("ACGT", "ACGT") == []


def test_length_difference_warns_and_reports_the_suffix():
    with pytest.warns(LengthDifferenceWarning, match="lengths differ"):
        insertion = DNATools.detect_mutations("ACGT", "AGCTA")
    assert insertion == [
        {"position": 2, "reference": "C", "sample": "G"},
        {"position": 3, "reference": "G", "sample": "C"},
        {"position": 5, "reference": "", "sample": "A"},
    ]

    with pytest.warns(LengthDifferenceWarning, match="lengths differ"):
        deletion = DNATools.detect_mutations("ACGTT", "ACGT")
    assert deletion == [{"position": 5, "reference": "T", "sample": ""}]


def test_internal_deletion_is_not_placed_as_a_gap():
    """A one-base internal deletion is reported from an ungapped scan.

    Reference ATGCAT compared with sample AGCAT is a deletion of T at
    base 2 under a gap-aware alignment. This method walks the shared
    prefix without gaps, so the shifted bases are substitutions and the
    extra reference base is the terminal suffix record.
    """
    with pytest.warns(LengthDifferenceWarning):
        mutations = DNATools.detect_mutations("ATGCAT", "AGCAT")
    assert mutations == [
        {"position": 2, "reference": "T", "sample": "G"},
        {"position": 3, "reference": "G", "sample": "C"},
        {"position": 4, "reference": "C", "sample": "A"},
        {"position": 5, "reference": "A", "sample": "T"},
        {"position": 6, "reference": "T", "sample": ""},
    ]


def test_detect_mutations_compares_non_acgt_characters():
    mutations = DNATools.detect_mutations("ACNT", "ACGT")
    assert mutations == [{"position": 3, "reference": "N", "sample": "G"}]


def test_parse_fasta_multiline_records_and_duplicate_headers(tmp_path: Path):
    fasta = tmp_path / "records.fasta"
    fasta.write_text(
        "\n".join(
            [
                ">first example",
                "ATGC",
                "TT",
                ">second",
                "AAAA",
                ">first example",
                "GG",
                "",
            ]
        ),
        encoding="utf-8",
    )
    # Sequence before the first header is ignored.
    with_preamble = tmp_path / "preamble.fasta"
    with_preamble.write_text(">only\nACGT\n", encoding="utf-8")

    parsed = DNATools.parse_fasta(str(fasta))
    assert parsed == {"first example": "GG", "second": "AAAA"}
    assert DNATools.parse_fasta(str(with_preamble)) == {"only": "ACGT"}


def test_parse_fasta_missing_file():
    with pytest.raises(FileNotFoundError):
        DNATools.parse_fasta("examples/does_not_exist.fasta")


def test_example_fasta_is_the_hbb_coding_sequence():
    records = DNATools.parse_fasta(str(EXAMPLE))
    assert len(records) == 1
    header, sequence = next(iter(records.items()))
    assert header.startswith("NM_000518.5:51-494")
    assert "NP_000509.1" in header
    assert "https://www.ncbi.nlm.nih.gov/nuccore/NM_000518.5" in header

    dna = DNATools(sequence)
    assert dna.is_valid_dna()
    assert dna.length() == 444
    assert dna.sequence.startswith("ATGGTGCATCTG")
    assert dna.sequence.endswith("AAGTATCACTAA")
    assert dna.gc_content() == 56.08
    assert {base: dna.sequence.count(base) for base in "ACGT"} == {
        "A": 88,
        "C": 113,
        "G": 136,
        "T": 107,
    }
    assert len(NP_000509_1) == 147
    assert dna.translate_to_protein() == NP_000509_1 + "*"
