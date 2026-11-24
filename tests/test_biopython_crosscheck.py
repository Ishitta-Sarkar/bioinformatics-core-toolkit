"""Cross-checks against Biopython 1.88.

GC content is compared after both sides are rounded to two decimal places,
which is the precision ``DNATools.gc_content`` returns. Translation is
compared on sequences whose length is a multiple of three. Biopython 1.88
also translates a trailing partial codon and emits BiopythonWarning; the
partial-codon test records that both libraries return the same residues.
"""

from __future__ import annotations

from pathlib import Path

import pytest
from Bio import BiopythonWarning
from Bio.Seq import Seq
from Bio.SeqUtils import gc_fraction

from src.dna_tools import DNATools

EXAMPLE = Path(__file__).resolve().parents[1] / "examples" / "sample_sequence.fasta"


def _example_sequence() -> str:
    records = DNATools.parse_fasta(str(EXAMPLE))
    return next(iter(records.values()))


def _assert_agrees(sequence: str) -> None:
    dna = DNATools(sequence)
    bio = Seq(dna.sequence)
    assert dna.transcribe() == str(bio.transcribe())
    assert dna.reverse_complement() == str(bio.reverse_complement())
    assert dna.gc_content() == round(gc_fraction(bio) * 100, 2)
    if len(dna.sequence) % 3 == 0:
        assert dna.translate_to_protein() == str(bio.translate())


def test_example_oligo_matches_biopython():
    _assert_agrees("ATGTTTAAAGGGCCCTAG")


def test_hbb_cds_matches_biopython():
    sequence = _example_sequence()
    _assert_agrees(sequence)
    assert DNATools(sequence).translate_to_protein().endswith("*")
    assert str(Seq(sequence).translate()) == DNATools(sequence).translate_to_protein()


def test_every_codon_matches_biopython():
    for a in "ACGT":
        for b in "ACGT":
            for c in "ACGT":
                codon = a + b + c
                assert DNATools(codon).translate_to_protein() == str(
                    Seq(codon).translate()
                )
                _assert_agrees(codon)


def test_fixed_sequences_match_biopython():
    sequences = [
        "A",
        "C",
        "G",
        "T",
        "GC",
        "AT",
        "ATGTT",
        "ACGT" * 40,
        "GGGCCCCAAAATTT",
        _example_sequence()[:100],
        _example_sequence()[50:200],
    ]
    for sequence in sequences:
        dna = DNATools(sequence)
        bio = Seq(dna.sequence)
        assert dna.gc_content() == round(gc_fraction(bio) * 100, 2)
        if dna.is_valid_dna() and dna.sequence:
            assert dna.transcribe() == str(bio.transcribe())
            assert dna.reverse_complement() == str(bio.reverse_complement())


def test_partial_codon_residues_match_biopython_and_biopython_warns():
    with pytest.warns(BiopythonWarning, match="Partial codon"):
        biopython_protein = str(Seq("ATGTT").translate())
    assert DNATools("ATGTT").translate_to_protein() == biopython_protein == "M"
