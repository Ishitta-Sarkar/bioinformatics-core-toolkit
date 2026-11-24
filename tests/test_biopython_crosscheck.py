"""Cross-checks against Biopython.

GC content is compared to ``round(gc_fraction(sequence) * 100, 2)`` from
``Bio.SeqUtils``. Transcription, reverse complement, and translation are
compared to ``Bio.Seq.Seq`` using the standard genetic code. FASTA records
with unique headers are compared to ``Bio.SeqIO``.

Biopython does not provide an ungapped substitution list that matches
``detect_mutations``. Those positions are checked by walking a pair of
``Bio.Seq.Seq`` objects across the shared prefix. Pairwise alignment is not
used: it would call indels, which this method deliberately does not do.
"""

import warnings
from pathlib import Path

import pytest
from Bio import BiopythonWarning, SeqIO
from Bio.Seq import Seq
from Bio.SeqUtils import gc_fraction

from src.dna_tools import DNATools

ROOT = Path(__file__).resolve().parents[1]
BASES = "ACGT"


def _cleaned(sequence: str) -> str:
    return DNATools(sequence).sequence


def test_every_codon_matches_biopython_standard_table():
    assert len(DNATools.codon_table) == 64
    for first in BASES:
        for second in BASES:
            for third in BASES:
                codon = first + second + third
                assert DNATools.codon_table[codon] == str(Seq(codon).translate())


def test_gc_transcription_complement_and_translation_match_biopython():
    fasta = DNATools.parse_fasta(ROOT / "examples" / "sample_sequence.fasta")
    hbb = next(iter(fasta.values()))
    raw_sequences = [
        "",
        "A",
        "AT",
        "ATG",
        "ATGTT",
        "GAA",
        "GAAAAA",
        "ATGTTTAAAGGGCCCTAG",
        " atg ttt aaa ggg ccc tag\n",
        "TAATAGTGA",
        hbb.lower(),
        hbb,
    ]

    for raw in raw_sequences:
        sequence = _cleaned(raw)
        tool = DNATools(raw)
        if sequence == "":
            expected_gc = 0
        else:
            expected_gc = round(gc_fraction(sequence) * 100, 2)
        assert tool.gc_content() == expected_gc
        assert tool.transcribe() == str(Seq(sequence).transcribe())
        assert tool.reverse_complement() == str(Seq(sequence).reverse_complement())
        with warnings.catch_warnings():
            warnings.simplefilter("ignore", BiopythonWarning)
            expected_protein = str(Seq(sequence).translate())
        assert tool.translate_to_protein() == expected_protein


def test_parse_fasta_matches_seqio_for_unique_headers(tmp_path: Path):
    path = tmp_path / "unique.fasta"
    path.write_text(
        ">first record\nATGC\nTG\n>second\naa\n",
        encoding="utf-8",
    )
    ours = DNATools.parse_fasta(path)
    theirs = {
        record.description: str(record.seq) for record in SeqIO.parse(path, "fasta")
    }
    assert ours == theirs

    example = ROOT / "examples" / "sample_sequence.fasta"
    ours = DNATools.parse_fasta(example)
    theirs = {
        record.description: str(record.seq) for record in SeqIO.parse(example, "fasta")
    }
    assert ours == theirs


def test_substitutions_match_biopython_seq_positions():
    pairs = [
        ("ACGT", "AGGT"),
        ("ACGT", "TGCA"),
        ("ATGTTTAAAGGGCCCTAG", "ATGTTTAAAGGGCCCTAG"),
        ("atgc", "ATGA"),
        ("ACGT", "AGGAAA"),
        ("", "AC"),
        ("ATGCAT", "AGCATT"),
    ]
    for reference, sample in pairs:
        expected = []
        ref_seq = Seq(reference.upper())
        sample_seq = Seq(sample.upper())
        for position, (ref_base, sample_base) in enumerate(
            zip(ref_seq, sample_seq, strict=False),
            start=1,
        ):
            if ref_base != sample_base:
                expected.append(
                    {
                        "position": position,
                        "reference": str(ref_base),
                        "sample": str(sample_base),
                    }
                )
        if len(reference) == len(sample):
            with warnings.catch_warnings():
                warnings.simplefilter("error", UserWarning)
                assert DNATools.detect_mutations(reference, sample) == expected
        else:
            with pytest.warns(
                UserWarning, match="does not call insertions or deletions"
            ):
                assert DNATools.detect_mutations(reference, sample) == expected
