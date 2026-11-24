"""The demo and the committed figure are produced from the example FASTA."""

from __future__ import annotations

from pathlib import Path

import pytest

from main import main
from scripts.plot_base_composition import render_svg
from src.dna_tools import DNATools, LengthDifferenceWarning

ROOT = Path(__file__).resolve().parents[1]
FIGURE = ROOT / "docs" / "hbb_base_composition.svg"
EXAMPLE = ROOT / "examples" / "sample_sequence.fasta"


def test_main_prints_the_hbb_record_and_both_mutation_examples(capsys):
    with pytest.warns(LengthDifferenceWarning):
        main()
    out = capsys.readouterr().out
    assert "NM_000518.5:51-494" in out
    assert "NP_000509.1" in out
    assert "Length: 444 bp" in out
    assert "CDS start: ATGGTGCATCTG" in out
    assert "GC content: 56.08%" in out
    assert "A: 88 (19.82%)" in out
    assert "C: 113 (25.45%)" in out
    assert "G: 136 (30.63%)" in out
    assert "T: 107 (24.10%)" in out
    assert (
        "MVHLTPEEKSAVTALWGKVNVDEVGGEALGRLLVVYPWTQRFFESFGDLSTPDAVMGNPKVKAHGKKVLGAFSDGLAHLDNLKGTFATLSELHCDKLHVDPENFRLLGNVLVCVLAHHFGKEFTPPVQAAYQKVVAGVANALAHKYH*"
        in out
    )
    assert "Position 4: 'C' -> 'A'" in out
    assert "Position 12: 'A' -> 'G'" in out
    assert "Position 6: 'CGT' -> ''" in out


def test_committed_figure_matches_the_fasta():
    records = DNATools.parse_fasta(str(EXAMPLE))
    sequence = next(iter(records.values()))
    rendered = render_svg(sequence)
    assert FIGURE.read_text(encoding="utf-8") == rendered
    assert "A 88 (19.82%)" in rendered
    assert "C 113 (25.45%)" in rendered
    assert "G 136 (30.63%)" in rendered
    assert "T 107 (24.10%)" in rendered
    assert "GC content 56.08%" in rendered


def test_version_is_recorded_in_the_citation_file():
    from src import __version__

    citation = (ROOT / "CITATION.cff").read_text(encoding="utf-8")
    assert f"version: {__version__}" in citation
