"""Provenance checks for the bundled HBB coding sequence and the demo."""

import importlib
import importlib.util
import subprocess
import sys
from pathlib import Path

from src.dna_tools import DNATools

ROOT = Path(__file__).resolve().parents[1]
FASTA_PATH = ROOT / "examples" / "sample_sequence.fasta"
SVG_PATH = ROOT / "examples" / "hbb_cds_base_composition.svg"

# Translation qualifier of NM_000518.3 (protein NP_000509.1), copied from the
# GenBank file downloaded from NCBI on 24 November 2026. The CDS translation
# in that file does not include the terminal stop.
NCBI_NP_000509_1 = (
    "MVHLTPEEKSAVTALWGKVNVDEVGGEALGRLLVVYPWTQRFFESFGDLSTPDAVMGNPKVKAHGKKVLGAFS"
    "DGLAHLDNLKGTFATLSELHCDKLHVDPENFRLLGNVLVCVLAHHFGKEFTPPVQAAYQKVVAGVANALAHKYH"
)


def _example_record() -> tuple[str, str]:
    records = DNATools.parse_fasta(FASTA_PATH)
    assert len(records) == 1
    return next(iter(records.items()))


def test_example_fasta_is_the_cited_hbb_cds():
    header, sequence = _example_record()
    assert "Human_Hemoglobin_Example" not in header
    assert "NM_000518.3:51-494" in header
    assert "NP_000509.1" in header
    assert "GeneID:3043" in header
    assert "accessed 2026-11-24" in header
    assert sequence.startswith("ATGGTGCACCTG")
    assert sequence.endswith("CACTAA")
    assert len(sequence) == 444
    assert set(sequence) <= set("ACGT")
    assert " " not in sequence

    dna = DNATools(sequence)
    assert dna.is_valid_dna() is True
    assert dna.length() == 444
    assert dna.gc_content() == 56.31
    assert {base: sequence.count(base) for base in "ACGT"} == {
        "A": 88,
        "C": 114,
        "G": 136,
        "T": 106,
    }
    assert len(NCBI_NP_000509_1) == 147
    assert dna.translate_to_protein() == NCBI_NP_000509_1 + "*"


def test_main_constant_matches_the_ncbi_protein():
    import main

    assert main.NP_000509_1 == NCBI_NP_000509_1


def test_base_composition_figure_is_generated_from_the_fasta():
    script = ROOT / "examples" / "make_base_composition_figure.py"
    spec = importlib.util.spec_from_file_location(
        "make_base_composition_figure", script
    )
    module = importlib.util.module_from_spec(spec)
    assert spec.loader is not None
    spec.loader.exec_module(module)

    _, sequence = _example_record()
    dna = DNATools(sequence)
    counts = {base: sequence.count(base) for base in "ACGT"}
    svg = module.render_hbb_base_composition_svg(counts, dna.gc_content())
    assert svg == SVG_PATH.read_text(encoding="utf-8")
    assert "<!-- bases A=88 C=114 G=136 T=106 gc=56.31 n=444 -->" in svg
    assert "<script" not in svg


def test_importing_main_does_not_run_the_demo(capsys):
    import main

    importlib.reload(main)
    captured = capsys.readouterr()
    assert captured.out == ""
    assert captured.err == ""


def test_main_script_prints_the_computed_demo():
    result = subprocess.run(
        [sys.executable, str(ROOT / "main.py")],
        check=False,
        capture_output=True,
        text=True,
    )
    assert result.returncode == 0
    assert result.stderr == ""
    output = result.stdout
    assert "NM_000518.3:51-494" in output
    assert "Length: 444" in output
    assert "Valid DNA: True" in output
    assert "GC content: 56.31 %" in output
    assert "AUGGUGCACCUG" in output
    assert NCBI_NP_000509_1 + "*" in output
    assert "Matches NP_000509.1 plus terminal stop: True" in output
    assert "not a biological variant" in output
    assert "Reference: ATGGTGCACCTG" in output
    assert "Sample:    ATGGTGAACCTG" in output
    assert "Position: 7 | Reference: C | Sample: A" in output
    assert "Reference length: 12" in output
    assert "Sample length: 11" in output
    assert "does not call insertions or deletions" in output
    assert "Human_Hemoglobin_Example" not in output
