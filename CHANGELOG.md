# Changelog

## 0.2.0 - 2026-11-24

### Changed

- Replaced `examples/sample_sequence.fasta`. The previous record was headed `>Human_Hemoglobin_Example` and contained the 18 bp synthetic oligo `ATGTTTAAAGGGCCCTAG`. The file is now the coding sequence of human HBB, NCBI Nucleotide [NM_000518.5](https://www.ncbi.nlm.nih.gov/nuccore/NM_000518.5) positions 51–494 (444 bp), protein [NP_000509.1](https://www.ncbi.nlm.nih.gov/protein/NP_000509.1), retrieved 24 November 2026.
- `detect_mutations` warns with `LengthDifferenceWarning` when the reference and sample lengths differ, and appends one record for the unmatched suffix. Equal-length substitution records are unchanged. The scan still does not insert gaps, so an internal indel shifts later bases into substitutions.
- `DNATools` and `detect_mutations` remove all whitespace, not only spaces and newlines.
- `main.py` reads `examples/sample_sequence.fasta` and prints the HBB summary plus two synthetic mutation examples.

### Added

- `pyproject.toml`, so `pip install -e .` installs the `src` package.
- `requirements.txt` (no runtime dependencies) and `requirements-dev.txt` (pytest 9.2.0, Biopython 1.88, Ruff 0.17.6).
- A pytest suite. Sequence operations on the HBB coding sequence, on all 64 codons, and on the worked 18 bp example are compared with Biopython 1.88.
- GitHub Actions workflow `.github/workflows/ci.yml` (Ruff and pytest on Python 3.10–3.13).
- `CITATION.cff`.
- Base-composition figure `docs/hbb_base_composition.svg`, written by `scripts/plot_base_composition.py`.
- `scripts/fetch_hbb_cds.py`, which re-downloads the coding sequence from NCBI E-utilities.

## 0.1.0 - 2026-06-16

Initial public version: GC content, transcription, reverse complement, translation with the standard genetic code, FASTA parsing, and pairwise substitution detection.
