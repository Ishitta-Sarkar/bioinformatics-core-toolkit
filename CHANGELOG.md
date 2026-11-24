# Changelog

All notable changes to this project are documented in this file.

The format is based on [Keep a Changelog](https://keepachangelog.com/en/1.1.0/),
and this project follows [Semantic Versioning](https://semver.org/spec/v2.0.0.html).

## [0.1.0] - 2026-11-24

First versioned release of the learning toolkit.

### Added

- `pyproject.toml`, so the package installs with `pip install -e .`.
- `requirements.txt` (no runtime dependencies) and `requirements-dev.txt`
  (pytest 9.2.0, Biopython 1.88, NumPy 2.2.6, Ruff 0.17.6).
- Pytest suite for GC content, transcription, reverse complement, translation,
  FASTA parsing, and substitution detection, cross-checked against Biopython.
- GitHub Actions workflow (`.github/workflows/ci.yml`) that lints with Ruff and
  runs the tests on Python 3.10, 3.11, and 3.12.
- `CITATION.cff`.
- Base-composition figure for the example HBB coding sequence
  (`examples/hbb_cds_base_composition.svg`).
- `LengthMismatchWarning`, raised by `detect_mutations` when the reference and
  sample lengths differ.

### Fixed

- `examples/sample_sequence.fasta` was labeled `Human_Hemoglobin_Example` but
  contained an 18 bp synthetic oligonucleotide (`ATGTTTAAAGGGCCCTAG`). It is
  now the coding sequence of NCBI RefSeq NM_000518.3 nucleotides 51–494
  (human HBB), which translates to NP_000509.1. Accessed 24 November 2026.
- `detect_mutations` now states that insertions and deletions are not called.
  A length difference emits `LengthMismatchWarning` instead of being ignored
  without comment. Reported positions are still ungapped substitutions on the
  shared prefix.

### Changed

- `main.py` reads `examples/sample_sequence.fasta` and prints the computed
  summary, a labeled constructed substitution, and the length-mismatch warning.
