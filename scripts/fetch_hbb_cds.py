"""Write examples/sample_sequence.fasta from NCBI Nucleotide NM_000518.5.

The download is the coding sequence reported by NCBI as
``rettype=fasta_cds_na`` for NM_000518.5 (CDS 51..494, protein
NP_000509.1). The script checks that header before it replaces the
example file.
"""

from __future__ import annotations

import urllib.request
from pathlib import Path

EUTILS_URL = (
    "https://eutils.ncbi.nlm.nih.gov/entrez/eutils/efetch.fcgi"
    "?db=nuccore&id=NM_000518.5&rettype=fasta_cds_na&retmode=text"
    "&tool=bioinformatics-core-toolkit&email=ishittas4@gmail.com"
)
ROOT = Path(__file__).resolve().parents[1]
OUT = ROOT / "examples" / "sample_sequence.fasta"
HEADER = (
    ">NM_000518.5:51-494 Homo sapiens hemoglobin subunit beta (HBB) "
    "coding sequence | protein NP_000509.1 | CCDS7753.1 | "
    "https://www.ncbi.nlm.nih.gov/nuccore/NM_000518.5 | "
    "retrieved 2026-11-24"
)
REQUIRED_HEADER_TOKENS = ("NM_000518.5", "NP_000509.1", "location=51..494")


def fetch_cds(url: str = EUTILS_URL) -> str:
    """Download the NM_000518.5 coding sequence and return the DNA string."""
    request = urllib.request.Request(
        url,
        headers={
            "User-Agent": "bioinformatics-core-toolkit (ishittas4@gmail.com)",
        },
    )
    with urllib.request.urlopen(request, timeout=60) as response:
        text = response.read().decode("utf-8")

    lines = [line.strip() for line in text.splitlines() if line.strip()]
    if not lines or not lines[0].startswith(">"):
        raise RuntimeError("NCBI did not return a FASTA record")

    header = lines[0]
    missing = [token for token in REQUIRED_HEADER_TOKENS if token not in header]
    if missing:
        raise RuntimeError(f"NCBI header missing {missing}: {header}")

    sequence = "".join(lines[1:]).upper()
    if not set(sequence) <= set("ACGT"):
        raise RuntimeError("CDS contains a base other than A, C, G, or T")
    if len(sequence) != 444:
        raise RuntimeError(f"Expected 444 bp from CDS 51..494, got {len(sequence)}")
    if not sequence.startswith("ATGGTGCATCTG") or not sequence.endswith("TAA"):
        raise RuntimeError("CDS does not have the HBB start codon and terminal TAA")
    return sequence


def write_fasta(sequence: str, path: Path = OUT) -> None:
    """Write ``sequence`` with the repository's provenance header."""
    wrapped = "\n".join(sequence[i : i + 60] for i in range(0, len(sequence), 60))
    path.write_text(f"{HEADER}\n{wrapped}\n", encoding="utf-8")


def main() -> None:
    sequence = fetch_cds()
    write_fasta(sequence)
    print(f"Wrote {OUT} ({len(sequence)} bp)")


if __name__ == "__main__":
    main()
