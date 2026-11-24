"""Demo of DNATools on the bundled HBB coding sequence.

Run from any working directory:

    python main.py
"""

from __future__ import annotations

import warnings
from pathlib import Path

from src.dna_tools import DNATools, LengthMismatchWarning

# NCBI NP_000509.1, from the CDS translation on NM_000518.3 (accessed
# 24 November 2026). translate_to_protein keeps the terminal stop codon.
NP_000509_1 = (
    "MVHLTPEEKSAVTALWGKVNVDEVGGEALGRLLVVYPWTQRFFESFGDLSTPDAVMGNPKVKAHGKKVLGAFS"
    "DGLAHLDNLKGTFATLSELHCDKLHVDPENFRLLGNVLVCVLAHHFGKEFTPPVQAAYQKVVAGVANALAHKYH"
)


def main() -> None:
    """Print sequence statistics, one constructed substitution, and a length warning."""
    fasta_path = Path(__file__).resolve().parent / "examples" / "sample_sequence.fasta"
    records = DNATools.parse_fasta(fasta_path)
    if len(records) != 1:
        raise SystemExit(f"Expected one FASTA record, found {len(records)}")

    header, sequence = next(iter(records.items()))
    dna = DNATools(sequence)
    protein = dna.translate_to_protein()

    print("Bioinformatics Core Toolkit")
    print("=" * 35)
    print("Record:", header)
    print("Length:", dna.length())
    print("Valid DNA:", dna.is_valid_dna())
    print("GC content:", dna.gc_content(), "%")
    print("RNA transcript (first 60 nt):", dna.transcribe()[:60])
    print("Protein translation:", protein)
    print("Protein length including stop:", len(protein))
    print("Matches NP_000509.1 plus terminal stop:", protein == NP_000509_1 + "*")

    # Built from the first 12 bp of the example record. The edited base is
    # a demonstration, not a reported HBB variant.
    print("\nConstructed substitution (not a biological variant)")
    print("-" * 35)
    reference = dna.sequence[:12]
    sample = reference[:6] + "A" + reference[7:]
    print("Reference:", reference)
    print("Sample:   ", sample)
    for mutation in DNATools.detect_mutations(reference, sample):
        print(
            "Position:",
            mutation["position"],
            "| Reference:",
            mutation["reference"],
            "| Sample:",
            mutation["sample"],
        )

    print("\nLength mismatch (substitutions on the shared prefix only)")
    print("-" * 35)
    short_sample = sample[:-1]
    print("Reference length:", len(reference))
    print("Sample length:", len(short_sample))
    with warnings.catch_warnings(record=True) as caught:
        warnings.simplefilter("always", LengthMismatchWarning)
        mutations = DNATools.detect_mutations(reference, short_sample)
    for item in caught:
        print("Warning:", item.message)
    if not caught:
        print("Warning: expected a length-mismatch warning, but none was issued")
    for mutation in mutations:
        print(
            "Position:",
            mutation["position"],
            "| Reference:",
            mutation["reference"],
            "| Sample:",
            mutation["sample"],
        )


if __name__ == "__main__":
    main()
