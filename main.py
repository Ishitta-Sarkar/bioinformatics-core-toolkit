"""Run the toolkit on the bundled HBB coding sequence."""

from pathlib import Path

from src.dna_tools import DNATools

EXAMPLE = Path(__file__).resolve().parent / "examples" / "sample_sequence.fasta"


def _print_mutations(mutations: list[dict[str, str | int]]) -> None:
    if not mutations:
        print("No differences in the compared span.")
        return
    for mutation in mutations:
        print(
            f"Position {mutation['position']}: "
            f"{mutation['reference']!r} -> {mutation['sample']!r}"
        )


def main() -> None:
    records = DNATools.parse_fasta(str(EXAMPLE))
    header, sequence = next(iter(records.items()))
    dna = DNATools(sequence)

    print("Bioinformatics Core Toolkit")
    print("=" * 60)
    print("Record:", header)
    print("Length:", dna.length(), "bp")
    print("Valid DNA:", dna.is_valid_dna())
    print("CDS start:", dna.sequence[:12])
    print(f"GC content: {dna.gc_content():.2f}%")
    for base in "ACGT":
        count = dna.sequence.count(base)
        percent = round(count / dna.length() * 100, 2)
        print(f"{base}: {count} ({percent:.2f}%)")
    print("Protein:", dna.translate_to_protein())
    print("Transcript (first 60 nt):", dna.transcribe()[:60])
    print("Reverse complement (first 60 nt):", dna.reverse_complement()[:60])

    print()
    print("Substitution scan on a synthetic equal-length pair")
    print("-" * 60)
    reference = "ATGCGTACCGTAGCTAGC"
    sample = "ATGAGTACCGTGGCTAGC"
    print("Reference:", reference)
    print("Sample:   ", sample)
    _print_mutations(DNATools.detect_mutations(reference, sample))

    print()
    print("Length difference on a synthetic pair")
    print("-" * 60)
    reference = "ACGTACGT"
    sample = "ACGTA"
    print("Reference:", reference)
    print("Sample:   ", sample)
    _print_mutations(DNATools.detect_mutations(reference, sample))


if __name__ == "__main__":
    main()
