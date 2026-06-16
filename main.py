from src.dna_tools import DNATools

dna = DNATools("ATGTTTAAAGGG")

print("Bioinformatics Core Toolkit")
print("=" * 35)

print("DNA Sequence:", dna.sequence)
print("Length:", dna.length())
print("GC Content:", dna.gc_content(), "%")
print("RNA Transcript:", dna.transcribe())
print("Reverse Complement:", dna.reverse_complement())
print("Protein Translation:", dna.translate_to_protein())

reference = "ATGCGTACCGTAGCTAGC"
sample = "ATGAGTACCGTGGCTAGC"

mutations = DNATools.detect_mutations(reference, sample)

print("\nMutation Detection")
print("-" * 35)

for mutation in mutations:
    print(
        "Position:",
        mutation["position"],
        "| Reference:",
        mutation["reference"],
        "| Sample:",
        mutation["sample"]
    )
