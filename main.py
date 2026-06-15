print("Bioinformatics Core Toolkit")
print("=" * 35)

dna_sequence = "ATGCGTACCGTAGCTAGC"

print("DNA Sequence:", dna_sequence)
print("Length:", len(dna_sequence))

gc_count = dna_sequence.count("G") + dna_sequence.count("C")
gc_content = (gc_count / len(dna_sequence)) * 100

print("GC Content:", round(gc_content, 2), "%")
print("RNA Transcript:", dna_sequence.replace("T", "U"))
