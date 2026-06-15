from src.dna_tools import DNATools

sequence = DNATools("ATGCGTACCGTAGCTAGC")

print("Length:", sequence.length())
print("GC Content:", sequence.gc_content())
print("RNA:", sequence.transcribe())
