class DNATools:
    codon_table = {
        "ATA": "I", "ATC": "I", "ATT": "I", "ATG": "M",
        "ACA": "T", "ACC": "T", "ACG": "T", "ACT": "T",
        "AAC": "N", "AAT": "N", "AAA": "K", "AAG": "K",
        "AGC": "S", "AGT": "S", "AGA": "R", "AGG": "R",
        "CTA": "L", "CTC": "L", "CTG": "L", "CTT": "L",
        "CCA": "P", "CCC": "P", "CCG": "P", "CCT": "P",
        "CAC": "H", "CAT": "H", "CAA": "Q", "CAG": "Q",
        "CGA": "R", "CGC": "R", "CGG": "R", "CGT": "R",
        "GTA": "V", "GTC": "V", "GTG": "V", "GTT": "V",
        "GCA": "A", "GCC": "A", "GCG": "A", "GCT": "A",
        "GAC": "D", "GAT": "D", "GAA": "E", "GAG": "E",
        "GGA": "G", "GGC": "G", "GGG": "G", "GGT": "G",
        "TCA": "S", "TCC": "S", "TCG": "S", "TCT": "S",
        "TTC": "F", "TTT": "F", "TTA": "L", "TTG": "L",
        "TAC": "Y", "TAT": "Y", "TAA": "*", "TAG": "*",
        "TGC": "C", "TGT": "C", "TGA": "*", "TGG": "W"
    }

    def __init__(self, sequence):
        self.sequence = sequence.upper().replace(" ", "").replace("\n", "")

    def is_valid_dna(self):
        return all(base in "ATGC" for base in self.sequence)

    def length(self):
        return len(self.sequence)

    def gc_content(self):
        if len(self.sequence) == 0:
            return 0

        gc = self.sequence.count("G") + self.sequence.count("C")
        return round((gc / len(self.sequence)) * 100, 2)

    def transcribe(self):
        if not self.is_valid_dna():
            return "Invalid DNA sequence"
        return self.sequence.replace("T", "U")

    def reverse_complement(self):
        if not self.is_valid_dna():
            return "Invalid DNA sequence"

        complement = {
            "A": "T",
            "T": "A",
            "G": "C",
            "C": "G"
        }

        return "".join(complement[base] for base in reversed(self.sequence))

    def translate_to_protein(self):
        if not self.is_valid_dna():
            return "Invalid DNA sequence"

        protein = ""

        for i in range(0, len(self.sequence) - 2, 3):
            codon = self.sequence[i:i + 3]
            protein += self.codon_table.get(codon, "X")

        return protein

    @staticmethod
    def parse_fasta(file_path):
        sequences = {}
        current_header = None
        current_sequence = []

        with open(file_path, "r") as file:
            for line in file:
                line = line.strip()

                if line.startswith(">"):
                    if current_header:
                        sequences[current_header] = "".join(current_sequence)

                    current_header = line[1:]
                    current_sequence = []
                else:
                    current_sequence.append(line)

            if current_header:
                sequences[current_header] = "".join(current_sequence)

        return sequences

    @staticmethod
    def detect_mutations(reference_sequence, sample_sequence):
        reference_sequence = reference_sequence.upper()
        sample_sequence = sample_sequence.upper()

        mutations = []

        min_length = min(len(reference_sequence), len(sample_sequence))

        for i in range(min_length):
            if reference_sequence[i] != sample_sequence[i]:
                mutations.append({
                    "position": i + 1,
                    "reference": reference_sequence[i],
                    "sample": sample_sequence[i]
                })

        return mutations
