"""Elementary DNA operations implemented in pure Python.

The constructor stores one sequence. Class methods read FASTA files and
compare two sequences. The genetic code is NCBI translation table 1.
"""

from __future__ import annotations

import warnings


class LengthDifferenceWarning(UserWarning):
    """The reference and sample lengths differ.

    ``detect_mutations`` compared the shared prefix and reported the
    unmatched suffix. It did not place a gap for an internal indel.
    """


def _normalize(sequence: str) -> str:
    """Uppercase ``sequence`` and remove whitespace."""
    if not isinstance(sequence, str):
        raise TypeError(f"sequence must be a str, got {type(sequence).__name__}")
    return "".join(sequence.split()).upper()


class DNATools:
    """DNA sequence operations for one cleaned sequence.

    The constructor uppercases the sequence and removes whitespace.
    ``transcribe``, ``reverse_complement``, and ``translate_to_protein``
    return the string ``"Invalid DNA sequence"`` when any character is
    outside A, C, G, and T. They do not raise.

    Attributes:
        sequence: The cleaned sequence.
        codon_table: NCBI standard genetic code. Stops are ``*``.
    """

    codon_table = {
        "ATA": "I",
        "ATC": "I",
        "ATT": "I",
        "ATG": "M",
        "ACA": "T",
        "ACC": "T",
        "ACG": "T",
        "ACT": "T",
        "AAC": "N",
        "AAT": "N",
        "AAA": "K",
        "AAG": "K",
        "AGC": "S",
        "AGT": "S",
        "AGA": "R",
        "AGG": "R",
        "CTA": "L",
        "CTC": "L",
        "CTG": "L",
        "CTT": "L",
        "CCA": "P",
        "CCC": "P",
        "CCG": "P",
        "CCT": "P",
        "CAC": "H",
        "CAT": "H",
        "CAA": "Q",
        "CAG": "Q",
        "CGA": "R",
        "CGC": "R",
        "CGG": "R",
        "CGT": "R",
        "GTA": "V",
        "GTC": "V",
        "GTG": "V",
        "GTT": "V",
        "GCA": "A",
        "GCC": "A",
        "GCG": "A",
        "GCT": "A",
        "GAC": "D",
        "GAT": "D",
        "GAA": "E",
        "GAG": "E",
        "GGA": "G",
        "GGC": "G",
        "GGG": "G",
        "GGT": "G",
        "TCA": "S",
        "TCC": "S",
        "TCG": "S",
        "TCT": "S",
        "TTC": "F",
        "TTT": "F",
        "TTA": "L",
        "TTG": "L",
        "TAC": "Y",
        "TAT": "Y",
        "TAA": "*",
        "TAG": "*",
        "TGC": "C",
        "TGT": "C",
        "TGA": "*",
        "TGG": "W",
    }

    def __init__(self, sequence: str) -> None:
        self.sequence = _normalize(sequence)

    def is_valid_dna(self) -> bool:
        """Return whether every base is A, C, G, or T.

        An empty sequence returns ``True``.
        """
        return all(base in "ATGC" for base in self.sequence)

    def length(self) -> int:
        """Return the number of characters in the cleaned sequence."""
        return len(self.sequence)

    def gc_content(self) -> float:
        """Return the GC percentage, rounded to two decimal places.

        The value is ``100 * (count(G) + count(C)) / length``, rounded
        with Python's ``round`` (ties to even). An empty sequence
        returns ``0``. Bases other than G and C count toward the length
        and do not count as GC.
        """
        if len(self.sequence) == 0:
            return 0

        gc = self.sequence.count("G") + self.sequence.count("C")
        return round((gc / len(self.sequence)) * 100, 2)

    def transcribe(self) -> str:
        """Return the sequence with T replaced by U.

        Returns ``"Invalid DNA sequence"`` when ``is_valid_dna`` is false.
        """
        if not self.is_valid_dna():
            return "Invalid DNA sequence"
        return self.sequence.replace("T", "U")

    def reverse_complement(self) -> str:
        """Return the reverse complement (A-T and C-G).

        Returns ``"Invalid DNA sequence"`` when ``is_valid_dna`` is false.
        """
        if not self.is_valid_dna():
            return "Invalid DNA sequence"

        complement = {
            "A": "T",
            "T": "A",
            "G": "C",
            "C": "G",
        }

        return "".join(complement[base] for base in reversed(self.sequence))

    def translate_to_protein(self) -> str:
        """Translate the sequence from the first base in non-overlapping codons.

        Stops are ``*``. A codon that is not in ``codon_table`` becomes
        ``X``. A trailing incomplete codon is omitted. Returns
        ``"Invalid DNA sequence"`` when ``is_valid_dna`` is false.
        """
        if not self.is_valid_dna():
            return "Invalid DNA sequence"

        protein = ""

        for i in range(0, len(self.sequence) - 2, 3):
            codon = self.sequence[i : i + 3]
            protein += self.codon_table.get(codon, "X")

        return protein

    @staticmethod
    def parse_fasta(file_path: str) -> dict[str, str]:
        """Read a FASTA file into a dictionary of header to sequence.

        The header is the text after ``>``, with surrounding whitespace
        removed. Sequence lines are concatenated. A repeated header
        replaces the earlier sequence. Sequence text that appears before
        the first header is ignored. The file is read as UTF-8.
        """
        sequences: dict[str, str] = {}
        current_header = None
        current_sequence: list[str] = []

        with open(file_path, encoding="utf-8") as file:
            for line in file:
                line = line.strip()

                if line.startswith(">"):
                    if current_header is not None:
                        sequences[current_header] = "".join(current_sequence)

                    current_header = line[1:]
                    current_sequence = []
                else:
                    current_sequence.append(line)

            if current_header is not None:
                sequences[current_header] = "".join(current_sequence)

        return sequences

    @staticmethod
    def detect_mutations(
        reference_sequence: str,
        sample_sequence: str,
    ) -> list[dict[str, str | int]]:
        """Compare two sequences from the 5' end and list differences.

        Both sequences are uppercased and whitespace is removed. Positions
        are 1-based. Each substitution in the shared prefix is a dictionary
        ``{"position", "reference", "sample"}`` whose ``reference`` and
        ``sample`` values are single characters.

        When the lengths differ, a :class:`LengthDifferenceWarning` is
        issued and one extra dictionary is appended. Its ``position`` is
        the first unmatched base. ``reference`` or ``sample`` holds the
        unmatched suffix and the other value is an empty string. A longer
        sample suffix is the sample tail that has no reference base under
        it. A longer reference suffix is the reference tail that has no
        sample base under it.

        The scan does not insert gaps. An internal insertion or deletion
        shifts every following base, and those shifted bases are reported
        as substitutions, followed by the terminal suffix record when the
        lengths differ. Alignment-based indel calling is not implemented.

        The comparison does not require the sequences to be valid DNA.
        ``N`` and other characters are compared as stored after cleaning.
        """
        reference_sequence = _normalize(reference_sequence)
        sample_sequence = _normalize(sample_sequence)

        mutations: list[dict[str, str | int]] = []
        shared = min(len(reference_sequence), len(sample_sequence))

        for i in range(shared):
            if reference_sequence[i] != sample_sequence[i]:
                mutations.append(
                    {
                        "position": i + 1,
                        "reference": reference_sequence[i],
                        "sample": sample_sequence[i],
                    }
                )

        if len(reference_sequence) != len(sample_sequence):
            warnings.warn(
                (
                    "sequence lengths differ "
                    f"(reference {len(reference_sequence)} bp, "
                    f"sample {len(sample_sequence)} bp). "
                    f"Substitutions were called on the shared prefix of "
                    f"{shared} bp. The unmatched suffix is reported at "
                    f"position {shared + 1}. This comparison does not "
                    "insert gaps, so an internal indel shifts every "
                    "following base into a substitution."
                ),
                LengthDifferenceWarning,
                stacklevel=2,
            )
            mutations.append(
                {
                    "position": shared + 1,
                    "reference": reference_sequence[shared:],
                    "sample": sample_sequence[shared:],
                }
            )

        return mutations
