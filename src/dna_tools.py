"""Core DNA sequence operations implemented in pure Python.

The operations cover validation, GC content, transcription, reverse
complement, translation with the standard genetic code, FASTA parsing,
and ungapped substitution detection. They are a learning implementation.
Production work should use Biopython.
"""

from __future__ import annotations

import warnings
from pathlib import Path


class LengthMismatchWarning(UserWarning):
    """Reference and sample lengths differ, so indels were not called.

    ``DNATools.detect_mutations`` only reports substitutions on the shared
    prefix. It does not turn a length difference into an insertion or a
    deletion.
    """


class DNATools:
    """Operations on one DNA sequence.

    Parameters
    ----------
    sequence : str
        DNA sequence. Whitespace and newlines are removed and the sequence
        is uppercased. Characters other than A, C, G, and T are kept so
        that ``is_valid_dna`` can reject them.

    Notes
    -----
    An empty sequence after cleaning is treated as valid DNA of length 0.
    GC content of an empty sequence is 0.
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
        if not isinstance(sequence, str):
            raise TypeError("sequence must be a string")
        self.sequence = sequence.upper().replace(" ", "").replace("\n", "")

    def is_valid_dna(self) -> bool:
        """Return True when every base is A, C, G, or T.

        An empty sequence returns True. Ambiguous IUPAC bases, including N,
        return False. Uracil returns False: this class expects DNA.
        """
        return all(base in "ATGC" for base in self.sequence)

    def length(self) -> int:
        """Return the number of characters in the cleaned sequence."""
        return len(self.sequence)

    def gc_content(self) -> float:
        """Return the GC percentage, rounded to two decimal places.

        GC content is ``(count(G) + count(C)) / length * 100``, rounded with
        the built-in ``round`` function (banker's rounding). An empty
        sequence returns ``0`` (an integer).

        ``is_valid_dna`` is not required. Any character that is not G or C,
        including N, counts toward the length and not toward the GC total.
        """
        if len(self.sequence) == 0:
            return 0

        gc = self.sequence.count("G") + self.sequence.count("C")
        return round((gc / len(self.sequence)) * 100, 2)

    def transcribe(self) -> str:
        """Return the DNA sequence with T replaced by U.

        Returns the string ``"Invalid DNA sequence"`` when the sequence
        contains a base other than A, C, G, or T. This is not an RNA
        ``Seq`` object and it does not reverse the strand.
        """
        if not self.is_valid_dna():
            return "Invalid DNA sequence"
        return self.sequence.replace("T", "U")

    def reverse_complement(self) -> str:
        """Return the reverse complement (A-T, C-G), 5' to 3' on the other strand.

        Returns the string ``"Invalid DNA sequence"`` when the sequence is
        not valid DNA. An empty sequence returns an empty string.
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
        """Translate the sequence from the first base with the standard genetic code.

        Codons are read in steps of three from the start of the sequence.
        A trailing incomplete codon is dropped. Stop codons are written as
        ``*`` and translation continues after a stop. The codon table is
        NCBI translation table 1 (standard), written with DNA bases (T, not U).
        Unknown codons would be written as ``X``; every codon of valid DNA
        is in the table, so ``X`` is only a fallback.

        Returns the string ``"Invalid DNA sequence"`` when the sequence is
        not valid DNA. This method does not search for open reading frames
        and does not require the sequence to start with ATG.
        """
        if not self.is_valid_dna():
            return "Invalid DNA sequence"

        protein = ""

        for i in range(0, len(self.sequence) - 2, 3):
            codon = self.sequence[i : i + 3]
            protein += self.codon_table.get(codon, "X")

        return protein

    @staticmethod
    def parse_fasta(file_path: str | Path) -> dict[str, str]:
        """Read a FASTA file into a dictionary of header to sequence.

        Parameters
        ----------
        file_path : str or pathlib.Path
            Path of a FASTA file, opened as UTF-8.

        Returns
        -------
        dict[str, str]
            Keys are header lines without the leading ``>``. Values are the
            sequence lines concatenated, with whitespace stripped from each
            line. Case is preserved. A repeated header keeps the last record.
            Sequence lines that appear before the first header are ignored.

        Raises
        ------
        FileNotFoundError
            If ``file_path`` does not exist.

        Notes
        -----
        The parser does not interpret ``;`` comment lines, quality scores, or
        multi-line headers. A ``;`` line is treated as sequence text.
        """
        sequences: dict[str, str] = {}
        current_header = None
        current_sequence: list[str] = []

        with open(file_path, encoding="utf-8") as file:
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
    def detect_mutations(
        reference_sequence: str,
        sample_sequence: str,
    ) -> list[dict[str, int | str]]:
        """Report ungapped substitutions between a reference and a sample.

        Parameters
        ----------
        reference_sequence, sample_sequence : str
            Sequences to compare. They are uppercased and otherwise compared
            as given: whitespace is not removed (unlike ``DNATools.__init__``).

        Returns
        -------
        list of dict
            One entry per differing position on the shared prefix. Each entry
            has ``position`` (1-based index into that prefix), ``reference``,
            and ``sample``. Positions past the shorter sequence are omitted.

        Raises
        ------
        TypeError
            If either argument is not a string.

        Warns
        -----
        LengthMismatchWarning
            If the two strings have different lengths. The return value is
            still only the substitutions on the shared prefix. The length
            difference is not reported as an insertion or a deletion.

        Notes
        -----
        This is not variant calling. Insertions, deletions, and other
        alignment gaps are ignored. Two sequences of equal length can still
        hide a compensating insertion and deletion, because those would look
        like a run of substitutions (or like nothing, if the comparison is
        not aligned). Align the sequences before treating these positions as
        biological variants.

        Examples
        --------
        >>> DNATools.detect_mutations("ACGT", "AGGT")
        [{'position': 2, 'reference': 'C', 'sample': 'G'}]
        """
        if not isinstance(reference_sequence, str) or not isinstance(
            sample_sequence, str
        ):
            raise TypeError("reference_sequence and sample_sequence must be strings")

        reference_sequence = reference_sequence.upper()
        sample_sequence = sample_sequence.upper()

        if len(reference_sequence) != len(sample_sequence):
            warnings.warn(
                (
                    f"Reference length ({len(reference_sequence)}) and sample length "
                    f"({len(sample_sequence)}) differ. detect_mutations reports "
                    "substitutions on the shared prefix only and does not call "
                    "insertions or deletions."
                ),
                LengthMismatchWarning,
                stacklevel=2,
            )

        mutations: list[dict[str, int | str]] = []
        shared_length = min(len(reference_sequence), len(sample_sequence))

        for i in range(shared_length):
            if reference_sequence[i] != sample_sequence[i]:
                mutations.append(
                    {
                        "position": i + 1,
                        "reference": reference_sequence[i],
                        "sample": sample_sequence[i],
                    }
                )

        return mutations
