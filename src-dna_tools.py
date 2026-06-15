class DNATools:

    def __init__(self, sequence):
        self.sequence = sequence.upper()

    def length(self):
        return len(self.sequence)

    def gc_content(self):
        gc = self.sequence.count("G") + self.sequence.count("C")
        return round((gc / len(self.sequence)) * 100, 2)

    def transcribe(self):
        return self.sequence.replace("T", "U")
