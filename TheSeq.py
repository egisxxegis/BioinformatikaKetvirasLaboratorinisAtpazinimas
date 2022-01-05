class TheSeq:
    def __init__(self, name: str, seq: str):
        self.name = name
        self.seq = seq.replace("\r", "").replace("\n", "").upper()
