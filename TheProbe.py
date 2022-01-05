from TheSymbol import TheSymbol


class TheProbe:
    def __init__(self, start: int, end: int, count_seq: [TheSymbol]):
        self.start = start
        self.end = end
        self.seq = "".join([x.symbol for x in count_seq[start:end]])
