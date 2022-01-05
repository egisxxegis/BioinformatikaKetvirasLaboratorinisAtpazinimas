from TheSymbol import TheSymbol


class TheProbeStart:

    source = None  # TheSymbol type

    def __init__(self, start: int, end: int):
        self.start = start
        self.end = end
        if self.source == "":
            raise Exception("TheProbeStart has no source set")
        self.min = self.find_min(self.start, self.end)

    @staticmethod
    def find_min(start, end):
        return min([x.count for x in TheProbeStart.source[start:end]])

    @staticmethod
    def set_source(source: [TheSymbol]):
        TheProbeStart.source = source

    def get_move_right_delta(self, steps=1):
        le_right = self.source[self.end + steps - 1].count
        if le_right < self.min:
            return le_right - self.min
        if self.min == self.source[self.start].count:
            return self.find_min(self.start+steps, self.end+steps) - self.min
        return 0
