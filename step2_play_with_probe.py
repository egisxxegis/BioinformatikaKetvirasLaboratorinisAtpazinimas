from TheProbe import TheProbe
from TheSymbol import TheSymbol

file_name = "data/probes_6.txt"

probe_max_mismatches = 3
probe_min_start_matches = 5
gap_symbol = "-"
probe_gap_means_evil_in_detection = True


def get_probes():
    file = open(file_name, "r")
    step = 0
    the_len = 0
    start = 0
    end = 0
    probes = []
    while line := file.readline():
        if line[0] == "#":
            continue
        if step == 0:
            the_len = int(line)
            step += 1
        elif step == 1:
            start = int(line) - 1
            step += 1
        elif step == 2:
            end = int(line)
            step += 1
        elif step == 3:
            seq = line.replace("\r", "").replace("\n", "")
            probes.append(TheProbe(start, end,  [TheSymbol("z", 999)] * end))
            probes[-1].seq = seq
            if len(probes) < the_len:
                step = 1
                continue
        continue
    return probes


def is_sequence_detected(seq: [chr], probe: TheProbe):
    seq2 = seq[probe.start:probe.end]
    # start matches?
    # print(seq2)
    # print(probe.seq)
    if not seq2[0:probe_min_start_matches] == probe.seq[0:probe_min_start_matches]:
        return False
    mismatches = 0
    for i in range(probe_min_start_matches, len(probe.seq)):
        if seq2[i] != probe.seq[i]:
            mismatches += 1
        if seq2[i] == gap_symbol and probe_gap_means_evil_in_detection:
            mismatches += probe_max_mismatches
    # mismatches exceeded?
    if mismatches > probe_max_mismatches:
        return False
    return True


if __name__ == "__main__":
    z_probes = get_probes()
    while True:
        the_seq = input("Input sequence for check\n").replace("\r", "").replace("\n", "").upper()
        z_results = [is_sequence_detected(the_seq, x) for x in z_probes]
        if z_results.count(True) > 0:
            print("----Positive result")
        else:
            print("--Negative result")
