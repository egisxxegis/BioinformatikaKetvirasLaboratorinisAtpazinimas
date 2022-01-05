import os
import utils
import math

from TheSeq import TheSeq
from TheSymbolsCounter import TheSymbolsCounter
from TheSymbol import TheSymbol
from TheProbeStart import TheProbeStart
from TheProbe import TheProbe

# settings
source_file = "data/maffta_5_done_with_external_tool.fasta.txt"
probe_max_region_bp = 200
probe_max_mismatches = 3
probe_min_start_matches = 5
probe_min_length = 19
probe_max_length = 20
probe_target_percentage = 100  # how much sequences we want to cover
gap_symbol = "-"
probe_gap_means_evil_in_detection = True
# probe_initial_sensitivity = 80  # 0 - 100. 100 reacts to slightest change. 0 ignores everything


# code below
def do_check_settings():
    if probe_min_length > probe_max_length:
        raise Exception("probe min length can not be longer than max length")
    if probe_min_length < probe_min_start_matches:
        raise Exception("probe min length can not be shorter than min start matches")
    if probe_min_length > probe_max_region_bp:
        raise Exception("probe min length can not be longer than probe max region bp")
    return True


def do_check_source():
    if not os.path.isfile(source_file):
        raise Exception(f"source file was not found at ''{os.getcwd()}{source_file}''")
    file = open(source_file, "r")
    head = file.readline()
    body = file.readline()
    if not head[0] == ">":
        raise Exception("is source file even in a FASTA format? No symbol ''>'' at the beginning found")
    body = body.strip("\t\r\n").upper()
    if len(body) != body.count("G") + body.count("A") + body.count("C") + body.count("T") + body.count(gap_symbol):
        raise Exception("are your sequences expressed in GACT nucleotides and gaps? Make sure they are")
    file.close()
    return True


def get_counts(seqs: [TheSeq], split_r=None):
    split_r = split_r if split_r is not None else [0, len(seqs[0].seq)]
    freqs = [TheSymbolsCounter()] * len(seqs[0].seq)
    # first time create
    for i in range(len(freqs)):
        freqs[i] = TheSymbolsCounter()
    # second time fill-in
    for seq in seqs:
        for i in range(split_r[0], split_r[1]):  # len seq.seq
            symbol = seq.seq[i]
            if symbol == "A":
                freqs[i].a_count += 1
            elif symbol == "G":
                freqs[i].g_count += 1
            elif symbol == "T":
                freqs[i].t_count += 1
            elif symbol == "C":
                freqs[i].c_count += 1
            elif symbol == gap_symbol:
                freqs[i].gap_count += 1
    return freqs


def get_greedy_count(freqs: [TheSymbolsCounter], split_r=None):
    split_r = split_r if split_r is not None else [0, len(freqs)]
    gfreqs = [TheSymbol("-", 420)] * len(freqs)
    for i in range(split_r[0], split_r[1]):  # len freqs
        i_freq = freqs[i]
        c = TheSymbol("A", i_freq.a_count)
        if i_freq.g_count > c.count:
            c = TheSymbol("G", i_freq.g_count)
        if i_freq.t_count > c.count:
            c = TheSymbol("T", i_freq.t_count)
        if i_freq.c_count > c.count:
            c = TheSymbol("C", i_freq.c_count)
        gfreqs[i] = c
    return gfreqs


# def get_greedy_probe_starts(counts: [TheSymbol], threshold_delta):
#     init_probes = []
#     TheProbeStart.set_source(counts)
#     probe = TheProbeStart(0, probe_min_start_matches)
#     steps = 1
#     for i in range(1, len(counts) - probe_min_start_matches + 1):
#         delta = probe.get_move_right_delta(steps)
#         if delta < threshold_delta:
#             init_probes.append(probe)
#             probe = TheProbeStart(i + (steps-1), i + (steps-1) + probe_min_start_matches)
#             steps = 1
#             continue
#         if delta < 0 and steps < 5:
#             steps += 1
#             continue
#         if steps >= 5:
#             init_probes.append(probe)
#             probe = TheProbeStart(i + (steps - 1), i + (steps - 1) + probe_min_start_matches)
#             steps = 1
#             continue
#
#         probe.start += steps
#         probe.end += steps
#         probe.min += delta
#         steps = 1
#     init_probes.append(probe)
#     return init_probes


def get_stamp_min(source: [TheSymbol], slots: int, split_r=None):
    split_r = split_r if split_r is not None else [0, len(source)]
    source = [x.count for x in source]
    stamp_min = [max(source)*2] * (len(source) - slots + 1)
    buffer = source[split_r[0]:split_r[0]+slots]
    the_min = min(buffer)
    stamp_min[split_r[0]] = the_min
    for i in range(split_r[0]+1, split_r[1] - slots + 1):  # len source
        le_right = source[i+slots-1]
        if buffer[0] == the_min:
            buffer = source[i:i+slots]
            the_min = min(buffer)
        elif le_right < the_min:
            the_min = le_right
            buffer = source[i:i+slots]
        else:  # min not lost
            buffer = source[i:i+slots]
        stamp_min[i] = the_min
    return stamp_min


def get_max_probe(min_stamp: [int], counts: [TheSymbol], split_r=None):
    split_r = split_r if split_r is not None else [None, len(min_stamp)]
    the_max = max(min_stamp[split_r[0]:split_r[1]-probe_min_length])
    ind = min_stamp.index(the_max)
    return TheProbe(ind, ind + probe_min_length, counts)


def is_sequence_detected(seq: [chr], probe: TheProbe):
    seq2 = seq[probe.start:probe.end]
    # start matches?
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


def write_probes(probes: [TheProbe]):
    file = open("data/probes_6.txt", "w")
    file.write(f"# first number - amount of probes; nextly, each probe is described with start, end, seq\n")
    file.write(f"{len(probes)}\n")
    for probe in probes:
        file.write(f"{probe.start+1}\n{probe.end}\n")
        file.write(f"{probe.seq}\n")
    file.close()


def main():

    probes = []

    print("step1 starting.")
    do_check_settings()
    print("settings checked.")
    do_check_source()
    print("file beginning checked.")
    print("starting analysis.")
    seqs = utils.read_fasta(source_file)

    target_sequences = math.ceil(len(seqs) / 100 * (100-probe_target_percentage))
    print("------")
    print("length of first sequence = " + str(len(seqs[0].seq)) + " .")
    print("number of sequences = " + str(len(seqs)))
    print(f"targeted amount of sequences to remain undetected = {target_sequences}")
    print("------")

    print("getting counts.")
    counts = get_counts(seqs)
    print("counts got. getting greedy count")
    g_count = get_greedy_count(counts)
    g_count_max = max([x.count for x in g_count])
    print(f"greedy count got. max count: {g_count_max} / {len(seqs)}")
    # g_probe_delta = g_count_max/100*(100-probe_initial_sensitivity)*-1
    # print(f"getting greedy probes init with delta = {g_probe_delta}")
    # g_probe_starts = get_greedy_probe_starts(g_count, g_probe_delta)
    # print(f"got {len(g_probe_starts)} init probes")
    print("------")
    # print("counts:")
    # print(f"{[x.count for x in g_count]}")
    # print("probes init:")
    # print(f"{[x.start for x in g_probe_starts]}")
    # for x in [x.start for x in g_probe_starts]:
    #     print(f"counts [{x}:{x+10}]: {[y.count for y in g_count[x:x+10]]}")
    # print(f"{[TheProbeStart.find_min(i, i+5) for i in range(len(g_count) - 5)]}")

    print("getting min greedy stamps")
    g_count_stamps = get_stamp_min(g_count, probe_min_length)
    # print(g_count_stamps)
    print("------")
    print("getting max probe")
    max_probe = get_max_probe(g_count_stamps, g_count)
    print(f"got probe at {max_probe.start}: {max_probe.seq}")

    probes.append(max_probe)

    # cherry picking
    max_end_index = max_probe.start + probe_max_region_bp
    count_all_seqs = len(seqs)

    while True:
        last_seq_len = len(seqs)
        seqs = [x for x in seqs if not is_sequence_detected(x.seq, probes[-1])]
        print(f"Last one probe helped to identify +{last_seq_len - len(seqs)} unique cases. \n"
              f"Identified {count_all_seqs-len(seqs)}/{count_all_seqs} sequences. "
              f"Still ({len(seqs) - target_sequences}) seqs to go for {probe_target_percentage} % rate")
        if len(seqs) <= target_sequences:
            print("All seqs detected. Work is done.")
            write_probes(probes)
            print("Probes wrote into file at /data")
            return

        the_slice = [max_probe.start, max_end_index]

        print("getting counts.")
        counts = get_counts(seqs, the_slice)
        print("counts got. getting greedy count")
        g_count = get_greedy_count(counts, the_slice)
        g_count_max = max([x.count for x in g_count[the_slice[0]:the_slice[1]]])
        print(f"greedy count got. max count: {g_count_max} / {len(seqs)}")
        # print([x.count for x in g_count])
        print("getting min greedy stamps")
        g_count_stamps = get_stamp_min(g_count, probe_min_length, the_slice)
        # print(g_count_stamps)
        print("------")
        print("getting max probe")
        probe = get_max_probe(g_count_stamps, g_count, the_slice)
        print(f"got probe at {probe.start}: {probe.seq}")
        probes.append(probe)
        continue


if __name__ == "__main__":
    main()
else:
    print("step1 was not launched to execute. Skipping execution.")
