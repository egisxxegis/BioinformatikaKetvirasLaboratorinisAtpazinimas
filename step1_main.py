import os
import utils
import math

from TheSeq import TheSeq
from TheSymbolsCounter import TheSymbolsCounter
from TheSymbol import TheSymbol

# settings
source_file = "data/maffta_5_done_with_external_tool.fasta.txt"
probe_max_region_bp = 200
probe_max_mismatches = 3
probe_min_start_matches = 5
probe_min_length = 16
probe_max_length = 20
probe_target_percentage = 80
gap_symbol = "-"


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


def get_counts(seqs: [TheSeq]):
    freqs = [TheSymbolsCounter()] * len(seqs[0].seq)
    # first time create
    for i in range(len(freqs)):
        freqs[i] = TheSymbolsCounter()
    # second time fill-in
    for seq in seqs:
        for i in range(len(seq.seq)):
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


def get_greedy_count(freqs: [TheSymbolsCounter]):
    gfreqs = [TheSymbol("-", 420)] * len(freqs)
    for i in range(len(freqs)):
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


def main():
    print("step1 starting.")
    do_check_settings()
    print("settings checked.")
    do_check_source()
    print("file beginning checked.")
    print("starting analysis.")
    seqs = utils.read_fasta(source_file)

    target_sequences = math.ceil(len(seqs) / 100 * probe_target_percentage)
    print("------")
    print("length of first sequence = " + str(len(seqs[0].seq)) + " .")
    print("number of sequences = " + str(len(seqs)))
    print(f"targeted amount of sequences = {target_sequences}")
    print("------")

    print("getting counts.")
    counts = get_counts(seqs)
    print("counts got. getting greedy count")
    g_count = get_greedy_count(counts)
    print(f"greedy count got. max count: {max([x.count for x in g_count])} / {len(seqs)}")


if __name__ == "__main__":
    main()
else:
    print("step1 was not launched to execute. Skipping execution.")
