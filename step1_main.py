import os
import utils
from TheSeq import TheSeq

# settings
source_file = "data/maffta_5_done_with_external_tool.fasta.txt"
probe_max_region_bp = 200
probe_max_mismatches = 3
probe_min_start_matches = 5
probe_min_length = 16
probe_max_length = 20
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


def main():
    print("step1 starting.")
    do_check_settings()
    print("settings checked.")
    do_check_source()
    print("file beginning checked.")
    print("starting analysis.")
    seqs = utils.read_fasta(source_file)
    print("length of first sequence = " + str(len(seqs[0].seq)) + " .")


if __name__ == "__main__":
    main()
else:
    print("step1 was not launched to execute. Skipping execution.")
