import string
import os
import random

from TheSeq import TheSeq
from TheProbe import TheProbe
from TheSymbol import TheSymbol


def read_fasta(source) -> [TheSeq]:
    seq_file = source
    seq_file = open(seq_file, "r")

    seqs = []
    head = ""
    body = ""
    head_found = False
    while line := seq_file.readline():
        if not head_found:
            if not line[0] == ">":
                continue
            head = line
            head_found = True
            continue
        else:
            if line[0] == ">":
                seqs.append(TheSeq(head, body))
                body = ""
                head = line
                continue
            body = body + line
    seqs.append(TheSeq(head, body))
    seq_file.close()
    return seqs


def get_probes(the_file_name):
    file = open(the_file_name, "r")
    step = 1
    the_len = 0
    reference_seq = ""
    reference_name = ""
    start = 0
    end = 0
    probes = []
    while line := file.readline():
        if line[0] == "#":
            continue
        line = line.replace("\r", "").replace("\n", "")
        if step == 0:  # legacy. ignore this step
            the_target_len = int(line)
            step += 1
        elif step == 1:
            the_len = int(line)
            step += 1
        elif step == 2:
            start = int(line) - 1
            step += 1
        elif step == 3:
            end = int(line)
            step += 1
        elif step == 4:
            seq = line
            probes.append(TheProbe(start, end,  [TheSymbol("z", 999)] * end))
            probes[-1].seq = seq
            if len(probes) < the_len:
                step = 2
                continue
            step += 1
        elif step == 5:
            reference_name = line
            step += 1
        elif step == 6:
            reference_seq = line
            step += 1
        continue
    return probes, TheSeq(reference_name, reference_seq)


def is_sequence_detected(seq: [chr], probe: TheProbe, the_gap_symbol, prb_min_start_mtch, prb_no_gap, prb_max_mismtch):
    seq2 = seq[probe.start:probe.end]
    # start matches?
    # print(seq2)
    # print(probe.seq)
    if not seq2[0:prb_min_start_mtch] == probe.seq[0:prb_min_start_mtch]:
        return False
    mismatches = 0
    for i in range(prb_min_start_mtch, len(probe.seq)):
        if seq2[i] != probe.seq[i]:
            mismatches += 1
        if seq2[i] == the_gap_symbol and prb_no_gap:
            mismatches += prb_max_mismtch
    # mismatches exceeded?
    if mismatches > prb_max_mismtch:
        return False
    return True


def get_two_full_free_filenames(folder):
    if folder[-1] != "/":
        folder = folder + "/"
    start = "zzz_7"
    end = ".txt"
    f1 = ""
    for i in range(20):
        f1 = ''.join(random.choices(string.ascii_uppercase + string.digits, k=10))
        if not os.path.exists(folder + start + f1 + end):
            break
        if i+1 == 20:
            f1 = ''.join(random.choices(string.ascii_uppercase + string.digits, k=30))

    f2 = ""
    for i in range(20):
        f2 = ''.join(random.choices(string.ascii_uppercase + string.digits, k=11))
        if not os.path.exists(folder + start + f2 + end):
            break
        if i+1 == 20:
            f2 = ''.join(random.choices(string.ascii_uppercase + string.digits, k=31))

    return folder + start + f1 + end, folder + start + f2 + end
