from TheSeq import TheSeq


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
