
class __TheSeq:
    def __init__(self, name: str, seq: str):
        self.name = name
        self.org_seq = seq
        self.seq = seq.replace("\r\n", "")


def main():
    seq_file = "data/result_searched_protein_aligned_3_var2.fasta.txt"
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
                seqs.append(__TheSeq(head, body))
                body = ""
                head = line
                continue
            body = body + line
    seqs.append(__TheSeq(head, body))
    seq_file.close()

    seq_sizes = []
    for i in range(len(seqs)):
        seq_sizes.append(len(seqs[i].seq))
    print(f"max seq len = {max(seq_sizes)}")
    threshold = int(max(seq_sizes) * 0.8)
    print(f"80% of max = {threshold}")

    seq_below = []
    for i in range(len(seqs)):
        if seq_sizes[i] < threshold:
            seq_below.append(i)

    print(f"{len(seq_below)} seqs' length is below 80% of max len")
    output_file = "data/no_trash_4.fasta.txt"
    print(f"creating new file without them. Filename: {output_file}")

    seq_file = open(output_file, "w")
    for i in range(len(seqs)):
        if i in seq_below:
            continue
        seq_file.write(seqs[i].name)
        seq_file.write(seqs[i].org_seq)
    seq_file.close()
    print(f"created new file with {len(seqs) - len(seq_below)} seqs")

    return 0


if __name__ == "__main__":
    main()

else:
    print("step0 was executed not as main. Execution did not happen")
