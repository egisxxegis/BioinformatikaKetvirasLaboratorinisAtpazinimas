from utils import *

file_name = "data/probes_6.txt"
maffted_output_file = "data/maffted_user_output_8.fasta.txt"

probe_max_mismatches = 3
probe_min_start_matches = 5
gap_symbol = "-"
probe_gap_means_evil_in_detection = False


def save_seq_to_fasta(the_file, seq_name, seq_seq):
    file = open(the_file, "w")
    if seq_name[0] != ">":
        seq_name = ">" + seq_name
    file.write(seq_name + "\n")
    file.write(seq_seq)
    file.close()


def get_maffted_seq(the_file, the_name):
    if the_name[0] != ">":
        the_name = ">" + the_name
    seqs = read_fasta(the_file)
    for seq in seqs:
        seq.name = seq.name.replace("\r", "").replace("\n", "")
        if seq.name == the_name:
            return seq
    return None


def get_seq_without_name_from_file(the_file):
    file = open(the_file, "r")
    seq = ""
    while the_line := file.readline():
        seq = seq + the_line
    return seq


if __name__ == "__main__":
    z_probes, reference_seq = get_probes(file_name)
    line = "-------------------------------"
    the_name_of_user_input = "-------The input of user-------- UNNAMED SAMPLE."
    while True:
        the_seq = input("Input sequence for check\n").replace("\r", "").replace("\n", "")
        if len(the_seq) > 10 and the_seq.lstrip(" ")[:6] == "--file":
            c_line = the_seq.strip("\n").strip("\r").strip(" ")
            input_file = the_seq.split("--file")[1].strip(" ")
            the_seq = get_seq_without_name_from_file(input_file)
            print("Got input from file.")
        the_seq = the_seq.upper()
        if len(the_seq) != len(reference_seq.seq):
            print("Your sequence doesn't match target len! " + f"({len(the_seq)} != {len(reference_seq.seq)})\n"
                  "Mafft'ing your sequence.")

            mafft_org, mafft_input = get_two_full_free_filenames("data/")
            save_seq_to_fasta(mafft_org, reference_seq.name, reference_seq.seq)
            save_seq_to_fasta(mafft_input, the_name_of_user_input, the_seq)

            cmd = f"mafft --maxiterate 1000 --localpair " \
                  f"--addfull {mafft_input} --keeplength {mafft_org} > {maffted_output_file}"
            print("Calling the command: \n" + cmd + "\n" + line)
            os.system(cmd)
            print(line + "\nCommand has finished.")

            os.remove(mafft_org)
            os.remove(mafft_input)

            maffted_seq = get_maffted_seq(maffted_output_file, the_name_of_user_input)
            the_seq = maffted_seq.seq

        # nextly
        z_results = [is_sequence_detected(the_seq,
                                          x,
                                          gap_symbol,
                                          probe_min_start_matches,
                                          probe_gap_means_evil_in_detection,
                                          probe_max_mismatches) for x in z_probes]
        if z_results.count(True) > 0:
            print("----Positive result")
        else:
            print("--Negative result")
