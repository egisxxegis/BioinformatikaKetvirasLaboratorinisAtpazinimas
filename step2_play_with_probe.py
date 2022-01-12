from utils import *

file_name = "data/probes_6.txt"

probe_max_mismatches = 3
probe_min_start_matches = 5
gap_symbol = "-"
probe_gap_means_evil_in_detection = True


if __name__ == "__main__":
    z_probes, reference_seq = get_probes(file_name)
    while True:
        the_seq = input("Input sequence for check\n").replace("\r", "").replace("\n", "").upper()
        if len(the_seq) != len(reference_seq.seq):
            print("Your sequence doesn't match target len! " + f"({len(the_seq)} != {len(reference_seq.seq)})\n"
                                                               "Results might be inaccurate.")
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
