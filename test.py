from TheSymbol import TheSymbol
from TheProbeStart import TheProbeStart
from TheProbe import TheProbe
import step1_main

test_counter = 0


def print_if_not_equal(title, left, right):
    global test_counter
    test_counter += 1
    if left != right:
        print("----- TEST. " + title + " failed: ")
        print(f"----------- {left}")
        print(f"-----------  !=")
        print(f"----------- {right}")
        print("-----------------------------------------")
        test_counter -= 1
        return


if __name__ == "__main__":
    the_src = [TheSymbol("A", 10), TheSymbol("A", 11), TheSymbol("A", 12), TheSymbol("G", 13), TheSymbol("G", 9),
               TheSymbol("G", 13), TheSymbol("T", 2), TheSymbol("G", 20), TheSymbol("G", 21), TheSymbol("T", 5),
               TheSymbol("G", 9), TheSymbol("T", 4), TheSymbol("T", 5), TheSymbol("T", 7), TheSymbol("T", 9),
               TheSymbol("C", 4), TheSymbol("T", 20), TheSymbol("T", 28), TheSymbol("T", 1), TheSymbol("T", 99),
               TheSymbol("T", 9), TheSymbol("T", 9), TheSymbol("T", 9), TheSymbol("T", 9), TheSymbol("T", 9)]
    TheProbeStart.set_source(the_src)
    the_in = [1, 4]
    the_out = 11
    print_if_not_equal("min search in probe", TheProbeStart.find_min(*the_in), the_out)

    src = TheProbeStart(*the_in)
    the_out = -2
    print_if_not_equal("probe delta negative", src.get_move_right_delta(), the_out)

    the_in = [0, 3]
    the_out = 1
    src = TheProbeStart(*the_in)
    print_if_not_equal("probe delta positive", src.get_move_right_delta(), the_out)

    the_in = [0, 3]
    the_out = -1
    src = TheProbeStart(*the_in)
    print_if_not_equal("probe delta negative with step2", src.get_move_right_delta(2), the_out)

    the_in = [5, 8]
    the_out = 0
    src = TheProbeStart(*the_in)
    print_if_not_equal("probe delta zero", src.get_move_right_delta(), the_out)

    the_in = [the_src, 15]
    the_out = [2] * 4 + [1] * 7
    print_if_not_equal("stamp min 15 slots for 25 len", step1_main.get_stamp_min(*the_in), the_out)

    the_in = [the_src[0:15], 15]
    the_out = [2]
    print_if_not_equal("stamp min 15 slots for 15 len", step1_main.get_stamp_min(*the_in), the_out)

    the_in = [[2] * 4 + [1] * 7, the_src]
    the_out = "AAAGGGTGGTGTTTTC"  # note the 16 char long. It comes from main file
    the_probe = step1_main.get_max_probe(*the_in)
    print_if_not_equal("creating max probe", the_probe.seq, the_out)

    the_in = [the_probe.start, the_probe.end]
    the_out = [0, 16]
    print_if_not_equal("probe start end indexes", the_in, the_out)

    the_probe.start = 4
    the_probe.end += 4
    the_in = ["AAAGAAAGGGTGGTGTTTTCAGAGAGAGAGAGAG", the_probe]
    the_out = True
    print_if_not_equal("probe seq ideal match", step1_main.is_sequence_detected(*the_in), the_out)

    the_in = ["AAAGCAAGGGTGGTGTTTTCAGAGAGAGAGAGAG", the_probe]
    the_out = False
    print_if_not_equal("probe seq no match because of start", step1_main.is_sequence_detected(*the_in), the_out)

    the_in = ["AAAGAAAGGCTGGTGTTTTCAGAGAGAGAGAGAG", the_probe]
    the_out = True
    print_if_not_equal("probe seq match with 1 mistake", step1_main.is_sequence_detected(*the_in), the_out)

    the_in = ["AAAGAAAGGGTCCTGTTATCAGAGAGAGAGAGAG", the_probe]
    the_out = True
    print_if_not_equal("probe seq match with 3 mistakes", step1_main.is_sequence_detected(*the_in), the_out)

    the_in = ["AAAGAAAGGGTAAAATTTTCAGAGAGAGAGAGAG", the_probe]
    the_out = False
    print_if_not_equal("probe seq no match with 4 mistakes", step1_main.is_sequence_detected(*the_in), the_out)

    the_in = ["AAAGAAAGGGTGGTG-TTTCAGAGAGAGAGAGAG", the_probe]
    the_out = True
    print_if_not_equal("probe seq match with 1 mistake - gap. It depends on main variable probe_gap_is_evil",
                       step1_main.is_sequence_detected(*the_in), the_out)

    the_in = step1_main.get_stamp_min(the_src, 15, [4, 21])
    the_out = [198, 198, 198, 198, 1, 1, 1, 198, 198, 198, 198]
    print_if_not_equal("spliced g_count and g_min_stamp", the_in, the_out)

    print(f"{test_counter} tests were successful")
