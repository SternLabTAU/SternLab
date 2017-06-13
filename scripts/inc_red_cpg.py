#! /usr/local/python_anaconda/bin/python3.4
import sys
from Bio.Alphabet import generic_dna
from Bio.Seq import Seq
from optparse import OptionParser


def main():
    parser = OptionParser("usage: %prog [options]")
    parser.add_option("-s", "--sequence", dest="sequence", help="DNA sequence you want to increase/reduce CpG")
    parser.add_option("-i", "--output_file_increased", dest="output_file_increased",
                      help="report file of increased CpG sequence")
    parser.add_option("-r", "--output_file_reduced", dest="output_file_reduced",
                      help="report file of reduced CpG sequence")
    (options, args) = parser.parse_args()
    input_seq = options.sequence
    output_file_increase = options.output_file_increased
    output_file_reduced = options.output_file_reduced

    do_increase_cpg(input_seq, output_file_increase)
    do_reduce_cpg(input_seq, output_file_reduced)


def do_increase_cpg(input_file, out_file):
    """
    :param input_file: The Sequence to increase CpG
    :param out_file: txt file of increased CpG sequence with synonymous substitution
    :return: Increased CpG sequence with synonymous substitution
    """

    with open(input_file, 'r') as entire_file:
        seq = entire_file.read()
        seq = seq.upper()
        seq = Seq(seq, generic_dna)
        original_seq = seq
        translated = seq.translate()
        for i in range(0, len(seq) - 1, 3):
            if (seq[i + 1:i + 3] == "CT") or (seq[i + 1:i + 3] == "CC") or (seq[i + 1:i + 3] == "CA"):
                cpg = "CG"
                changed = False
                good_change = ""
                translated_new_seq = Seq(str(seq[:i + 1]) + cpg + str(seq[i + 3:])).translate()

                if translated == translated_new_seq:
                    changed = True
                    good_change = cpg
                    seq = Seq(str(seq[:i + 1]) + good_change + str(seq[i + 3:]))

    with open(out_file, 'w') as out_file1:
        print(">Original Sequence" + '\n' + original_seq, file=out_file1)
        print("No. of CpG:" + str(original_seq.count("CG")), file=out_file1)
        print(">Sequence after CpG was increased:" + '\n' + str(seq), file=out_file1)
        print("No. of CpG:" + str(seq.count("CG")), file=out_file1)
        print(">Original Protein sequence" + '\n' + (Seq(str(original_seq)).translate()), file=out_file1)
        print(">Protein sequence after CpG was increased" + '\n' + str(seq.translate()), file=out_file1)
    return seq


def do_reduce_cpg(input_file, out_file):
    """
    :param input_file: The Sequence to reduce CpG
    :param out_file: txt file of reduced CpG sequence with synonymous substitution
    :return: Reduced CpG sequence with synonymous substitution
    """

    with open(input_file, 'r') as entire_file:
        seq = entire_file.read()
        seq = seq.upper()
        seq = Seq(seq, generic_dna)
        original_seq = seq
        translated = seq.translate()
        for i in range(len(seq) - 1):
            if (seq[i:i + 2] == "CG"):
                bases = "ACGT"
                proposed_dinucs = [x + y for x in bases for y in bases]
                proposed_dinucs.remove("CG")
                changed = False
                good_change = ""

                for proposed_change in proposed_dinucs:
                    translated_new_seq = Seq(str(seq[:i]) + proposed_change + str(seq[i + 2:])).translate()
                    if (translated == translated_new_seq):
                        changed = True
                        good_change = proposed_change
                        break
                if not changed and translated[
                            i // 3] == 'R':  # special case, arginine "CGC" must be replaced into "AGA"
                    seq = Seq(str(seq[:i]) + "AGA" + str(seq[i + 3:]))

                if changed:
                    seq = Seq(str(seq[:i]) + good_change + str(seq[i + 2:]))
    with open(out_file, 'w') as out_file2:
        print(">Original Sequence" + '\n' + original_seq, file=out_file2)
        print("No. of CpG:" + str(original_seq.count("CG")), file=out_file2)
        print(">Sequence after CpG was reduced:" + '\n' + str(seq), file=out_file2)
        print("No. of CpG:" + str(seq.count("CG")), file=out_file2)
        print(">Original Protein sequence" + '\n' + (Seq(str(original_seq)).translate()), file=out_file2)
        print(">Protein sequence after CpG was reduced" + '\n' + str(seq.translate()), file=out_file2)
    return seq


def exp_cpg(seq):
    seq = seq.upper()
    c_count = seq.count('C')
    g_count = seq.count('G')
    len_seq = len(seq)
    return c_count*g_count/len_seq

if __name__ == "__main__":
    main()
