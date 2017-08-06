#! /usr/local/python_anaconda/bin/python3.4

from optparse import OptionParser
from file_utilities import check_filename, check_dirname
import pandas as pd
import os


def main():
    parser = OptionParser("usage: %prog [options]")
    parser.add_option("-b", "--blast", dest="blast_file", help="blast result file")
    parser.add_option("-o", "--output", dest="output", help="output file path")

    (options, args) = parser.parse_args()

    blast_file = options.blast_file
    output = options.output

    blast_file = check_filename(blast_file)
    output = check_filename(output, Truefile=False)





    blast_df = pd.read_csv(blast_file, sep="\t", header=None)
    blast_df.columns = ["sseqid",  "qstart",  "qend", "sstart", "ssend", "sstrand", "length", "btop"]

    ids =blast_df.sseqid.unique()
    df = pd.DataFrame(columns=["sseqid", "first_start", "first_end", "second_start", "second_end"])

    for id in ids:
        temp = blast_df[blast_df.sseqid == id]
        if len(temp) == 1:
            continue
        if len(temp) >2:
            continue

        first_start = temp.iloc[0][1]
        first_end = temp.iloc[0][2]
        second_start = temp.iloc[1][1]
        second_end = temp.iloc[1][2]

        if second_start < first_end:
            first_start = temp.iloc[1][1]
            first_end = temp.iloc[1][2]
            second_start = temp.iloc[0][1]
            second_end =temp.iloc[0][2]



        if first_end >= second_start and first_end <= second_end:
            continue
        if first_start >= second_start and first_start <= second_end:
            continue
        if first_start == second_start or first_end == second_end:
            continue
        if second_start >= first_start and second_start <= first_end:
            continue
        if first_end > second_start:
            print(first_start, first_end, second_start, second_end)


        df = df.append({"sseqid":id, "first_start":first_start, "first_end":first_end,
                        "second_start":second_start, "second_end":second_end}, ignore_index=True)


        #print(first_start, first_end, second_start, second_end)

    df.to_csv(output, index=False)

def filter_reads_that_aligned_more_than_once(input_sam, output_sam):
    """
    aligned >1 times
    :param input_sam: input sam file path
    :param output_sam: output sam file path
    :return: nothing
    """
    input_sam = check_filename(input_sam)
    output_sam = check_filename(output_sam, Truefile=False)

    os.system("grep 'XS:i:' %s > %s" %(input_sam, output_sam))


if __name__ == "__main__":
    main()