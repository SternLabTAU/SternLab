#! /usr/local/python_anaconda/bin/python3.4

from optparse import OptionParser
import glob
from file_utilities import check_filename, check_dirname
from freqs_utilities import merge_freqs_files, filter_freqs_for_indel_analysis, add_mutation_to_freq_file, mutation_accomulation_analysis
from pandas_utilities import merge_dfs
import pandas as pd


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
        if first_end > second_start:
            print(first_start, first_end, second_start, second_end)


        df = df.append({"sseqid":id, "first_start":first_start, "first_end":first_end,
                        "second_start":second_start, "second_end":second_end}, ignore_index=True)


        #print(first_start, first_end, second_start, second_end)

    df.to_csv(output, index=False)



if __name__ == "__main__":
    main()