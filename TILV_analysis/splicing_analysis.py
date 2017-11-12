#! /usr/local/python_anaconda/bin/python3.4

from optparse import OptionParser
from file_utilities import check_filename, check_dirname
import pandas as pd
import os
import glob


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

def get_splice_junctions_information(csv, segments_refs):
    df = pd.read_csv(csv)
    segments_refs = glob.glob(segments_refs + "/*fasta")

    df = df.iloc[:, 0:5]
    df.columns = ["sample", "segment", "count", "start_intron", "end_intron"]
    df = df.sort_index(by=["start_intron", "end_intron"])
    for segment_file in segments_refs:
        segment = int(segment_file.split(".fasta")[0].split("segment")[-1])

        temp = df.loc[df["segment"] == segment]

        print(segment_file)
        sequence = open(segment_file, "r").readlines()[1]

        count_true = 0
        count_false = 0
        for index, row in temp.iterrows():
            start = row[3]
            end = row[4]
            sample = row[0]
            #print(start, end)
            count = row[2]
            start_seq = sequence[start-3:start+3]
            end_seq = sequence[end-3:end+3]
            consensus = False


            if "GT" in start_seq and "AG" in end_seq:
                consensus = True
                count_true += count

            else:
                count_false += count
            if start == 279 and end == 1460:
                print(sample, start, end, start_seq, end_seq, count, consensus)

        print("counts:", count_true, count_false)


    return(df)


def check_duplicity_of_blast_results(blast_perfix):
    blast_files = glob.glob(blast_perfix)

    for blast_file in blast_files:
        print(blast_file)
        try:
            blast_df = pd.read_csv(blast_file, sep="\t", header=None)
        except:
            continue
        blast_df.columns = ["sseqid", "qstart", "qend", "sstart", "ssend", "sstrand", "length", "btop"]
        df = pd.DataFrame(columns=["sseqid", "first_start", "first_end", "second_start", "second_end", "strand",
                                   "repeat", "repeat_size", "repeat_start", "repeat_end"])

        ids = blast_df.sseqid.unique()
        ok = 0
        not_ok = 0
        minus = 0

        for id in ids:
            strand = "plus/plus"
            repeat = "no"
            repeat_size = 0
            repeat_start = 0
            repeat_end = 0
            temp = blast_df[blast_df.sseqid == id]
            if len(temp) == 1:
                continue
            if len(temp) > 2:
                continue

            if temp.iloc[0][5] == "minus" or temp.iloc[1][5] == "minus":
                minus += 1
                strand="plus/minus"
                continue



            first_start_query = temp.iloc[0][3]
            first_end_query = temp.iloc[0][4]
            second_start_query = temp.iloc[1][3]
            second_end_query = temp.iloc[1][4]



            if (second_start_query < first_start_query):
                first_start_query = temp.iloc[1][3]
                first_end_query = temp.iloc[1][4]
                second_start_query = temp.iloc[0][3]
                second_end_query = temp.iloc[0][4]




            #print(id, first_start_query, first_end_query, second_start_query, second_end_query)



            if first_end_query <= second_start_query:
                 ok += 1
            else:
                not_ok += 1
                repeat = "yes"
                repeat_size = first_end_query - second_start_query
                repeat_start = second_start_query
                repeat_end = first_end_query



            first_start = temp.iloc[0][1]
            first_end = temp.iloc[0][2]
            second_start = temp.iloc[1][1]
            second_end = temp.iloc[1][2]



            if second_start < first_end:
                first_start = temp.iloc[1][1]
                first_end = temp.iloc[1][2]
                second_start = temp.iloc[0][1]
                second_end = temp.iloc[0][2]

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

            df = df.append({"sseqid": id, "first_start": first_start, "first_end": first_end,
                            "second_start": second_start, "second_end": second_end,
                            "strand": strand, "repeat": repeat, "repeat_size": repeat_size,
                            "repeat_start": repeat_start,
                            "repeat_end": repeat_end}, ignore_index=True)
        print(ok, not_ok, minus)




        output = blast_file.split(".blast")[0] + "_information.csv"
        df.to_csv(output)

if __name__ == "__main__":
    main()