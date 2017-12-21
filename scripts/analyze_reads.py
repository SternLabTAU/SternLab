

"""
@Author: odedkushnir

"""

import pandas as pd
import glob
from optparse import OptionParser
from Bio.Blast import NCBIXML
from pbs_runners import blast_runner
from pbs_jobs import check_pbs
import matplotlib as plt
from nt_freq_counter import *
import subprocess
import os
import re


def main():
    # for Cluster
    # parser = OptionParser("usage: %prog [options]")
    # parser.add_option("-t", "--tmp_dir", dest="tmp_dir", help="the path of the tmp directory")
    # (options, args) = parser.parse_args()
    # tmp = options.tmp_dir
    # print(title)
    input_dir = "/sternadi/nobackup/volume1/okushnir/HeLaRVB14/pipeline/tmp/"
    out_file = "/sternadi/nobackup/volume1/okushnir/HeLaRVB14/pipeline/"

    #for Local
    # tmp = "/Users/odedkushnir/Google Drive/Studies/PhD/Python Scripts/test/tmp/"
    # tmp = "/Volumes/STERNADILABTEMP$/volume1/okushnir/Cirseq/CV/test/tmp/"

    # blast_df = analyze_reads(tmp)

    # seq = "GGGTGTGTGTGGCTTGAGGGTGTGTGTGGCTTGAGGGTGTGTGTGGCTTGAGGGTGTGTGTGGCTTGAGGGTGTGTGTGGCTTGAGGGTGTGTGTGGCTTGAGGG" \
    #       "TGTGTGTGGCTTGAGGGTGTGTGTGGCTTGAGGGTGTGTGTGGCTGGAGGGTGTGTGTGGCTTGTGGGTGTGTGTGACTTGAGTGTGTGTGTGGGTTGGGGGTGT" \
    #       "GTGTGGGGGTAGGTTGTGTGTGGGTTGTGGGTGTGTGGGGGNT" # NO Hits

    # seq = "GGGTAACAGAAGTGCTTGGACTACCAACTAGCTCAATAGACTCTTCGCACCATGTCTGTATTAGAGCGTCCCATGGGTTTCCCCATGGGCAGGCCGCCAACGCAGCC" \
    #       "ACCGCCACGGTCGCCCGTGGGGAATGCGGTGACTCATCGACCTGATCTACACTGGTTTTTCGAAGTAGTTGGCCGGATAACGAACGCTTTCTCCTTCAACCGCGTGA" \
    #       "GCAGTCTATTGATACTCAGTCCGGGGTAACAGAAGTGNG" # Human coxsackievirus B3

    # title = blast_seq(seq)
    # print(title)
    # input_dir = "/Volumes/STERNADILABTEMP$/volume1/okushnir/RCseq/CV/20171029_q23r2_evalue-1_newpipeline/tmp/"
    # out_file = "/Volumes/STERNADILABTEMP$/volume1/okushnir/RCseq/CV/20171029_q23r2_evalue-1_newpipeline/table.csv"

    make_table(input_dir, out_file, [0])

def extract_location(name, start, end):
    return(name[start:end])


def parse_fasta(file_name):
    sequences = {}
    with open(file_name, 'r') as f: # open fasta file for reading
        for line in f:# loop over file lines
            if line.startswith('>'):# if header line
                seq_id = line.split('>')[1].strip()
            else: # if sequence line
                seq = line.strip()
                sequences[seq_id] = seq
    return sequences


def blast_seq(seq):
    print("Blasting...")
    # local blast
    with open("my_fas.fas", "w") as my_fasta:
        my_fasta.write(">new seq\n"+seq)
    job_id = blast_runner("my_fas.fas", outfile="my_blast.xml", hitlist_size=1) #pbs job id

    # try:
    #     process = subprocess.check_output("qstat | grep " + str(job_id), shell=True)
    #     while process != "":
    #         process = subprocess.check_output("qstat | grep " + str(job_id), shell=True)
    #         sleep(0.05)
    # except (subprocess.CalledProcessError):
    #     process = ""
    # if process == "":
    #     print("Blasted!")

    status = check_pbs(job_id)
    print(status)
    if status == "Done!":
        xml_file = open("my_blast.xml", "r")
        blast_record = NCBIXML.read(xml_file)
        xml_file.close()
        try:
            for alignment in blast_record.alignments:
                for hsp in alignment.hsps:
                    title = (str(alignment.title).split("|")[4])
                    return title
        except (RuntimeError, TypeError, NameError, ValueError):
            return None

    # www blast
    # blast_handle = NCBIWWW.qblast("blastn", "nt", seq, hitlist_size=1)
    # print(blast_handle)
    # blast_record = NCBIXML.read(blast_handle)
    # print(blast_record)
    # for alignment in blast_record.alignments:
    #     for hsp in alignment.hsps:
    #         title = (str(alignment.title).split("|")[4])
    #         print(title)
            # return title


def sum_total_edges(df, column):
    df[column+"_len"] = df.apply(lambda x: len(x[column]), axis=1)
    df[column+"_len_sum"] = df[column+"_len"].sum()
    sum_of_sum = df[column+"_len_sum"][0]
    return df, sum_of_sum


def analyze_reads(tmp_cirseq_dir, filter_by = 0):
    fasta_files = glob.glob(tmp_cirseq_dir + "*.fasta")
    records = {}
    sum_of_sum5 = 0
    sum_of_sum3 = 0
    for file in fasta_files:
        records = parse_fasta(file)
        fasta_df = pd.DataFrame.from_dict(records, orient='index')
        fasta_df.index.name = 'id'
        fasta_df.columns = ['seq']
        blast_file = file + ".blast"
        blast_arr_names = ["sseqid", "qstart", "qend", "sstart1", "send1", "sstrand", "length", "btop", "sstart2", "send2"]
        blast_df = pd.DataFrame()

        data = pd.read_csv(blast_file, sep="\t", header=None, names=blast_arr_names)
        data["sstart2"] = data["sstart1"]
        data["send2"] = data["send1"]
        grouped_df = data.groupby("sseqid").agg({'sstart1': 'min', 'sstart2': 'max', 'send1': 'min', 'send2': 'max'})
        grouped_df['sstart'] = grouped_df.min(axis=1)
        grouped_df['send'] = grouped_df.max(axis=1)
        blast_df = pd.DataFrame.append(blast_df, grouped_df)
        blast_df = blast_df.join(fasta_df)
        blast_df['edge5'] = blast_df.apply(lambda x: extract_location(x["seq"], 0, x["sstart"]), axis=1)
        blast_df['edge3'] = blast_df.apply(lambda x: extract_location(x["seq"], x["send"], -1), axis=1)
        blast_df['edge5_100'] = blast_df['edge5'].apply(lambda x: len(x) > filter_by)
        blast_df['edge3_100'] = blast_df['edge3'].apply(lambda x: len(x) > filter_by)
        blast_df = blast_df[blast_df.edge3_100 != False]
        blast_df = blast_df[blast_df.edge5_100 != False]
        del blast_df['edge3_100']
        del blast_df['edge5_100']
        blast_df = sum_total_edges(blast_df, "edge5")[0]
        blast_df = sum_total_edges(blast_df, "edge3")[0]
        sum_of_sum5 += sum_total_edges(blast_df, "edge5")[1]
        sum_of_sum3 += sum_total_edges(blast_df, "edge3")[1]

        #blast
        # blast_df["blast5"] = blast_df.apply(lambda row: blast_seq(row["edge5"]), axis=1)
        # blast_df["blast3"] = blast_df.apply(lambda row: blast_seq(row["edge3"]), axis=1)

        blast_df.to_csv(blast_file + ".edges.csv", sep=',', encoding='utf-8')
        # plt.hist(blast_df["edge3"], bins=50)
        # plt.hist(blast_df["edge5"], bins=50)
        # plt.savefig(tmp_cirseq_dir + 'plot.png')
    with open(tmp_cirseq_dir+"edges_sum.txt", "w") as edges_len_sum:
        edges_len_sum.write("sum_of_sum5:%i \n" % sum_of_sum5)
        edges_len_sum.write("sum_of_sum3:%i \n" % sum_of_sum3)
    return sum_of_sum5, sum_of_sum3, blast_df


def tmp_fasta_to_df(tmp_cirseq_dir):
    fasta_files = glob.glob(tmp_cirseq_dir + "*.fasta")
    records = {}
    for file in fasta_files:
        records = parse_fasta(file)
        df = pd.DataFrame.from_dict(records, orient='index')
        df.index.name = 'id'
        df.columns = ['seq']
        return df

def make_table(input_dir, out_file, idx):
    """
    :param input_dir: tmp dir of CirSeq pipeline
    :param out_file:
    :return:
    """
    total_reads = int(subprocess.check_output("cat {}toFastaAndSplit.*| awk -F 'is' '{{print $2}}'| awk -F ',' '{{sum+=$1}}END{{print sum}}'".format(input_dir), shell=True))
    parts_num = int(subprocess.check_output("cat {}list_parts_fasta_files.txt | wc -l".format(input_dir), shell=True))
    mapped_ref = int(subprocess.check_output("cat {}*blast | awk '{{print $1}}'|sort | uniq | wc -l".format(input_dir),
                                         shell=True))
    not_mapped_ref = int(total_reads) - mapped_ref
    percent_mapped_ref = round(mapped_ref/int(total_reads)*100, 2)
    mapped_once = int(subprocess.check_output("grep -P '^1\t' {}*stats -h | awk '{{sum+=$2}}END{{print sum}}'"
                                          .format(input_dir), shell=True))
    supposed_freqs = mapped_ref - mapped_once
    actualy_contribute = int(subprocess.check_output("grep 'reads contributing to frequency counts' -h {}*stats | awk "
                                                     "'{{sum+=$1}}END{{print sum}}'".format(input_dir), shell=True))
    double_mapping = int(subprocess.check_output("grep 'mapped bases removed due to double mapping' {}*stats | awk -F "
                                                 "':' '{{sum+=$2}}END{{print sum}}'".format(input_dir), shell=True))
    bases_called = int(subprocess.check_output("grep 'num bases called' {}*stats | awk -F = "
                                               "'{{sum+=$2}}END{{print sum}}'".format(input_dir), shell=True))
    bases_used = nt_freq_counter(input_dir)[0]
    total_bases = nt_freq_counter(input_dir)[1]
    percent_used_bases = (bases_used/total_bases*100)
    percent_used_bases = round(percent_used_bases, 2)
    nuc_edge5 = analyze_reads(input_dir)[0]

    nuc_edge3 = analyze_reads(input_dir)[1]


    data = pd.DataFrame({'Total Reads': total_reads, 'Number of parts': parts_num, 'Mapped to reference': mapped_ref,
                         "didnt map to reference": not_mapped_ref, "% mapped to reference": percent_mapped_ref,
                         "Reads that were mapped once": mapped_once, "Reads that are supposed to contribute to freqs":
                             supposed_freqs, "Reads that actualy contribute to freqs": actualy_contribute,
                         "dbl bases removed": double_mapping, "Number of bases called": bases_called,
                         "Number of bases used": bases_used, "Total number of bases": total_bases,
                         "% used number of bases": percent_used_bases, "No of bases didnst map to ref in 5'": nuc_edge5,
                         "No of bases didnt map to ref in 3'": nuc_edge3}, index=idx)
    print(data)
    data.to_csv(out_file + "cirseq_stas.csv", sep=',', encoding='utf-8')

if __name__ == "__main__":
    main()