

"""
@Author: odedkushnir

"""

import pandas as pd
import numpy as np
import re
import glob
import Bio.Seq as Seq
from Bio import SeqIO
from optparse import OptionParser
from Bio.Blast import NCBIWWW
from Bio.Blast import NCBIXML
from pbs_runners import blast_runner
import matplotlib as plt


def main():
    # for Cluster
    parser = OptionParser("usage: %prog [options]")
    parser.add_option("-t", "--tmp_dir", dest="tmp_dir", help="the path of the tmp directory")
    (options, args) = parser.parse_args()
    tmp = options.tmp_dir

    #for Local
    # tmp = "/Users/odedkushnir/Google Drive/Studies/PhD/Python Scripts/test/tmp/"

    blast_df = analyze_reads(tmp)
    print(blast_df)

    # seq = "GGGTAACAGAAGTGCTTGGACTACCAACTAGCTCAATAGACTCTTCGCACCATGTCTGTATTAGAGCGTCCCATGGGTTTCCCCATGGGCAGGCCGC" \
    #            "CAACGCAGCCACCGCCACGGTCGCCCGTGGGGAATGCGGTGACTCATCGACCTGATCTACACTGGTTTTTCGAAGTAGTTGGCCGGATAACGAACGCT"\
    #            "TTCTCCTTCAACCGCGTGAGCAGTCTATTGATACTCAGTCCGGGGTAACAGAAGTGNG"
    # title = blast_seq(seq)
    # print(title)


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
    # blast_handle = blast_runner(seq, hitlist_size=1)
    blast_handle = NCBIWWW.qblast("blastn", "nt", seq, hitlist_size=1)
    # with open("my_blast.xml", "w") as out_handle:
    #     out_handle.write(blast_handle.read())
    blast_record = NCBIXML.read(blast_handle)
    for alignment in blast_record.alignments:
        for hsp in alignment.hsps:
            title = (str(alignment.title).split("|")[4])
            return title


def analyze_reads(tmp_cirseq_dir):
    fasta_files = glob.glob(tmp_cirseq_dir + "*.fasta")
    records = {}
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
        blast_df['edge5_100'] = blast_df['edge5'].apply(lambda x: len(x) > 100)
        blast_df['edge3_100'] = blast_df['edge3'].apply(lambda x: len(x) > 100)
        blast_df = blast_df[blast_df.edge3_100 != False]
        blast_df = blast_df[blast_df.edge5_100 != False]
        del blast_df['edge3_100']
        del blast_df['edge5_100']
        blast_df["blast5"] = blast_df.apply(lambda row: blast_seq(row["edge5"]), axis=1)
        blast_df["blast3"] = blast_df.apply(lambda row: blast_seq(row["edge3"]), axis=1)
        blast_df.to_csv(blast_file + ".edges.csv", sep=',', encoding='utf-8')
        plt.hist(blast_df["edge3"], bins=50)
        plt.hist(blast_df["edge5"], bins=50)
        plt.savefig(tmp_cirseq_dir + 'plot.png')
        return blast_df

if __name__ == "__main__":
    main()