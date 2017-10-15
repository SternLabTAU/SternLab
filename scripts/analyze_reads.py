

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
from Bio.Blast import NCBIWWW as ncbi


def main():
    # for Cluster
    parser = OptionParser("usage: %prog [options]")
    parser.add_option("-t", "--tmp_dir", dest="tmp_dir", help="the path of the tmp directory")
    (options, args) = parser.parse_args()
    tmp = options.tmp_dir

    #for Local
    # tmp = "/Volumes/STERNADILABTEMP$/volume1/okushnir/Cirseq/CV/test/tmp/"

    analyze_reads(tmp)


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


def analyze_reads(tmp_cirseq_dir):
    fasta_files = glob.glob(tmp_cirseq_dir + "/*.fasta")
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
        blast_df.to_csv(blast_file + ".edges.csv", sep=',', encoding='utf-8')


if __name__ == "__main__":
    main()