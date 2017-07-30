
"""
@Author: daniellem1

"""

'''preprocess .freqs files in order to get for each genome position it's num of reads'''

import os.path
import os
import pandas as pd
import random
import matplotlib
matplotlib.use('Agg')
import matplotlib.pyplot as plt
import seaborn as sns
import argparse
import glob


def main():
#RV
    # freqs_file_path = "Z:/volume1/okushnir/Cirseq/RV/20170322_output_all_23_qscore/RVB14p2.freqs"
    # out_dir = "Z:/volume1/okushnir/Cirseq/RV/20170322_output_all_23_qscore/plots/"
    # coverage_graph(freqs_file_path, out_dir)
#CV
    freqs_file_path = r"W:\volume1\okushnir\Cirseq\CV\20170719_q30r2_edited\CVB3-p2.freqs"
    out_dir = r"W:\volume1\okushnir\Cirseq\CV\20170719_q30r2_edited\plots"
    coverage_graph(freqs_file_path, out_dir)

def coverage_graph(freqs, out_dir):
    # show a unified graph otherwise

    data = parse_reads(freqs)
    pos = data[0]
    reads = data[1]
    graph = plt.plot(pos, reads, color="DarkOrchid")

    plt.xlabel("Position In The Genome[bp]", fontsize=20)
    plt.ylabel("Number Of Reads", fontsize=20)
    plt.title("Coverage", fontsize=30)
    plt.xticks(fontsize=16)
    plt.yticks(fontsize=16)
    sns.set_style("darkgrid")
    # plt.legend()
    plt.xlim(0, 6550)
    plt.ylim(1, 100000000)
    plt.yscale("log")
    plt.tight_layout()
    plt.savefig(out_dir + "\coverage_0_106.png", dpi=680)
    plt.close('all')

    return graph


def parse_reads(freqs):

    ''' this method returns a vector of reads corresponding to genome positions.
input:
        freqs file
output:
        an integer vector containing for each position in the genome it's num of reads.
'''

    path = freqs
    df = pd.read_csv(path, sep='\t')

    # remove all duplicates from Pos except the first occurrence
    # remove all x.number duplicates
    df[["Pos"]] = df[["Pos"]].astype(int)
    df = df.drop_duplicates("Pos")

    pos = df["Pos"]  # a vector of all positions
    reads = df["Read_count"]
    med_reads = reads.median()
    print(med_reads)
    return pos, reads


if __name__ == "__main__":
    main()







