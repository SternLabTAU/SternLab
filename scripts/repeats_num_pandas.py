#! /usr/local/python-anaconda-2.7//bin/python

import sys, itertools
sys.path.insert(0, '/sternadi/home/volume1/okushnir/scripts/')
import matplotlib
matplotlib.use('Agg') # Must be before importing matplotlib.pyplot or pylab!
import pandas as pd
import seaborn as sns
import matplotlib.pyplot as plt
import numpy as np
from optparse import OptionParser
# from file_utilities import check_dirname
import re
import glob
import os.path

sns.set_context("talk")


def main():
    # for Cluster Run
    parser = OptionParser("usage: %prog [options]")
    parser.add_option("-d", "--dir", dest="dir", help="dir of temp files of cirseq pipeline run")
    parser.add_option("-o", "--output", dest="output", help="output folder to save")
    (options, args) = parser.parse_args()
    in_dir = options.dir
    out_dir = options.output
    # in_dir = check_dirname(in_dir)
    # out_dir = check_dirname(out_dir)

    # for local run
    #CV
    # in_dir = "C:/Users/Oded/Google Drive/Studies/PhD/test/tmp"
    # out_dir = "C:/Users/Oded/Google Drive/Studies/PhD/test"

    #RV
    # in_dir = "Z:/volume1/okushnir/Cirseq/RV/20170322_output_all_23_qscore/tmp"
    # out_dir = "Z:/volume1/okushnir/Cirseq/RV/20170322_output_all_23_qscore"

    if os.path.isfile(out_dir + '/repeat_summery.npy'):
            repeats_dict = np.load(out_dir + '/repeat_summery.npy').item() #dict for read vs. repeat  , encoding='latin1'
    else:
            print("Counting the repeats number")
            repeats_dict = get_repeats_num(in_dir, out_dir) #dict for read vs. repeat
    repeat_read_plot(repeats_dict, out_dir)
    if os.path.isfile(out_dir + '/length_df.csv'):
        read_and_repeat_length = pd.DataFrame.from_csv(out_dir + '/length_df.csv', sep=',',
                                                       encoding='utf-8')  # pandas_df for read and repeat length
    else:
        print("Getting the read and repeat length")
        read_and_repeat_length = get_read_and_repeat_length(in_dir, out_dir)  # pandas_df for read and repeat length
    repeat_len_plot(read_and_repeat_length, out_dir)
    read_len_plot(read_and_repeat_length, out_dir)

"""Functions"""
def get_repeats_num(tmp_cirseq_dir, out_dir):
    files = glob.glob(tmp_cirseq_dir + "/*.fasta.blast.freqs.stats")
    repeat_summery = {}
    for file in files:
        pattern = re.compile("(\d+\t{1}\d+\n{1})", re.MULTILINE)
        text = open(file, "r").read()
        reads = pattern.findall(text)
        for r in reads:
            key = int(r.split("\t")[0])
            value = int(r.split("\t")[1].split("\n")[0])
            if key in repeat_summery:
                repeat_summery[key] += value
            else:
                repeat_summery[key] = value
    np.save(out_dir + '/repeat_summery.npy', repeat_summery)
    print("repeat_summery.npy is in your folder")
    return repeat_summery


def get_read_and_repeat_length(in_dir, out_dir):
    files = glob.glob(in_dir + "/*.fasta.blast")
    arr_names = ["sseqid", "qstart", "qend", "sstart", "send", "sstrand", "length", "btop"]
    wdf = pd.DataFrame()
    for file in files:
        data = pd.read_csv(file, sep="\t",  header=None, names=arr_names)
        grouped_and_filtered = data.groupby("sseqid").filter(lambda x: min(x["qend"]) - max(x["qstart"]) > 0).groupby("sseqid")["length"].agg([np.count_nonzero, np.sum, np.max])
    wdf = pd.DataFrame.append(wdf, grouped_and_filtered)
    wdf.to_csv(out_dir + '/length_df.csv', sep=',', encoding='utf-8')
    print("length_df.csv is in your folder")
    return wdf

"""Plots"""
def repeat_read_plot(repeat_summery, out_dir):
    x = []
    y = []
    for key, val in repeat_summery.items():
        x.append(key)
        y.append(val)
    plt.bar(x, y, color="DarkOrchid", align="center")
    ax = plt.gca()
    ax.set_xlabel("Number of Repeats", fontsize=16)
    ax.set_ylabel("Number of Reads", fontsize=16)
    ax.set_xticks(list(range(1, 11)))
    plt.title('Amount of Reads per Repeat', fontsize=22)
    ax.yaxis.get_major_formatter().set_powerlimits((0, 1))
    plt.ylim(0.0, 2.0*10**7, 0.2*10**7)
    sns.set_style("darkgrid")
    path = out_dir + "/plots/reads_vs_repeats_with_lim.png"
    plt.savefig(path, dpi=680)
    plt.close('all')


def repeat_len_plot(data, out_dir):
    sns.boxplot(x="count_nonzero", y="amax", data=data, color='DarkOrchid')
    plt.xlabel("Number of Repeats", fontsize=16)
    plt.ylabel("Repeat Length (bp)", fontsize=16)
    plt.title('Repeat Length', fontsize=22)
    sns.set_style("darkgrid")
    path = out_dir + "/plots/repeat_length.png"
    plt.savefig(path, dpi=300)
    plt.close()


def read_len_plot(data, out_dir):
    sns.boxplot(x="count_nonzero", y="sum", data=data, color='DarkOrchid')
    plt.xlabel("Number of Repeats", fontsize=16)
    plt.ylabel("Read Length (bp)", fontsize=16)
    plt.title('Read Length', fontsize=22)
    sns.set_style("darkgrid")
    path = out_dir + "/plots/read_length.png"
    plt.savefig(path, dpi=300)
    plt.close()
if __name__ == "__main__":
    main()
