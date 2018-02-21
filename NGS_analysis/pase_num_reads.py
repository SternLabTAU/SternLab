import pandas as pd
import subprocess
import numpy as np
import os
import argparse


def main(args):

    fastq_dir = args.fastq
    freq_dir = args.freq

    fastq_files = [os.path.join(fastq_dir, f) for f in os.listdir(fastq_dir) if 'fastq' in f]
    freq_files = [os.path.join(freq_dir, f) for f in os.listdir(freq_dir) if 'freq' in f]

    sample_reads = []
    sample_coverage = []
    median_cov = []
    reads = []

    for f in fastq_files:
        total_reads = get_num_reads(f)
        reads.append(total_reads)
        s = os.path.basename(f).split('.')[0]
        sample_reads.append(s)

    for f in freq_files:
        coverage = get_median_coverage()
        median_cov.append(coverage)
        s = os.path.basename(f).split('.')[0]
        sample_coverage.append(s)


    result = pd.DataFrame({'Median_coverage': median_cov, 'Sample_coverage':sample_coverage, 'Total_reads':reads,
                           'Sample_reads':sample_reads})
    result.to_csv(os.path.join(args.out, 'reads_and_coverage.csv'), index=False)






def get_num_reads(fastq):
    """
    returns number of reads in a fastq file according counts of @M
    :param fastq: full file path to a fastq file, not zipped
    :return: number of reads
    """

    reads = subprocess.getoutput("cat {} | grep @M | wc -l".format(fastq))

    return reads

def get_median_coverage(freq):
    """
    returns median coverage
    :param freq: a full path to a freq file
    :return: median coverage
    """

    df = pd.read_csv(freq, sep='\t')
    return np.median(df['Read_count'])



if __name__ == "__main__":
    parser = argparse.ArgumentParser()
    parser.add_argument("-fq", "--fastq", type=str, help="input dir with fastq file, unzipped ", required=True)
    parser.add_argument("-fr", "--freq", type=str, help="input dir with freq files", required=True)
    parser.add_argument("-o", "--out", type=str, help="output directory to save the result", required=True)
    args = parser.parse_args()
    main(args)

