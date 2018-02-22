import pandas as pd
import subprocess
import numpy as np
import os
import argparse
from tqdm import tqdm


def main(args):

    fastq_dir = args.fastq
    freq_dir = args.freq
    type= args.type

    fastq_files = []
    freq_files = []

    sample_reads = []
    sample_coverage = []
    median_cov = []
    reads = []

    for root, dirs, files in os.walk(freq_dir):
        new_files = [os.path.join(root,f) for f in files if '.freq' in f and '_' not in f]
        freq_files.extend(new_files)

    for root, dirs, files in os.walk(fastq_dir):
        new_files = [os.path.join(root,f) for f in files if f.endswith('fastq')]
        fastq_files.extend(new_files)

    for f in tqdm(fastq_files):
        total_reads = get_num_reads(f, type)
        reads.append(total_reads)
        s = os.path.basename(f).split('.')[0]
        sample_reads.append(s)

    for f in tqdm(freq_files):
        coverage = get_median_coverage(f)
        median_cov.append(coverage)
        s = f.split('/')[-2]
        sample_coverage.append(s)


    reads = pd.DataFrame({'Total_reads':reads,'Sample_reads':sample_reads})
    coverage = pd.DataFrame({'Median_coverage': median_cov, 'Sample_coverage':sample_coverage})
    reads.to_csv(os.path.join(args.out, 'reads.csv'), index=False)
    coverage.to_csv(os.path.join(args.out, 'covergae.csv'), index=False)



def get_num_reads(fastq, type):
    """
    returns number of reads in a fastq file according counts of @M
    :param fastq: full file path to a fastq file, not zipped
    :param type: the type of sequencing machine. highly important becuse of the @X grep
    :return: number of reads
    """
    if type == 'n':
        reads = subprocess.getoutput("cat {} | grep @N | wc -l".format(fastq))

    else:
        assert( type == 'm')
        reads = subprocess.getoutput("cat {} | grep @M | wc -l".format(fastq))

    return reads

def get_median_coverage(freq):
    """
    returns median coverage
    :param freq: a full path to a freq file
    :return: median coverage
    """

    df = pd.read_csv(freq, sep='\t')
    df = df.drop_duplicates(['Pos'])
    return np.median(df['Read_count'])



if __name__ == "__main__":
    parser = argparse.ArgumentParser()
    parser.add_argument("-fq", "--fastq", type=str, help="input dir with fastq file, unzipped ", required=True)
    parser.add_argument("-fr", "--freq", type=str, help="input dir with freq files", required=True)
    parser.add_argument("-t", "--type", type=str, help="m for miseq and n for nextseq", required=True)
    parser.add_argument("-o", "--out", type=str, help="output directory to save the result", required=True)
    args = parser.parse_args()
    main(args)

