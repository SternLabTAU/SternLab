import pandas as pd
import numpy as np
import os
from pbs_runners import pipeline_runner
import glob
from MS2_analysis import add_mutation_type

'''
This file includes all the analysis for the FITS paper major revision.
'''


def main(args):
    pass



def run_mahoney_pipeline():

    reference = r'/sternadi/home/volume1/shared/data/ref_genomes/Mahoney.WT1.fasta'
    passages = [3,4,5,6,7,8]

    for p in passages:

        input_dir = os.path.join('/sternadi/datasets/volume1/WT1', 'p{}'.format(p))
        output_dir = os.path.join(input_dir, 'updated_pipeline')

        pipeline_runner(input_dir, output_dir, reference, NGS_or_Cirseq=1, gaps='N')


def merge_all_freqs(out=None):
    """
    merge all mahoney freq files
    :param out: a path to save the new data frame
    :return: a data frame containing all frequency files with additional passage column
    """
    all_freqs = []

    freqs = glob.glob('/sternadi/datasets/volume1/WT1/p*/updated_pipeline/*.freqs')

    # filter for each frequency file insertions and add passage column
    for f in freqs:
        df = pd.read_csv(f, sep='\t')
        df['passage'] = int(f.split('WT1/p')[-1][0])
        df = df[df['Pos'].astype(int) == df['Pos']]
        df['Pos'] = df['Pos'].astype(int)
        df['Mutation'] = df['Ref'] + df['Base']

        all_freqs.append(df)

    conct = pd.concat(all_freqs)

    if out != None:
        conct.to_csv(out, index=False)

    return conct

def add_mutations_information(df):
    """
    add mutation type to the data
    :param df: frequency file as a csv
    :return: the data frame with all the mutation types
    """

    coding = list(range(743, 7370))
    offset = 743
    default = 'non-coding'

    with open(r'/sternadi/home/volume1/shared/data/ref_genomes/Mahoney.WT1.fasta', 'r') as o:
        reference = o.read().replace('\n', '').split('genome')[-1]

    df['Type'] = df.apply(lambda row: add_mutation_type.get_mutation_type(row['Pos'], row['Base'], reference, offset)\
        if row['Pos'] in coding else default, axis=1)

    return df


def get_interquartile_range(df):
    """
    This method filter df to include only mutations who's frequencies are in the interquartile range (25%-75%)
    :param df: a mutation data frame
    :return: data frame filtered as described
    """

    passages = set(df.Time.values)
    all_filtered = []
    cnt = 0

    for passage in passages:
        lower_boundary = df['Freq'][df['passage'] == passage].quantile(0.25)
        upper_boundary = df['Freq'][df['passage'] == passage].quantile(0.75)

        filtered = df[(df['passage'] == passage) & (df['Freq'] >= lower_boundary) & (df['Freq'] <= upper_boundary)]
        all_filtered.append(filtered)
        cnt += filtered.shape[0]

    result = pd.concat(all_filtered)
    if cnt != result.shape[0]:
        raise Exception('Error: Dimensions of concatenated data frames do not fit\n')

    return result


def filter_input_file(df, syn=True, interquartile=True, SHAPE=True, out=None):
    """
    filter frequency file by filtering schemes
    :param df: a data frame with sequencing data.
    :param syn: indicator if to filter synonymous mutations only. default true
    :param interquartile: indicator if to filter interquartile range frequencies mutations only. default true
    :param SHAPE: TODO
    :return:
    """
    


