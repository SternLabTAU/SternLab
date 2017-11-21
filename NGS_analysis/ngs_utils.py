import pandas as pd
import os
import argparse
import matplotlib.pyplot as plt
import seaborn as sns




def merge_freqs(files, out, alias):
    """
    merges all frequency files
    :param files: list of full path freq files
    :param out: output folder to save the results
    :param alias: the output filename
    :return: saves one merged csv
    """

    all_files = []

    for f in files:
        df = pd.read_csv(f, sep='\t')
        base = os.path.basename(f)
        df['Sample'] = base
        all_files.append(df)

    result = pd.concat(all_files)
    result.to_csv(os.path.join(out, '{}.csv'.format(alias)), index=False)

    return result

def get_all_freq_files(input_directory):
    """
    get all frequency files in a directory
    :param input_directory: a directory containing freq files
    :return: a list of full path freq files
    """

    return [os.path.join(input_directory, f) for f in os.listdir(input_directory) if 'freq' in f]



