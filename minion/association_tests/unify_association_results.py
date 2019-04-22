#! /usr/local/python_anaconda/bin/python3.4


import pandas as pd
import os
import argparse
from tqdm import tqdm

def association_results_to_df(args):
    '''
    This function gets a directory with numbered directories, each one containing
    the files with the association test results. It creates a new csv file with 
    the chi square statistic and the p-value for every pair of positions.
    '''
    data = []
    results_directories = [args.input_results_directory + '/' + d for d in os.listdir(args.input_results_directory) if d.isnumeric()]
    for file_dir in tqdm(results_directories):
        files = [file_dir + '/' + f for f in os.listdir(file_dir) if '.csv' not in f]
        for f in files:
            [i,j] = f.split('/')[-1].split('_')
            i, j = int(i), int(j)
            with open(f) as fi:
                g = fi.read()
                pvalue = float(g.split('\n')[1])
                chi2 = float(g.split('\n')[0])
                data.append((i,j,pvalue,chi2))
                data.append((j,i,pvalue,chi2))
    df = pd.DataFrame(data, columns=['pos1', 'pos2', 'pvalue', 'chi2'])
    df.to_csv(args.output_csv, index=False)
    return

    
if __name__ == "__main__":
    parser = argparse.ArgumentParser()
    parser.add_argument("-i", "--input_results_directory", type=str, help="path to blasts df csv", required=True)
    parser.add_argument("-o", "--output_csv", type=str, help="a path to an output directory to save results", required=True)
    args = parser.parse_args()
    if not vars(args):
        parser.print_help()
        parser.exit(1)
    association_results_to_df(args)