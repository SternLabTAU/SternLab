#! /usr/local/python_anaconda/bin/python3.4

import pandas as pd
from tqdm import tqdm
import argparse


def window_mean(pos1, df, window_size):
    a = df[df.pos1.isin(range(int(pos1) - int(window_size/2), int(pos1) + (int(window_size/2) + 1)))]
    window_count = len(a)
    return  a.chi2.mean(), a.chi2.std(), window_count

if __name__ == "__main__":
    parser = argparse.ArgumentParser()
    parser.add_argument('-i', '--input_path', type=str, help='path to an input csv of chi scores, created by unify_association_results.py', required=True)
    parser.add_argument("-o", "--output_path", type=str, help="path to an output csv to save results to", required=True)
    parser.add_argument('-w', '--window_size', type=int, help="sliding window size for z-test", required=True)
    args = parser.parse_args()
    if not vars(args):
        parser.print_help()
        parser.exit(1)

    df = pd.read_csv(args.input_path)
    df = df.drop_duplicates()
    df = df[(df.pos1 - df.pos2).abs() > 15]  
    
    means_and_stds = df[['pos1']].drop_duplicates()
    tqdm.pandas(desc="z score progress")
    means_and_stds['pos1_mean'], means_and_stds['pos1_std'], means_and_stds['window_count'] = zip(*means_and_stds.progress_apply(lambda x: window_mean(x.pos1, df, args.window_size), axis=1))

    merged = pd.merge(df, means_and_stds, on='pos1', how='inner')
    merged['zscore'] = (merged.chi2 - merged.pos1_mean) / merged.pos1_std
    merged.to_csv(args.output_path, index=False)