import pandas as pd
import os
import argparse
#from tqdm import tqdm
import numpy as np
import glob
from sklearn.metrics import auc

def get_accuracy(filepath, threshold=0.8):
    motif = 'P4__G'

    df = pd.read_csv(filepath)

    truePositive = df[(df['name'] == motif) & (df['gammaMean'] >= threshold)].shape[0]
    falsePositive = df[(df['name'] != motif) & (df['gammaMean'] >= threshold)].shape[0]
    trueNegative = df[(df['name'] != motif) & (df['gammaMean'] < threshold)].shape[0]
    falseNegative = df[(df['name'] == motif) & (df['gammaMean'] < threshold)].shape[0]

    return truePositive, falsePositive, falseNegative, trueNegative


def folderAnalyze(folder):
    os.chdir(folder)

    x = [0.0]
    y = [0.0]
    files = glob.glob("*CT_bestChainsSummary.csv")

    for i in np.linspace(start=1, stop=0):
        results = []
        for f in files:
            results.append(get_accuracy(os.path.join(folder, f), threshold=i))
        tp = [r[0] for r in results]
        fp = [r[1] for r in results]
        fn = [r[2] for r in results]
        tn = [r[3] for r in results]

        x.append(sum(fp) / (sum(tn) + sum(fp)))
        y.append(sum(tp) / (sum(tp) + sum(fn)))

    curr_auc = auc(x, y)
    beta = int(folder.split('/')[-3].split('_')[-1])
    noised = folder.split('/')[-4]
    non_syn = int(folder.split('/')[-2].split('_')[-1])
    df = pd.DataFrame({'beta': beta, 'noised': noised, 'non_syn': non_syn, 'auc':curr_auc}, index=[0])

    return df



def main(args):

    all_result_summarys = []

    for root, dirs, files in os.walk(args.folder):
        if "resultSummary" in dirs:
            all_result_summarys.append(os.path.join(root, [f for f in dirs if 'resultSummary' in f][0]))

    dfs = []
    for f in all_result_summarys:
        try:
            df = folderAnalyze(f)
            dfs.append(df)
        except:
            print(f)

    all_results = pd.concat(dfs)
    all_results.to_csv(args.out, index=False)




if __name__ == "__main__":
    parser = argparse.ArgumentParser()
    parser.add_argument("-f", "--folder", type=str, help="the folder to analyze")
    parser.add_argument("-o", "--out", type=str, help="output folder")

    args = parser.parse_args()

    main(args)