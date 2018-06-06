import os
import pandas as pd
import re
from scipy.stats import chisqprob
from statsmodels.sandbox.stats.multicomp import multipletests
from optparse import OptionParser

def main():
    parser = OptionParser("usage: %prog[options]")
    parser.add_option("-g", "--gene_folder", dest="GENES_FOLDER", help="The folder that contains all branch type subfolders")
    parser.add_option("-o", "--output", dest="OUTPUT", help="Output path")
    parser.add_option("-c", "--branch_choice", dest="CHOICE", help="Choose the desired clade: 13, 6, 13up or both")

    (options, args) = parser.parse_args()

    GENES_FOLDER = options.GENES_FOLDER
    OUTPUT = options.OUTPUT
    CHOICE = options.CHOICE

    create_df(GENES_FOLDER, OUTPUT, CHOICE)


def create_df(folder, output_path, branch_num):
    """
    This function receives the folder of the results and creates a data frame of each gene and its lnL value and
    its test statistic D.
    All of the result files have one line that contains lnl so I will extract it
    :param folder: The folder that contains branch type subfolders
    :param output_path: The path where the result file will be created
    :param branch_num: "6", "13", "both" or "13up" according to the desired branch to be selected a foreground"""
    data = pd.DataFrame(columns = ["Chromosome", "Gene", "null", "alternative", "D", "Pval", "Adj-Pval"])
    cnt = 0
    for f in os.listdir(folder):
        if f != "excluded":
            f1 = open(folder + "/" + f + "/" + f + "-" + branch_num + "-n-Result", "r")
            f2 = open(folder + "/" + f + "/" + f + "-" + branch_num + "-a-Result", "r")
            f1_data = f1.read()
            f2_data = f2.read()
            lnl = re.compile("lnL.*")
            lnl_n = float(lnl.findall(f1_data)[0].split(" = ")[1])
            lnl_a = float(lnl.findall(f2_data)[0].split(" = ")[1])
            chromosome = f.split("-")[0]
            gene = f.split("-")[1]
            D = 2 * (lnl_a - lnl_n)
            pval = chisqprob(D, 1)
            new_line = pd.DataFrame({"Gene": [gene],
                                     "Chromosome": [chromosome],
                                     "alternative": [lnl_a],
                                     "null": [lnl_n],
                                     "D": [D],
                                     "Pval": [pval]},
                                    columns=["Chromosome", "Gene", "null", "alternative", "D", "Pval", "Adj-Pval"])
            data = data.append(new_line)
            cnt += 1

    # adjust the p-values using FDR method
    pvals = list(data["Pval"])
    new_pvals = multipletests(pvals, method="fdr_bh")[1]
    data["Adj-Pval"] = new_pvals

    data.to_csv(output_path)
    print(str(cnt) + " genes done")

if __name__ == "__main__":
    main()