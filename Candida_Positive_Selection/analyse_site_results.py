from os import path
import os
import pandas as pd
import re
from scipy.stats import chisqprob
from statsmodels.sandbox.stats.multicomp import multipletests
from optparse import OptionParser


def main():
    parser = OptionParser("usage: %prog[options]")
    parser.add_option("-i", "--gene_folder", dest="GENES_FOLDER", help="The folder that contains all candida genes subfolders")
    parser.add_option("-o", "--output", dest="OUTPUT", help="Output path")
    parser.add_option("-j", "--PAML_output", dest="PAML_OUTPUT", help="The folder that contains all PAML output files")

    (options, args) = parser.parse_args()

    GENES_FOLDER = options.GENES_FOLDER
    OUTPUT = options.OUTPUT
    PAML_OUTPUT = options.PAML_OUTPUT

    check_output_fies(PAML_OUTPUT)
    create_df(GENES_FOLDER, OUTPUT)


def create_df(folder, output_path):
    """
    This function receives the folder of the results and creates a data frame of each gene and its lnL value and
    its the test statistic D.
    All of the result files have one line that contains lnl so I will extract it
    :param folder: The folder that contains all of the gene subfolders
    :param output_path: Path to result file"""
    data = pd.DataFrame(columns = ["Chromosome", "Gene", "m8", "m8a", "D", "Pval", "Adj-Pval"])
    cnt = 0
    for f in os.listdir(folder):
        if f != "excluded":
            f1 = open(folder + "/" + f + "/" + f + "-M8-Result", "r")
            f2 = open(folder + "/" + f + "/" + f + "-M8a-Result", "r")
            f1_data = f1.read()
            f2_data = f2.read()
            lnl = re.compile("lnL.*")
            lnl_M8 = float(lnl.findall(f1_data)[0].split(" = ")[1])
            lnl_M8a = float(lnl.findall(f2_data)[0].split(" = ")[1])
            chromosome = f.split("-")[0]
            gene = f.split("-")[1]
            D = 2 * (lnl_M8 - lnl_M8a)
            pval = chisqprob(D, 1)
            new_line = pd.DataFrame({"Gene": [gene],
                                     "Chromosome": [chromosome],
                                     "m8": [lnl_M8],
                                     "m8a": [lnl_M8a],
                                     "D": [D],
                                     "Pval": [pval]},
                                    columns=["Chromosome", "Gene", "m8", "m8a", "D", "Pval", "Adj-Pval"])
            data = data.append(new_line)
            cnt +=1

    # adjust the p-values using FDR method
    pvals = list(data["Pval"])
    new_pvals = multipletests(pvals, method="fdr_bh")[1]
    data["Adj-Pval"] = new_pvals

    data.to_csv(output_path)


def check_output_fies(folder):
    """Goes over all PAML output files and check they end well and have no ??? error in them"""
    cnt = 0
    files = []
    for filename in os.listdir(folder):
        if "cml" in filename:
            cnt += 1
            file = open(folder + filename, "r")
            text = file.readlines()
            last = text[-1]
            end = 0

            for line in text:
                if "I give up" in line:
                    files += filename
                if "Time used:" in line:
                    end += 1
                if "???" in line:
                    print("Error ??? in: " + filename)
                    files += filename
            if end != 2:
                print("Error 'not done' in" + filename)
                files += filename

    print("checked " + str(cnt) + " files")
    return files


def NEB_BEB(result_file):
    """
    Given a PAML result file, checks if the NEB and BEB analysis brought up the same locations
    :param result_file: The path to the result file (M8)
    :return: True or False
    """
    BEB = set()
    for line in reversed(list(open(result_file))):
        if line == "\n":
            continue
        elif line == "\tProb(w>1)  mean w\n":
            break
        else:
            line = line.split(" ")
            start = 0
            while line[start] == "":
                start += 1
                BEB.add(line[start])
    NEB = None
    for line in list(open(result_file)):
        if line == "\tProb(w>1)  mean w\n":
            NEB = set()
        elif NEB != None and "lnL" in line:
            break
        elif NEB != None and line != "\n":
            line = line.split(" ")
            start = 0
            while line[start] == "":
                start += 1
                NEB.add(line[start])
    return NEB,BEB


def compare_NEB_to_BEB(folder):
    """
    Will iterate over all of the gene folders and check if their BEB and NEB results are equal
    :param folder: The folder that contains all of the gene subfolders
    """
    cnt =0
    for f in os.listdir(folder):
        NEB, BEB = NEB_BEB(folder + f + "/" + f + "-M8-Result")
        if NEB != BEB:
            print(f + "\n" + str(NEB.issubset(BEB)) + "\n")
        else:
            cnt +=1
    print(str(cnt))

if __name__ == "__main__":
    main()