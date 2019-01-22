#! /usr/local/python_anaconda/bin/python3.4

import os
import glob
from file_utilities import check_dirname, check_filename
import pandas as pd
import re

HYPHY_PROGRAM = "/Users/taliakustin/Software/hyphy-2.3.13/bin/HYPHYMP"
HYPHY_BF = "/Users/taliakustin/Software/hyphy-2.3.13/res/TemplateBatchFiles/AnalyzeDiNucData_talia.bf"


def run_hyphy(aln, tree, output, model, hyphy_program=HYPHY_PROGRAM, hyphy_bf=HYPHY_BF):
    aln = check_filename(aln)
    tree = check_filename(tree)
    output = check_filename(output, Truefile=False)
    model = check_filename(model)
    os.system("%s %s %s %s %s %s" % (hyphy_program, hyphy_bf, aln, tree, output, model))



cols= ['R_TGTT', 'R_TCTT', 'R_AAAC', 'R_TCTG', 'R_CTTT',
       'R_CTGT', 'R_GAGC', 'R_GAGT', 'R_GAGG', 'R_CGTG',
       'R_CCTC', 'R_CCGC', 'R_CGCT', 'R_CGGG', 'R_GATA',
       'R_TATC', 'R_GTTT', 'R_TATG', 'R_TATT', 'R_GGTG',
       'R_GCGT', 'R_GCGG', 'R_GCTC', 'R_GGGT', 'R_CCCT',
       'R_ACAT', 'R_ACAG', 'R_ACCC', 'R_ACGC', 'R_AATA',
       'R_AAAT', 'R_AAAG', 'R_AACA', 'R_AAGA', 'R_AGAT',
       'R_CAGA', 'R_CACT', 'R_CCCG', 'R_CATA', 'R_CACG',
       'R_AGGG', 'R_AGTG', 'R_AGCG', 'R_ATGT', 'R_ATTT',
       'R_ATCT', 'R_CACC', 'AIC', 'lnL', 'rate_Normalizer',
       'filename', 'family', 'protein', 'group']



def extract_hyphy_results(files=[]):
    df = pd.DataFrame(columns=cols)
    if files==[]:
        print("no files!!")
    for f in files:
         with open(f, "r") as file:
             family = f.split("/")[-1].split("_")[0]
             protein = f.split("/")[-1].split(family + "_")[1].split(".")[0]
             group = f.split("/")[-1].split(family + "_")[1].split(".hyphy.txt")[0]
             res = {}
             res["filename"] = f
             res["protein"] = protein
             res["group"] = group
             res["family"] = family
             data = file.readlines()
             for l in data:
                     if "AIC" in l:
                             aic = float(l.split("=")[-1].strip().split("\n")[0])
                             res["AIC"] = aic
                     elif "Log Likelihood" in l:
                         lnL = float(l.split("=")[-1].strip().split("\n")[0].split(";")[0])
                         res["lnL"] = lnL
                     elif "rate_Normalizer" in l:
                         rate_Normalizer = float(l.split("=")[-1].strip().split("\n")[0])
                         res["rate_Normalizer"] = rate_Normalizer
                     elif "R_" in l:
                         name = l.split("=")[0]
                         value = float(l.split("=")[-1].strip().split("\n")[0])
                         res[name] = value
             df = df.append(res, ignore_index=True)