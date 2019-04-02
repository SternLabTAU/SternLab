#! /usr/local/python_anaconda/bin/python3.4

import os
import glob
from file_utilities import check_dirname, check_filename
import pandas as pd
import re
from PAML_utilities import retrive_kappa_from_paml_output
from phyVirus.get_baltimore import get_baltimore_classifiaction

def write_params_file(param_file, input_aln, input_tree, output_res, output_log, output_tree, gtr_output=None):

    param_file = check_filename(param_file, Truefile=False)
    input_aln = check_filename(input_aln)
    input_tree = check_filename(input_tree)
    output_res = check_filename(output_res, Truefile=False)
    output_log = check_filename(output_log, Truefile=False)
    output_tree =  check_filename(output_tree, Truefile=False)

    kappa = 2.0
    if gtr_output != None:
        gtr_output = check_filename(gtr_output)
        kappa_from_gtr = retrive_kappa_from_paml_output(gtr_output)
        if kappa_from_gtr != None:
            kappa = kappa_from_gtr

    with open(input_aln, "r") as aln_handle:
        aln_data = aln_handle.readline()
        first_seq = aln_data.split(">")[1].strip()

    params = open(param_file, "w")
    params.write("# input\n")
    params.write("_inSeqFile %s\n" % input_aln)
    params.write("_inTreeFile %s\n" % input_tree)
    params.write("_inQuerySeq %s\n" % first_seq)
    params.write("_rootAt %s\n" % first_seq)
    params.write("_useQueryFreqsAtRoot 1\n")
    params.write("\n")
    params.write("# output\n")
    params.write("_outResFile %s\n" % output_res)
    params.write("_logFile %s\n" % output_log)
    params.write("_outTreeFile %s\n" % output_tree)
    params.write("\n")
    params.write("# advanced: remove # to enable advanced parameter. Parameters are described in the ParaSel manual at https://www.sternadi.com/parasel\n")
    params.write("# _modelName hky\n")
    params.write("_fixedS 1\n")
    params.write("_initS 1.0\n")
    params.write("_fixedProbS 1\n")
    params.write("_initProbS 0.01\n")
    params.write("_initKappa %f\n" % kappa)
    params.write("_fixedKappa 1\n")
    params.write("_initAlpha 0\n")
    params.write("_fixedAlpha 1\n")
    params.write("_initBeta 1\n")
    params.write("_fixedBeta 0\n")
    params.write("_initTau 1\n")
    params.write("_fixedTau 1\n")
    params.write("_bblOpt 0\n")
    params.write("_doMutationMapping 1\n")
    params.write("#_threshold 0.3\n")

    params.close()
    return params



def return_branch_type(row):
    """
    check for each node if it's internal or external
    :param row: dataframe row
    :return: internal or external
    """
    #if the node name is N## than it's an internal node
    if re.match("^N\d+", row["node_name"]) == None:
        return "external"
    return "internal"

def analyze_dirSel_mutation_map(file, ratios_output=None):
    file = check_filename(file)
    if ratios_output == None:
        ratios_output = file.split(".dirSel.results.mutation.map")[0] + ".dirSel.ratios"
    else:
        ratios_output = check_filename(ratios_output, Truefile=False)
    base = file.split("/")[-1].split(".dirSel.results.mutation.map")[0]
    family = base.split("_")[0]
    baltimore = get_baltimore_classifiaction(family)
    mapping = {0:"A", 1:"C", 2:"G", 3:"T"}
    df = pd.read_csv(file, sep="\t", index_col=False)
    if df.shape[0] == 0:
        print("%s - has not mutation mapping" % file)
        return
    df["branch"] = df.apply(return_branch_type, axis=1)
    df["count"] = 1
    summary = df.groupby(["letter1", "letter2", "branch"], as_index=False).agg({"count": "count"})
    for num,let in mapping.items():
        summary.loc[summary.letter1 == num, "letter1"] = let
        summary.loc[summary.letter2 == num, "letter2"] = let
    summary["basename"] = base
    summary["family"] = family
    summary["baltimore"] = baltimore

    ratios = pd.DataFrame()
    for branch in ["external", "internal"]:
        for nuc1 in mapping.values():
            count_from_nuc1 = sum(summary.loc[(summary.letter1 == nuc1) & (summary.branch == branch), "count"].values)
            for nuc2 in mapping.values():
                if nuc1 == nuc2:
                    continue
                mut = nuc1 + nuc2
                mut_count = list(summary.loc[(summary.letter1 == nuc1) &
                                    (summary.letter2 == nuc2) & (summary.branch == branch), "count"].values)
                if mut_count == []:
                    mut_count.append(0)
                ratio = float(mut_count[0]) /  count_from_nuc1


                ratios = ratios.append({"basename":base, "family":family, "baltimore":baltimore, "branch":branch,
                                        "substitution":mut, "ratio":ratio}, ignore_index=True)
    ratios.to_csv(ratios_output)