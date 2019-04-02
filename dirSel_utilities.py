#! /usr/local/python_anaconda/bin/python3.4

import os
import glob
from file_utilities import check_dirname, check_filename
import pandas as pd
import re
from PAML_utilities import retrive_kappa_from_paml_output


def write_params_file(param_file, input_aln, input_tree, output_res, output_log, output_tree, gtr_output=None):

    param_file = check_filename(param_file, Truefile=False)
    input_aln = check_filename(input_aln)
    input_tree = check_filename(input_tree)
    output_res = check_filename(output_res, Truefile=False)
    output_log = check_filename(output_log, Truefile=False)
    output_tree =  check_filename(output_tree, Truefile=False)

    if gtr_output == None:
        kappa = 2.0
    else:
        gtr_output = check_filename(gtr_output)
        kappa = retrive_kappa_from_paml_output(gtr_output)


    with open(input_aln, "r") as aln_handle:
        aln_data = aln_handle.readline()
        first_seq = aln_data.split(">")[1].strip()

    print(kappa)
    print(first_seq)

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
