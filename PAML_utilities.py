#! /usr/local/python_anaconda/bin/python3.4

import os
import glob
from file_utilities import check_dirname, check_filename
import pandas as pd
import re


def write_ctl_file(ctl, seq, tree, out, model, fix_alpha="1", alpha="0"):
    """
    writes a ctl file for baseml PAML
    :param ctl: ctl file output path
    :param seq: sequence file path (phylip)
    :param tree: tree file path
    :param out: output file path
    :param model: model to run
    :param fix_alpha: fix_alpha value
    :param alpha: alpha value
    :return: ctl file path
    """
    ctl_file = open(ctl, "w")

    noisy = "3"
    verbose = "1"
    runmode = "0"
    Mgene = "0"
    ndata = "1"  # was 100
    clock = "0"  # was 1
    fix_kappa = "0"
    kappa = "2.5"  # was 5
    Malpha = "0"
    ncatG = "5"
    fix_rho = "1"
    rho = "0"
    nparK = "0"
    nhomo = "0"
    getSE = "0"
    RateAncestor = "0"
    Small_Diff = "1e-6"
    cleandata = "1"
    icode = "0"
    fix_blength = "-1"
    method = "0"

    ctl_file.write("\tseqfile = %s\n" % seq)
    ctl_file.write("\ttreefile = %s\n" % tree)
    ctl_file.write("\toutfile = %s\n\n" % out)
    ctl_file.write("\tnoisy = %s\n" % noisy)
    ctl_file.write("\tverbose = %s\n" % verbose)
    ctl_file.write("\trunmode = %s\n" % runmode)
    ctl_file.write("\tmodel = %s\n" % model)
    ctl_file.write("\tMgene = %s\n" % Mgene)
    ctl_file.write("*\tndata = %s\n" % ndata)
    ctl_file.write("\tclock = %s\n" % clock)
    ctl_file.write("\tfix_kappa = %s\n" % fix_kappa)
    ctl_file.write("\tkappa = %s\n" % kappa)
    ctl_file.write("\tfix_alpha = %s\n" % fix_alpha)
    ctl_file.write("\talpha = %s\n" % alpha)
    ctl_file.write("\tMalpha = %s\n" % Malpha)
    ctl_file.write("\tncatG = %s\n" % ncatG)
    ctl_file.write("\tfix_rho = %s\n" % fix_rho)
    ctl_file.write("\trho = %s\n" % rho)
    ctl_file.write("\tnparK = %s\n" % nparK)
    ctl_file.write("\tnhomo = %s\n" % nhomo)
    ctl_file.write("\tgetSE = %s\n" % getSE)
    ctl_file.write("\tRateAncestor = %s\n" % RateAncestor)
    ctl_file.write("\tSmall_Diff = %s\n" % Small_Diff)
    ctl_file.write("*\tcleandata = %s\n" % cleandata)
    ctl_file.write("*\ticode = %s\n" % icode)
    ctl_file.write("*\tfix_blength = %s\n" % fix_blength)
    ctl_file.write("\tmethod = %s\n" % method)

    ctl_file.close()
    return ctl


def mlbs_to_df(output, mlbs = [], dirname = None):
    """
    analyzes mlb file to dataframe - extracts lnL, base frequencies and substiution matrics
    :param output: output csv file path
    :param mlbs: list of mlb files
    :param dirname: dirname that has mlb files
    :return: output file path
    """
    if mlbs == [] and dirname == None:
        raise Exception("you need to provide mlb or dirname that contains mlbs")
    if mlbs != [] and dirname != None:
        raise Exception("you need to provide only one - mlb or dirname")

    if dirname != None:
        dirname = check_dirname(dirname)
        mlbs = glob.glob(dirname + "/*.mlb")
    if mlbs != []:
        mlbs = [check_filename(m) for m in mlbs]

    output = check_filename(output, Truefile=False)

    df = pd.DataFrame(columns = ["mlb_file_name",  "basename", "family", "protein", "group", "model", "lnL",
                                 "freq_T", "freq_C", "freq_A", "freq_G",
                                 "TC", "TA", "TG", "CT", "CA", "CG", "AT",
                                 "AC", "AG", "GT", "GC", "GA"])
                        
    lnL_1 = re.compile("lnL.*")
    lnL_2 = re.compile("\-?\d*\.\d*")
    base_1 = re.compile("Base frequencies.*")
    base_2 = re.compile("0.\d+")
    rate_1 = re.compile("Rate matrix Q.*\n.*\n.*\n.*\n.*", re.IGNORECASE)
    rate_2 = re.compile("\d+.\d+")
    for mlb_file_name in mlbs:
        basename = mlb_file_name.split("/")[-1].split(".mlb")[0]
        print(mlb_file_name)
        family = mlb_file_name.split("/")[-1].split("_")[0]
        protein = mlb_file_name.split("/")[-1].split(family+"_")[1].split(".")[0]
        filename = mlb_file_name.split("/")[-1].split(family+"_")[-1].split(".mlb")[0]
        #if "_gtr" in filename or "_unrest" in filename:
        #    filename = filename.split("_gtr")[0]
        #    filename = filename.split("_unrest")[0]
        model = mlb_file_name.split(".mlb")[-1]
        

        mlb = open(mlb_file_name, "r").read()
        L = lnL_1.findall(mlb)
        if len(L) != 1:
            L = None
        elif "nan" in L[0]:
            L = None
        else:
            L = float(lnL_2.findall(L[0])[0])
        
        B = base_1.findall(mlb)
        if len(B) != 1:
            freq_T = None; freq_C = None
            freq_A = None; freq_G = None
        elif "nan" in B[0]:
            freq_T = None; freq_C = None
            freq_A = None; freq_G = None
        else:
            B =  base_2.findall(B[0])
            freq_T = float(B[0])
            freq_C = float(B[1])
            freq_A = float(B[2])
            freq_G = float(B[3])
            
        R = rate_1.findall(mlb)
        if len(R) != 1:
            TC = None; TA = None; TG = None;
            CT = None; CA = None; CG = None;
            AT = None; AC = None; AG = None;
            GT = None; GC = None; GA = None
        elif len(R) >= 1 and "nan" in R[0]:
            TC = None; TA = None; TG = None;
            CT = None; CA = None; CG = None;
            AT = None; AC = None; AG = None;
            GT = None; GC = None; GA = None

        else:
            R = R[0].split("\n")
            first = R[1]
            first = rate_2.findall(first)
            TC = first[1]; TA = first[2]; TG = first[3]
            second = R[2]
            second = rate_2.findall(second)
            CT = second[0]; CA = second[2]; CG = second[3]
            third = R[3]
            third = rate_2.findall(third)
            AT = third[0]; AC = third[1]; AG = third[3]
            fourth = R[4]
            fourth = rate_2.findall(fourth)
            GT = fourth[0]; GC = fourth[1]; GA = fourth[2]
            
         
        df = df.append({"mlb_file_name":mlb_file_name, "basename":basename, "family":family, "protein":protein, "group":filename, "model":model, "lnL":L,
                                 "freq_T":freq_T, "freq_C":freq_C, "freq_A":freq_A,
                                 "freq_G":freq_G, "TC":TC, "TA":TA, "TG":TG, 
                                 "CT":CT, "CA":CA, "CG":CG, "AT":AT,
                                 "AC":AC, "AG":AG, "GT":GT, "GC":GC, "GA":GA},
                                    ignore_index = True)

    df.to_csv(output)
    return(output)


def write_ctl_codeml_file(ctl, seq, tree, out):
    """
    writes a ctl file for codeml PAML
    :param ctl: ctl file output path
    :param seq: sequence file path (phylip)
    :param tree: tree file path
    :param out: output file path
    :param model: model to run
    :param fix_alpha: fix_alpha value
    :param alpha: alpha value
    :return: ctl file path
    """
    ctl = check_filename(ctl, Truefile=False)
    seq = check_filename(seq)
    tree = check_filename(tree)
    out =check_filename(out, Truefile=False)


    ctl_file = open(ctl, "w")

    noisy = "3"
    verbose = "1"
    runmode = "0"
    seqtype = "1"
    CodonFreq = "2"
    ndata = "1"  # was 100
    clock = "0"  # was 1
    aaDist = "0"
    model = "0"
    NSsites = "3"
    icode = "0"
    Mgene = "0"
    fix_kappa = "0"
    kappa = "0.3"
    fix_omega = "0"
    omega = "1.3"
    fix_alpha = "1"
    alpha = "0"
    Malpha = "0"
    ncatG = "3"
    fix_rho = "1"
    rho = "0"
    getSE = "0"
    RateAncestor = "0"
    Small_Diff = ".5e-6"
    cleandata = "1"
    fix_blength = "0"
    method = "0"

    ctl_file.write("\tseqfile = %s\n" % seq)
    ctl_file.write("\ttreefile = %s\n" % tree)
    ctl_file.write("\toutfile = %s\n\n" % out)
    ctl_file.write("\tnoisy = %s\n" % noisy)
    ctl_file.write("\tverbose = %s\n" % verbose)
    ctl_file.write("\trunmode = %s\n" % runmode)
    ctl_file.write("\tseqtype = %s\n" % seqtype)
    ctl_file.write("\tCodonFreq = %s\n" % CodonFreq)
    ctl_file.write("*\tndata = %s\n" % ndata)
    ctl_file.write("\tclock = %s\n" % clock)
    ctl_file.write("\taaDist = %s\n" % aaDist)
    ctl_file.write("\tmodel = %s\n" % model)
    ctl_file.write("\tNSsites = %s\n" % NSsites)
    ctl_file.write("\ticode = %s\n" % icode)
    ctl_file.write("\tMgene = %s\n" % Mgene)
    ctl_file.write("\tfix_kappa = %s\n" % fix_kappa)
    ctl_file.write("\tkappa = %s\n" % kappa)
    ctl_file.write("\tfix_omega = %s\n" % fix_omega)
    ctl_file.write("\tomega = %s\n" % omega)
    #ctl_file.write("\tfix_alpha = %s\n" % fix_alpha)
    #ctl_file.write("\talpha = %s\n" % alpha)
    #ctl_file.write("\tMalpha = %s\n" % Malpha)
    ctl_file.write("\tncatG = %s\n" % ncatG)
    #ctl_file.write("\tfix_rho = %s\n" % fix_rho)
    #ctl_file.write("\trho = %s\n" % rho)
    ctl_file.write("\tgetSE = %s\n" % getSE)
    ctl_file.write("\tRateAncestor = %s\n" % RateAncestor)
    ctl_file.write("\tSmall_Diff = %s\n" % Small_Diff)
    ctl_file.write("\tcleandata = %s\n" % cleandata)
    ctl_file.write("\tfix_blength = %s\n" % fix_blength)
    #ctl_file.write("\tmethod = %s\n" % method)


    ctl_file.close()
    return ctl

def candida_write_ctl_codeml_file(ctl, seq, tree, out, fix_omega):
    """
    writes a ctl file for codeml PAML
    :param ctl: ctl file output path
    :param seq: sequence file path (phylip)
    :param tree: tree file path
    :param out: output file path
    :return: ctl file path
    """
    ctl = check_filename(ctl, Truefile=False)
    seq = check_filename(seq)
    tree = check_filename(tree)
    out =check_filename(out, Truefile=False)


    ctl_file = open(ctl, "w")

    noisy = 9
    verbose = 1
    runmode = 0
    seqtype = 1
    CodonFreq = 2
    clock = 0
    aaDist = 0
    model = 0
    NSsites = 8
    icode = 8
    Mgene = 0
    fix_kappa = 0
    kappa = 2
    fix_omega = fix_omega
    omega = 1
    fix_alpha = 1
    alpha = 0
    Malpha = 0
    ncatG = 3
    fix_rho = 1
    rho = 0
    getSE = 0
    RateAncestor = 0
    Small_Diff = .5e-6
    method = 0

    ctl_file.write("\tseqfile = %s\n" % seq)
    ctl_file.write("\ttreefile = %s\n" % tree)
    ctl_file.write("\toutfile = %s\n\n" % out)
    ctl_file.write("\tnoisy = %s\n" % noisy)
    ctl_file.write("\tverbose = %s\n" % verbose)
    ctl_file.write("\trunmode = %s\n" % runmode)
    ctl_file.write("\tseqtype = %s\n" % seqtype)
    ctl_file.write("\tCodonFreq = %s\n" % CodonFreq)
    #ctl_file.write("*\tndata = %s\n" % ndata)
    ctl_file.write("\tclock = %s\n" % clock)
    ctl_file.write("\taaDist = %s\n" % aaDist)
    ctl_file.write("\tmodel = %s\n" % model)
    ctl_file.write("\tNSsites = %s\n" % NSsites)
    ctl_file.write("\ticode = %s\n" % icode)
    ctl_file.write("\tMgene = %s\n" % Mgene)
    ctl_file.write("\tfix_kappa = %s\n" % fix_kappa)
    ctl_file.write("\tkappa = %s\n" % kappa)
    ctl_file.write("\tfix_omega = %s\n" % fix_omega)
    if fix_omega == 1:
        ctl_file.write("\tomega = %s\n" % omega)
    ctl_file.write("\tfix_alpha = %s\n" % fix_alpha)
    ctl_file.write("\talpha = %s\n" % alpha)
    ctl_file.write("\tMalpha = %s\n" % Malpha)
    ctl_file.write("\tncatG = %s\n" % ncatG)
    ctl_file.write("\tfix_rho = %s\n" % fix_rho)
    ctl_file.write("\trho = %s\n" % rho)
    ctl_file.write("\tgetSE = %s\n" % getSE)
    ctl_file.write("\tRateAncestor = %s\n" % RateAncestor)
    ctl_file.write("\tSmall_Diff = %s\n" % Small_Diff)
    #ctl_file.write("\tcleandata = %s\n" % cleandata)
    #ctl_file.write("\tfix_blength = %s\n" % fix_blength)
    ctl_file.write("\tmethod = %s\n" % method)


    ctl_file.close()
    return ctl


def candida_write_ctl_codeml_branch_site_file(ctl, seq, tree, out, fix_omega):
    """
    writes a ctl file for codeml PAML
    :param ctl: ctl file output path
    :param seq: sequence file path (phylip)
    :param tree: tree file path
    :param out: output file path
    :param fix_omega: 1 for the null model and 0 for the alternative
    :return: ctl file path
    """
    ctl = check_filename(ctl, Truefile=False)
    seq = check_filename(seq)
    tree = check_filename(tree)
    out =check_filename(out, Truefile=False)


    ctl_file = open(ctl, "w")

    noisy = 9
    verbose = 1
    runmode = 0
    seqtype = 1
    CodonFreq = 2
    clock = 0
    aaDist = 0
    model = 2
    NSsites = 2
    icode = 8
    Mgene = 0
    fix_kappa = 0
    kappa = 2
    fix_omega = fix_omega
    omega = 1       # will be written only if fix_omega == 1
    fix_alpha = 1
    alpha = 0
    Malpha = 0
    ncatG = 3
    fix_rho = 1
    rho = 0
    getSE = 0
    RateAncestor = 0
    Small_Diff = .5e-6
    method = 0

    ctl_file.write("\tseqfile = %s\n" % seq)
    ctl_file.write("\ttreefile = %s\n" % tree)
    ctl_file.write("\toutfile = %s\n\n" % out)
    ctl_file.write("\tnoisy = %s\n" % noisy)
    ctl_file.write("\tverbose = %s\n" % verbose)
    ctl_file.write("\trunmode = %s\n" % runmode)
    ctl_file.write("\tseqtype = %s\n" % seqtype)
    ctl_file.write("\tCodonFreq = %s\n" % CodonFreq)
    #ctl_file.write("*\tndata = %s\n" % ndata)
    ctl_file.write("\tclock = %s\n" % clock)
    ctl_file.write("\taaDist = %s\n" % aaDist)
    ctl_file.write("\tmodel = %s\n" % model)
    ctl_file.write("\tNSsites = %s\n" % NSsites)
    ctl_file.write("\ticode = %s\n" % icode)
    ctl_file.write("\tMgene = %s\n" % Mgene)
    ctl_file.write("\tfix_kappa = %s\n" % fix_kappa)
    ctl_file.write("\tkappa = %s\n" % kappa)
    ctl_file.write("\tfix_omega = %s\n" % fix_omega)
    if fix_omega == 1:
        ctl_file.write("\tomega = %s\n" % omega)
    ctl_file.write("\tfix_alpha = %s\n" % fix_alpha)
    ctl_file.write("\talpha = %s\n" % alpha)
    ctl_file.write("\tMalpha = %s\n" % Malpha)
    ctl_file.write("\tncatG = %s\n" % ncatG)
    ctl_file.write("\tfix_rho = %s\n" % fix_rho)
    ctl_file.write("\trho = %s\n" % rho)
    ctl_file.write("\tgetSE = %s\n" % getSE)
    ctl_file.write("\tRateAncestor = %s\n" % RateAncestor)
    ctl_file.write("\tSmall_Diff = %s\n" % Small_Diff)
    #ctl_file.write("\tcleandata = %s\n" % cleandata)
    #ctl_file.write("\tfix_blength = %s\n" % fix_blength)
    ctl_file.write("\tmethod = %s\n" % method)

    ctl_file.close()
    return ctl

def candida_write_ctl_codeml_clade_model(ctl, seq, tree, out, is_null):
    """
    writes a ctl file for codeml PAML
    :param ctl: ctl file output path
    :param seq: sequence file path (phylip)
    :param tree: tree file path
    :param out: output file path
    :param is_null: 1 for the null model and 0 for the alternative
    :return: ctl file path
    """
    ctl = check_filename(ctl, Truefile=False)
    seq = check_filename(seq)
    tree = check_filename(tree)
    out =check_filename(out, Truefile=False)


    ctl_file = open(ctl, "w")

    noisy = 9
    verbose = 1
    runmode = 0
    seqtype = 1
    CodonFreq = 2
    clock = 0
    aaDist = 0
    model = 0 if is_null == 1 else 3
    NSsites = 22 if is_null == 1 else 2
    icode = 8
    Mgene = 0
    fix_kappa = 0
    kappa = 2
    fix_omega = 0
    omega = 0.5
    fix_alpha = 1
    alpha = 0
    Malpha = 0
    ncatG = 3
    fix_rho = 1
    rho = 0
    getSE = 0
    RateAncestor = 0
    Small_Diff = .5e-6
    method = 0

    ctl_file.write("\tseqfile = %s\n" % seq)
    ctl_file.write("\ttreefile = %s\n" % tree)
    ctl_file.write("\toutfile = %s\n\n" % out)
    ctl_file.write("\tnoisy = %s\n" % noisy)
    ctl_file.write("\tverbose = %s\n" % verbose)
    ctl_file.write("\trunmode = %s\n" % runmode)
    ctl_file.write("\tseqtype = %s\n" % seqtype)
    ctl_file.write("\tCodonFreq = %s\n" % CodonFreq)
    #ctl_file.write("*\tndata = %s\n" % ndata)
    ctl_file.write("\tclock = %s\n" % clock)
    ctl_file.write("\taaDist = %s\n" % aaDist)
    ctl_file.write("\tmodel = %s\n" % model)
    ctl_file.write("\tNSsites = %s\n" % NSsites)
    ctl_file.write("\ticode = %s\n" % icode)
    ctl_file.write("\tMgene = %s\n" % Mgene)
    ctl_file.write("\tfix_kappa = %s\n" % fix_kappa)
    ctl_file.write("\tkappa = %s\n" % kappa)
    ctl_file.write("\tfix_omega = %s\n" % fix_omega)
    ctl_file.write("\tomega = %s\n" % omega)
    ctl_file.write("\tfix_alpha = %s\n" % fix_alpha)
    ctl_file.write("\talpha = %s\n" % alpha)
    ctl_file.write("\tMalpha = %s\n" % Malpha)
    ctl_file.write("\tncatG = %s\n" % ncatG)
    ctl_file.write("\tfix_rho = %s\n" % fix_rho)
    ctl_file.write("\trho = %s\n" % rho)
    ctl_file.write("\tgetSE = %s\n" % getSE)
    ctl_file.write("\tRateAncestor = %s\n" % RateAncestor)
    ctl_file.write("\tSmall_Diff = %s\n" % Small_Diff)
    #ctl_file.write("\tcleandata = %s\n" % cleandata)
    #ctl_file.write("\tfix_blength = %s\n" % fix_blength)
    ctl_file.write("\tmethod = %s\n" % method)

    ctl_file.close()
    return ctl


def write_ctl_codeml_deleterious_load_file(ctl, seq, tree, out, model):
    """
    writes a ctl file for codeml PAML
    :param ctl: ctl file output path
    :param seq: sequence file path (phylip)
    :param tree: tree file path
    :param out: output file path
    :param model: model to run - 0 for 1 rate, 2 for branch specific rate
    :param fix_alpha: fix_alpha value
    :param alpha: alpha value
    :return: ctl file path
    """
    ctl = check_filename(ctl, Truefile=False)
    seq = check_filename(seq)
    tree = check_filename(tree)
    out =check_filename(out, Truefile=False)


    ctl_file = open(ctl, "w")

    noisy = "3"
    verbose = "1"
    runmode = "0"
    seqtype = "1"
    CodonFreq = "2"
    clock = "0"
    aaDist = "0"
    NSsites = "0"
    icode = "0"
    Mgene = "0"
    fix_kappa = "0"
    kappa = "0.3"
    fix_omega = "0"
    omega = "1.3"
    fix_alpha = "1"
    alpha = "0"
    ncatG = "3"
    getSE = "0"
    Small_Diff = ".5e-6"
    method = "0"

    ctl_file.write("\tseqfile = %s\n" % seq)
    ctl_file.write("\ttreefile = %s\n" % tree)
    ctl_file.write("\toutfile = %s\n\n" % out)
    ctl_file.write("\tnoisy = %s\n" % noisy)
    ctl_file.write("\tverbose = %s\n" % verbose)
    ctl_file.write("\trunmode = %s\n" % runmode)
    ctl_file.write("\tseqtype = %s\n" % seqtype)
    ctl_file.write("\tCodonFreq = %s\n" % CodonFreq)
    ctl_file.write("\tclock = %s\n" % clock)
    ctl_file.write("\taaDist = %s\n" % aaDist)
    ctl_file.write("\tmodel = %s\n" % model)
    ctl_file.write("\tNSsites = %s\n" % NSsites)
    ctl_file.write("\ticode = %s\n" % icode)
    ctl_file.write("\tMgene = %s\n" % Mgene)
    ctl_file.write("\tfix_kappa = %s\n" % fix_kappa)
    ctl_file.write("\tkappa = %s\n" % kappa)
    ctl_file.write("\tfix_omega = %s\n" % fix_omega)
    ctl_file.write("\tomega = %s\n" % omega)
    ctl_file.write("\tfix_alpha = %s\n" % fix_alpha)
    ctl_file.write("\talpha = %s\n" % alpha)
    ctl_file.write("\tncatG = %s\n" % ncatG)
    ctl_file.write("\tgetSE = %s\n" % getSE)
    ctl_file.write("\tSmall_Diff = %s\n" % Small_Diff)
    ctl_file.write("\tmethod = %s\n" % method)


    ctl_file.close()
    return ctl


def retrive_codeml_delitirious_restults(files, output="./codeml_1-2ratio_results.csv"):
    """

    :param files: list of fasta / phy files.
        the mlc files are in 1 and 2 ratio directories in directory as the fasta/phy files
    :param output: output file. default ./codeml_1-2ratio_results.csv
    :return: dataframe with results
    """
    files = [check_filename(f) for f in files]
    output = check_filename(output, Truefile=False)
    dir = os.path.dirname(files[0]).split("codon_aln")[0]
    # regex patterns for retrieving results
    omega_1ratio_pattern = re.compile("omega.*")
    dN_length_pattern = re.compile("tree length for dN:.*")
    dS_length_pattern = re.compile("tree length for dS:.*")
    omega_branches_2ratio_pattern = re.compile(".*for branches:.*")
    lnL_pattern = re.compile("lnL.*")

    codeml_results = pd.DataFrame()
    # go over files
    for f in files:
        if "unaligned" in f:
            continue
        basename = f.split("/")[-1].split(".fasta")[0].split(".phy")[0].split(".aln")[0]
        print("%s/1ratio/%s*.cml" % (dir, basename))
        print("%s/2ratio/%s*.cml" % (dir, basename))
        mlc_1ratio = glob.glob("%s/1ratio/%s*.cml" % (dir, basename))[0]
        mlc_2ratio = glob.glob("%s/2ratio/%s*.cml" % (dir, basename))[0]
        omega = None
        dN_length_1ratio = None
        dS_length_1ratio = None
        number_of_seqs = None
        length_of_seqs_codons = None
        dN_length_2ratio = None
        dS_length_2ratio = None
        external_branchs = None
        internal_branchs = None
        lnL_1ratio = None
        lnL_2ratio = None

        # get information from 1 ratio result file
        res_1ratio = open(mlc_1ratio, "r").read()
        if "Time" in res_1ratio:
            pattern = re.compile("omega.*")
            omega = float(omega_1ratio_pattern.findall(res_1ratio)[0].split("=")[-1].strip())
            dN_length_1ratio = float(dN_length_pattern.findall(res_1ratio)[0].split(":")[-1].strip())
            dS_length_1ratio = float(dS_length_pattern.findall(res_1ratio)[0].split(":")[-1].strip())
            sequence_info = [i for i in res_1ratio.split("\n")[0].split(" ") if i != ""]
            number_of_seqs = float(sequence_info[0])
            length_of_seqs_codons = float(sequence_info[1]) / 3
            lnL_1ratio = float(lnL_pattern.findall(res_1ratio)[0].split("):")[-1].strip().split(" ")[0])

        res_2ratio = open(mlc_2ratio, "r").read()
        if "Time" in res_2ratio:
            dN_length_2ratio = float(dN_length_pattern.findall(res_2ratio)[0].split(":")[-1].strip())
            dS_length_2ratio = float(dS_length_pattern.findall(res_2ratio)[0].split(":")[-1].strip())
            omega_branches_2ratio = omega_branches_2ratio_pattern.findall(res_2ratio)[0].split(" ")[-2:]
            external_branchs = float(omega_branches_2ratio[0])
            internal_branchs = float(omega_branches_2ratio[1])
            lnL_2ratio = float(lnL_pattern.findall(res_2ratio)[0].split("):")[-1].strip().split(" ")[0])
        #write results to dataframe
        codeml_results = codeml_results.append(
            {"filename": basename, "num_of_seqs": number_of_seqs, "seg_len_codon": length_of_seqs_codons,
             "lnL_1ratio": lnL_1ratio, "lnL_2ratio": lnL_2ratio,
             "w": omega, "wi": internal_branchs, "we": external_branchs,
             "dN_length_1ratio": dN_length_1ratio, "dS_length_1ratio": dS_length_1ratio,
             "dN_length_2ratio": dN_length_2ratio, "dS_length_2ratio": dS_length_2ratio}, ignore_index=True)
    print(output)
    codeml_results.to_csv(output)
    return codeml_results


def retrive_kappa_from_paml_output(mlb):
    """
    return kappa
    :param mlb: mlb file
    :return: kappa
    """
    mlb = check_filename(mlb)
    with open(mlb, "r") as mlb_handle:
        mlb_data = mlb_handle.read()
        loc = mlb_data.find("Average Ts/Tv")
        try:
            kappa = float(mlb_data[loc:loc+50].split("=")[1].split("\n")[0].strip())
            return kappa
        except:
            return None