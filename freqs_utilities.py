#! /usr/local/python_anaconda/bin/python3.4

import pandas as pd
from file_utilities import check_filename
from Bio.Seq import Seq
from scipy.stats import ttest_ind
import numpy as np
import os
from functools import reduce
from optparse import OptionParser



def remove_problematic_lines(freqs_file, text_to_remove="#Created with join V5.1 May 15 2015\n#"):
    """
    removes problematic lines from freqs files and saves the files in the same location
    :param freqs_file: input freqs file path
    :param text_to_remove: text to remove from file (default: "#Created with join V5.1 May 15 2015\n#")
    :return:
    """
    freqs_file = check_filename(freqs_file)
    freqs = open(freqs_file, "r").read().split("#Created with join V5.1 May 15 2015\n#")[-1]
    freqs_out = open(freqs_file, "w")
    freqs_out.write(freqs)
    freqs_out.close()
    save_freqs_with_different_delimiter(freqs_file)
    return freqs_file

# 
def save_freqs_with_different_delimiter(freqs_file):
    freqs_file = check_filename(freqs_file)
    f = pd.read_csv(freqs_file, sep="\t")
    f.to_csv(freqs_file, index=False)


def merge_freqs_files(freqs_files, output):
    """
    merges freqs files of TILV - the file names are in the following format: P%i-S%i % (passage, segment)
    :param freqs_files: list of input freqs file paths
    :param output: output merged freqs file path
    :return:
    """
    freqs_files = [check_filename(f) for f in freqs_files]
    output = check_filename(output, Truefile=False)
    all = pd.DataFrame()
    # 
    for freqs_file in freqs_files:
        time = freqs_file.split("/")[-1].split(".")[0].split("-")[0]
        segment = freqs_file.split("/")[-1].split(".")[0].split("-")[1].split("S")[1].split("_")[0]
        f = pd.read_csv(freqs_file)
        f.loc[:, "Time"] = time
        f.loc[:, "Segment"] = segment
        f.loc[:, "Real_Pos"] = f.Pos
        f.loc[:, "Pos"] = f.loc[:, "Pos"].astype("int")
        f.loc[:, "Pos"] = f.loc[:, "Pos"].astype("str")
        f.loc[:, "Pos"] = f[["Pos", "Segment"]].apply(lambda x: ''.join(x), axis=1)
        all = pd.concat([all, f])
    all.to_csv(output, index=False)
    return output, all

def change_ref_to_consensus(freqs_file):
    freqs = pd.read_csv(freqs_file, sep="\t")
    freqs_tmp = pd.read_csv(freqs_file, sep="\t")
    freqs_tmp = freqs_tmp.drop_duplicates("Pos")
    freqs_tmp = freqs_tmp[['Pos', 'Base']]
    transformed_freq = freqs.set_index(freqs.Pos).join(
        freqs_tmp.set_index(freqs_tmp.Pos), rsuffix='_r')
    transformed_freq.Ref = transformed_freq.Base_r
    transformed_freq = transformed_freq.loc[(transformed_freq.Pos >= 1456) & (transformed_freq.Pos <= 4488)][
        ['Pos', 'Base', 'Freq', 'Ref', 'Read_count', 'Rank', 'Prob']]

    return transformed_freq



def add_mutation_to_freq_file(output, freqs_file = None, freqs = None, forced_rf_shift = 0):
    # assumes that position 1 is the beginning of the CDS
    # removes positions that at the beginning or at the end that are not part of a full codon
    if freqs_file == None and type(freqs) == "NoneType":
        raise Exception("Need to specify or freqs file path or a freqs pandas object")
    elif freqs_file != None and freqs != None:
        print(freqs_file, freqs)
        print(type(freqs))
        raise Exception("Need to specify EITHER freqs file path OR a freqs pandas object - only one!")
    elif freqs_file != None:
        freqs = pd.read_csv(freqs_file, sep="\t")
    freqs = freqs[freqs.Pos % 1 == 0] #removes insertions
    freqs = freqs[freqs.Base != "-"] #removes deletions
    freqs.reset_index(drop=True, inplace=True)

    first_pos = int(freqs.loc[1].Pos) #gets the first position in the right frameshift
    first_pos += forced_rf_shift
    if first_pos == 1:
        start_from = first_pos
    elif first_pos % 3 == 1:
        start_from = first_pos
    elif first_pos % 3 == 2:
        start_from = first_pos + 2
    elif first_pos % 3 == 0:
        start_from = first_pos + 1

    freqs["Mutation_type"] = None
    freqs["wt_aa"] = None
    freqs["mut_aa"] = None
    freqs["wt_codon"] = None
    freqs["mut_codon"] = None
    freqs["Mutation"] = None

    for pos in range(start_from, int(max(freqs.Pos)), 3): # add mutation information
        temp = freqs.loc[freqs['Pos'].isin([pos, pos+1, pos+2])]
        if len(temp) != 12: # 12 - is 3 positions * 4 base mutations [A, C, G, T] (ordered by Rank)
            continue

        first = temp.iloc[0].Ref
        second = temp.iloc[4].Ref
        third = temp.iloc[8].Ref
        wt_codon = "".join([first, second, third])
        wt_aa = str(Seq(wt_codon).translate())

        pos = temp.iloc[0].Pos
        for n in range(0, 12):
            ref_base = temp.iloc[n].Ref
            mut_base = temp.iloc[n].Base
            if n <= 3:
                mut_codon = "".join([mut_base, second, third])
                current_pos = pos
            elif n > 3 and n <= 7:
                mut_codon = "".join([first, mut_base, third])
                current_pos = pos + 1
            elif n > 7  :
                mut_codon = "".join([first, second, mut_base])
                current_pos = pos + 2

            mut_aa = str(Seq(mut_codon).translate())

            if wt_codon == mut_codon:
                mutation_type = "consensus"
            elif wt_aa == mut_aa:
                mutation_type = "synonymous"
            elif wt_aa != "*" and mut_aa == "*":
                mutation_type = "stop"
            else:

                mutation_type = "missense"

            freqs.loc[(freqs["Pos"] == current_pos) & (freqs["Base"] == mut_base), "Mutation_type"] = mutation_type
            freqs.loc[(freqs["Pos"] == current_pos) & (freqs["Base"] == mut_base), "wt_aa"] = wt_aa
            freqs.loc[(freqs["Pos"] == current_pos) & (freqs["Base"] == mut_base), "mut_aa"] = mut_aa
            freqs.loc[(freqs["Pos"] == current_pos) & (freqs["Base"] == mut_base), "wt_codon"] = wt_codon
            freqs.loc[(freqs["Pos"] == current_pos) & (freqs["Base"] == mut_base), "mut_codon"] = mut_codon
            freqs.loc[(freqs["Pos"] == current_pos) & (freqs["Base"] == mut_base), "Mutation"] = ref_base + mut_base

    freqs = freqs[freqs.Mutation_type.notnull()] #removes Nones - rows at the beginning and the end
    freqs.to_csv(output, index=False, sep='\t')
    return freqs


def filter_freqs_for_regression_analysis(output, freqs_file = None, freqs = None, filter = 1, control_cutoff=0.000185):
    if freqs_file == None and type(freqs) == "NoneType":
        raise Exception("Need to specify or freqs file path or a freqs pandas object")
    elif freqs_file != None and freqs != None:
        raise Exception("Need to specify or freqs file path OR a freqs pandas object - only one!")
    elif freqs_file != None:
        freqs = pd.read_csv(freqs_file)
    freqs = freqs[freqs.Mutation_type=="stop"]

    freqs = freqs[freqs.Read_count >= 10000] #read count over 10,000
    freqs = freqs[freqs.Real_Pos % 1 == 0] #removes insertions
    freqs = freqs[freqs.Base != "-"] #removes deletions
    freqs = freqs[freqs.Base != freqs.Ref] #removes consensus positions
    freqs = freqs[(freqs.Prob >= 0.95)]  # reads with low probabilities
    """

    control_cutoff = 0.000185
    transition_cutoff = 0.00020485
    transversion_cutoff = 0.0001355
    transitions = ["AG", "GA", "CT", "TC"]
    transversion = ["AC", "CA", "AT", "TA", "GC", "CG", "GT", "TG"]
    transition_transversion = False
    if filter == 1:
        P2_lower_boundry = freqs["Freq"][freqs.Time == "P2"].quantile(0.25)
        P2_upper_boundry = freqs["Freq"][freqs.Time == "P2"].quantile(0.75)
        P30_lower_boundry = freqs["Freq"][freqs.Time == "P30"].quantile(0.25)
        P30_upper_boundry = freqs["Freq"][freqs.Time == "P30"].quantile(0.75)
    elif filter == 2:
        P2_lower_boundry = control_cutoff
        P2_upper_boundry = freqs["Freq"][freqs.Time == "P2"].quantile(0.75)
        P30_lower_boundry = control_cutoff
        P30_upper_boundry = freqs["Freq"][freqs.Time == "P30"].quantile(0.75)
    elif filter == 3:
        P2_lower_boundry = control_cutoff
        P2_upper_boundry = freqs["Freq"][freqs.Time == "P2"].quantile(0.90)
        P30_lower_boundry = control_cutoff
        P30_upper_boundry = freqs["Freq"][freqs.Time == "P30"].quantile(0.90)
    elif filter == 4:
        P2_lower_boundry = control_cutoff
        P2_upper_boundry = 1
        P30_lower_boundry = control_cutoff
        P30_upper_boundry = 1
    elif filter == 5:
        transition_transversion = True
        P2_upper_boundry = freqs["Freq"][freqs.Time == "P2"].quantile(0.75)
        P30_upper_boundry = freqs["Freq"][freqs.Time == "P30"].quantile(0.75)
    elif filter == 6:
        transition_transversion = True
        P2_upper_boundry = freqs["Freq"][freqs.Time == "P2"].quantile(0.90)
        P30_upper_boundry = freqs["Freq"][freqs.Time == "P30"].quantile(0.90)
    elif filter == 7:
        transition_transversion = True
        P2_upper_boundry = 1
        P30_upper_boundry = 1
    else:
        raise Exception ("unknown filter")
    
    if not transition_transversion:
        filtered = freqs[(freqs.Time == "P2") & (freqs.Freq >= P2_lower_boundry) & (freqs.Freq <= P2_upper_boundry) |
                     (freqs.Time == "P30") & (freqs.Freq >= P30_lower_boundry) & (freqs.Freq <= P30_upper_boundry)]
    else:
        filtered = freqs[(freqs.Time == "P2")  & (freqs.Freq <= P2_upper_boundry) |
                         (freqs.Time == "P30") & (freqs.Freq <= P30_upper_boundry)]
        filtered = filtered[((filtered.Mutation.isin(transitions)) & (filtered.Freq >= transition_cutoff)) |
                            ((filtered.Mutation.isin(transversion)) & (filtered.Freq >= transversion_cutoff))]
    filtered.is_copy = False
    """
    filtered = freqs
    position_counts = filtered.groupby(["Base", "Pos"]).size().reset_index(name="Count")


    bases = list(position_counts.Base[position_counts.Count ==2])
    positions = list(position_counts.Pos[position_counts.Count ==2])

    positions_in_both_times = ["".join([bases[i], str(positions[i])]) for i in range(len(positions))]


    filtered.loc[:,"NEW"] = filtered.loc[:, "Base"] +  filtered.loc[:,"Pos"].astype("str")

    filtered = filtered[filtered.NEW.isin(positions_in_both_times)]
    filtered.to_csv(output, index=False)

    return output, filtered



def calculate_control_cutoff(freq_control_file, output):
    freq_control_file = check_filename(freq_control_file)
    output = check_filename(output, Truefile=False)
    control = pd.read_csv(freq_control_file)
    deletions = control[control.Base == "-"]
    deletions = deletions[deletions.Read_count >= 10000]
    deletions = deletions[(deletions.Prob >= 0.95) | (deletions.Prob <= 0)]  # reads with low probabilities
    control = control[control.Pos % 1 == 0] #removes insertions
    control = control[control.Base != "-"] #removes deletions
    control["Mutation"] = control.Ref + control.Base
    control = control[control.Read_count >= 10000]
    control = control[(control.Prob >= 0.95) | (control.Prob <= 0)]  # reads with low probabilities
    control = control[control.Ref != control.Base]



    all_control = control.Freq.quantile(0.95)
    transition_control = control.Freq[control.Mutation.isin(["AG", "GA", "CT", "TC"])].quantile(0.95)
    transversion_control = control.Freq[-control.Mutation.isin(["AG", "GA", "CT", "TC"])].quantile(0.95)
    deletions_cutoff = deletions.Freq.quantile(0.95)
    print(deletions_cutoff)
    text = "all control: %f\ntransition control: %f\ntransversion control:%f" % (all_control, transition_control,
                                                                                 transversion_control)
    out = open(output, "w")
    out.write(text)
    out.close()
    return all_control, transition_control, transversion_control



def filter_freqs_for_stop_analysis(output, freqs_file = None, freqs = None):
    if freqs_file == None and type(freqs) == "NoneType":
        raise Exception("Need to specify or freqs file path or a freqs pandas object")
    elif freqs_file != None and freqs != None:
        raise Exception("Need to specify or freqs file path OR a freqs pandas object - only one!")
    elif freqs_file != None:
        freqs = pd.read_csv(freqs_file)

    print(freqs_file)
    freqs = freqs[freqs.Real_Pos % 1 == 0] #removes insertions
    freqs = freqs[freqs.Base != "-"] #removes deletions

    quantile_05 = freqs.Freq[(freqs.Time == "P2") & (freqs.Ref == freqs.Base)].quantile(0.05)
    positions_to_remove = list(freqs.Real_Pos[(freqs.Time == "P2") &(freqs.Ref == freqs.Base) & (freqs.Freq < quantile_05)])
    freqs = freqs[-freqs.Real_Pos.isin(positions_to_remove)]
    freqs = freqs[freqs.Read_count >= 10000]  # read count over 10,000
    freqs = freqs[(freqs.Prob >= 0.95) | (freqs.Prob <= 0)]  # reads with low probabilities
    freqs = freqs[freqs.Base != freqs.Ref]  # removes consensus positions

    stop = freqs[freqs.Mutation_type == "stop"]
    stop.to_csv(output, index=False)

    return output, stop


def t_test_frequencies(freqs_1, freqs_2):
    #one sided - freqs_1 is less than freqs_2
    t, p = ttest_ind(freqs_1.apply(np.log10), freqs_2.apply(np.log10), equal_var=False)
    p = p/2
    print(t, p)
    significant = p<0.05 and t < 0
    if significant:
        print("Significant: %f" % p)
    else:
        print("Non Significant: %f" % p)


def mutation_accomulation_analysis(freqs_file = None, freqs = None, indel=False):
    if freqs_file == None and type(freqs) == "NoneType":
        raise Exception("Need to specify or freqs file path or a freqs pandas object")
    elif freqs_file != None and freqs != None:
        raise Exception("Need to specify or freqs file path OR a freqs pandas object - only one!")
    elif freqs_file != None:
        freqs = pd.read_csv(freqs_file)

    maximal_position = freqs.Real_Pos.max()

    control_cutoff = 0.000185
    transition_cutoff = 0.00020485
    transversion_cutoff = 0.0001355
    transitions = ["AG", "GA", "CT", "TC"]
    transversion = ["AC", "CA", "AT", "TA", "GC", "CG", "GT", "TG"]

    if indel:
        control_cutoff = 0

    print("Stop mutation test - whole segment")
    print("all cutoff")
    stop_P2 = freqs.Freq[(freqs.Time == "P2") & (freqs.Freq >= control_cutoff)]
    stop_P30 = freqs.Freq[(freqs.Time == "P30") & (freqs.Freq >= control_cutoff)]
    t_test_frequencies(stop_P2, stop_P30)

    if not indel:
        print("transition/transversion cutoff")
        stop_P2 = freqs.Freq[((freqs.Time == "P2") & (freqs.Freq >= transition_cutoff) & (freqs.Mutation.isin(transitions))) |
                            ((freqs.Time == "P2") & (freqs.Freq >= transversion_cutoff) & (freqs.Mutation.isin(transversion)))]
        stop_P30 = freqs.Freq[((freqs.Time == "P30") & (freqs.Freq >= transition_cutoff) & (freqs.Mutation.isin(transitions))) |
                            ((freqs.Time == "P30") & (freqs.Freq >= transversion_cutoff) & (freqs.Mutation.isin(transversion)))]
        t_test_frequencies(stop_P2, stop_P30)

    freqs = freqs[freqs.Real_Pos <= maximal_position*0.75]
    print("Stop mutation test - the first 75% of the sequence")
    print("all cutoff")
    stop_P2 = freqs.Freq[(freqs.Time == "P2") & (freqs.Freq >= control_cutoff)]
    stop_P30 = freqs.Freq[(freqs.Time == "P30") & (freqs.Freq >= control_cutoff)]

    t_test_frequencies(stop_P2, stop_P30)
    if not indel:
        print("transition/transversion cutoff")
        stop_P2 = freqs.Freq[((freqs.Time == "P2") & (freqs.Freq >= transition_cutoff) & (freqs.Mutation.isin(transitions))) |
                            ((freqs.Time == "P2") & (freqs.Freq >= transversion_cutoff) & (freqs.Mutation.isin(transversion)))]

        stop_P30 = freqs.Freq[((freqs.Time == "P30") & (freqs.Freq >= transition_cutoff) & (freqs.Mutation.isin(transitions))) |
                            ((freqs.Time == "P30") & (freqs.Freq >= transversion_cutoff) & (freqs.Mutation.isin(transversion)))]
        t_test_frequencies(stop_P2, stop_P30)



    print()

def filter_freqs_for_indel_analysis(output, freqs_file = None, freqs = None):
    if freqs_file == None and type(freqs) == "NoneType":
        raise Exception("Need to specify or freqs file path or a freqs pandas object")
    elif freqs_file != None and freqs != None:
        raise Exception("Need to specify or freqs file path OR a freqs pandas object - only one!")
    elif freqs_file != None:
        freqs = pd.read_csv(freqs_file)


    quantile_05 = freqs.Freq[(freqs.Time == "P2") & (freqs.Ref == freqs.Base)].quantile(0.05)
    positions_to_remove = list(freqs.Real_Pos[(freqs.Time == "P2") &(freqs.Ref == freqs.Base) & (freqs.Freq < quantile_05)])
    freqs = freqs[-freqs.Real_Pos.isin(positions_to_remove)]


    freqs = freqs[(freqs.Base == "-") | (freqs.Ref == "-")]


    positions_to_remove = []
    positions_P2 = np.array(freqs.Real_Pos[freqs.Time=="P2"])

    positions_P30 = np.array(freqs.Real_Pos[freqs.Time == "P30"])


    print(len(freqs))

    to_ignore = []

    for p in positions_P2:
        if np.any(to_ignore == p):
            continue
        if p % 3 == 1 and np.any(positions_P2 == p+1) and np.any(positions_P2 == p+2):
            positions_to_remove.append(p)
            to_ignore.append([p+1, p+2])
    to_ignore = []

    for p in positions_P30:
        if np.any(to_ignore == p):
            continue
        if p % 3 == 1 and np.any(positions_P30 == p + 1) and np.any(positions_P30 == p + 2):
            positions_to_remove.append(p)
            to_ignore.append([p + 1, p + 2])

    freqs = freqs[-freqs.Real_Pos.isin(positions_to_remove)]

    indel_cutoff = 9.269999999999992e-05

    freqs = freqs[freqs.Read_count >= 10000]  # read count over 10,000
    freqs = freqs[(freqs.Prob >= 0.95) | (freqs.Prob <= 0)]  # reads with low probabilities
    freqs = freqs[freqs.Freq >= indel_cutoff]


    freqs.to_csv(output, index=False)
    print(len(freqs))
    print()
    return output, freqs


def get_area_of_highest_coverage_and_stop_mutation_potential(freqs_file, size_of_output=1000, output=None):
    freqs_file = check_filename(freqs_file)
    freqs = pd.read_csv(freqs_file)
    if output == None:
        output = os.path.splitext(freqs_file)[0] + "_coverage_and_stop_mutation_potential.csv"


    df = pd.DataFrame(columns = ["start_pos", "coverage_in_pos", "stop_mutation_count_pos",
                                 "GA_stop_mutation_pos", "non_synonymous_transition_mutation_pos",
                                 "median_coverage_1000", "stop_mutation_1000", "GA_stop_mutation_1000",
                                 "non_synonymous_transition_mutation_1000",
                                 "median_coverage_100", "stop_mutation_100", "GA_stop_mutation_100",
                                 "non_synonymous_transition_mutation_100",])

    for i in range(1, freqs["Pos"].max() - 99):
        position = freqs[(freqs["Pos"]==i)]
        median_coverage_pos, stop_mutation_count_pos, GA_stop_mutation_count_pos, ns_transition_mutation_count_pos = \
            get_median_coverage_and_stop_mutation_for_df(position)


        stop_mutation_in_position = (position["Mutation.Type"] == "Premature Stop Codon").sum()
        coverage_in_position = (position["Read_count"]).mean()

        temp_1000 = freqs[(freqs["Pos"] >= i) & (freqs["Pos"] <= i +1000)]
        median_coverage_1000, stop_mutation_count_1000, GA_stop_mutation_count_1000, ns_transition_mutation_count_1000 = \
            get_median_coverage_and_stop_mutation_for_df(temp_1000)

        temp_100 = freqs[(freqs["Pos"] >= i) & (freqs["Pos"] <= i +100)]
        median_coverage_100, stop_mutation_count_100, GA_stop_mutation_count_100, ns_transition_mutation_count_100 = \
            get_median_coverage_and_stop_mutation_for_df(temp_100)


        df = df.append({"start_pos":i, "coverage_in_pos":median_coverage_pos,
                        "stop_mutation_count_pos":stop_mutation_count_pos,
                        "GA_stop_mutation_pos":GA_stop_mutation_count_pos,
                        "non_synonymous_transition_mutation_pos":ns_transition_mutation_count_pos,
                        "median_coverage_1000":median_coverage_1000,
                        "stop_mutation_1000":stop_mutation_count_1000,
                        "GA_stop_mutation_1000":GA_stop_mutation_count_1000,
                        "non_synonymous_transition_mutation_1000":ns_transition_mutation_count_1000,
                        "median_coverage_100":median_coverage_100,
                        "stop_mutation_100":stop_mutation_count_100,
                        "GA_stop_mutation_100":GA_stop_mutation_count_100,
                        "non_synonymous_transition_mutation_100": ns_transition_mutation_count_100}, ignore_index=True)


    print(output)
    df.to_csv(output, index="False")

    return (freqs, df)


def get_median_coverage_and_stop_mutation_for_df(df):
    median_coverage = df["Read_count"].median()
    stop_mutation_count = (df["Mutation.Type"] == "Premature Stop Codon").sum()
    GA_stop_mutation_count =  len(df[(df["Mutation.Type"] =="Premature Stop Codon") & (df["Ref"] =="G") & (df["Base"] =="A")])
    ns_transition_mutation_count = len(df[(df["mutation_in_conseved_area"] == True)])
    return (median_coverage, stop_mutation_count, GA_stop_mutation_count, ns_transition_mutation_count)


def compare_positions_between_freqs(dict_of_freqs, out_path=False, positions_to_compare=False):
    '''
    Compare and display frequencies and read counts per position between different freqs files.
    :param dict_of_freqs: a dictionary of freqs files to compare in the format of {name:freqs_path}
    :param out_path: path to save output csv, optional
    :param positions_to_compare: a list of positions to compare, optional
    :return: merged dataframe
    '''
    dfs = []
    for i in dict_of_freqs:
        df = pd.read_csv(dict_of_freqs[i], sep='\t')
        df = df[['Ref', 'Pos', 'Base', 'Freq', 'Read_count']]
        df.rename(columns={'Read_count':'Read_count_' + i, 'Freq':'Freq_' + i}, inplace=True)
        dfs.append(df)
    df_final = reduce(lambda left,right: pd.merge(left,right,on=['Ref', 'Pos', 'Base']), dfs)
    if positions_to_compare:
        df_final = df_final[(df_final.Pos.isin(positions_to_compare))]
    if out_path:
        df_final.to_csv(out_path, index=False)
    return df_final

    
def estimate_insertion_freq(df, extra_columns=[]):
    '''
    This function gets a freqs file(s) dataframe, calculates the frequency of insertions by using the read count of the
    previous base and returns a dataframe including this.
    :param df: a dataframe of freqs file(s).
    :param extra_columns: if df contains more than the basic freqs columns, for example a column of Sample_id etc., 
    provide a list of the extra columns to be included.
    :return: df with extra columns describing insertion frequency.
    '''
    read_counts = df[(df.Ref != '-')][ extra_columns + ['Pos', 'Read_count']].drop_duplicates()
    read_counts.rename(columns={'Read_count':'estimated_read_count', 'Pos':'rounded_pos'}, inplace=True)
    insertions = df[(df.Ref == '-')]
    not_insertions = df[(df.Ref != '-')]
    insertions['rounded_pos'] = insertions.Pos.astype(int).astype(float)
    insertions = pd.merge(insertions, read_counts, how='left', on= extra_columns + ['rounded_pos'])
    insertions['estimated_freq'] = insertions.Freq * insertions.Read_count / insertions.estimated_read_count
    df = pd.concat([insertions, not_insertions])
    return df.sort_values(extra_columns + ['Pos'])

def main():
    parser = OptionParser("usage: %prog [options]\nTry running %prog --help for more information")
    parser.add_option("-f", "--freqs", dest="freqs", help="frequency file")
    parser.add_option("-o", "--output_freqs", dest="output_file", help="output freqs file with mutations")
    (options, args) = parser.parse_args()
    freq_file = options.freqs
    output_freq_file = options.output_file

    print('Handling file:' + freq_file)  # TODO- remove
    add_mutation_to_freq_file_with_cons_as_ref(output_freq_file, freq_file)


def add_mutation_to_freq_file_with_cons_as_ref(output_freq_file, freq_file):
    transformed_freq = change_ref_to_consensus(freq_file)
    add_mutation_to_freq_file(output_freq_file, freqs= transformed_freq)


if __name__ == "__main__":
    main()
