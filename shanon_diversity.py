#! /usr/local/python_anaconda/bin/python3.4

import os
import glob
from file_utilities import check_dirname, check_filename
import pandas as pd
import re
import math


def creats_shannon_index_lst_nuc(df):
    """creates a list with shannon indexs for a given freqs.

    :param input_file: freq file that was genrated from the pipline in csv format only!
    :return: list of shannon indexs

    """
    sample = list(set(list(df["sample"])))[0]
    segment = list(set(list(df["segment"])))[0]
    experiment =  list(set(list(df["experiment"])))[0]

    lst_of_frqs = []
    lst_by_4 = []
    pos = list(df["Pos"])[0]
    shannon_df = pd.DataFrame(columns = ["Pos", "segment", "sample", "experiment", "shanon"])
    for index, row in df.iterrows():
        pos = row["Pos"]
        lst_by_4.append(row["Freq"])

        if len(lst_by_4) == 4:
            if lst_by_4 == [0.0, 0.0, 0.0, 0.0]:
                lst_by_4 = []
                continue
            shannon = calc_shannon_entropy(lst_by_4)
            shannon_df = shannon_df.append({"Pos": pos, "segment": segment,
                                            "sample": sample,
                                            "experiment": experiment,
                                            "shanon": shannon}, ignore_index=True)
            lst_by_4 = []


    return shannon_df



def creats_shannon_index_lst_aa(df):
    """creates a list with shannon indexs for a given freqs.

    :param input_file: freq file that was genrated from the pipline in csv format only!
    :return: list of shannon indexs

    """
    sample = list(set(list(df["sample"])))[0]
    segment = list(set(list(df["segment"])))[0]
    experiment =  list(set(list(df["experiment"])))[0]

    lst_of_frqs = []
    lst_by_12 = {}
    count = 0
    pos = list(df["Pos"])[0]
    shannon_df = pd.DataFrame(columns = ["Pos", "segment", "sample", "experiment", "shanon"])
    for index, row in df.iterrows():
        count +=1
        pos = row["Pos"] / 3
        aa = row["mut_aa"]
        if aa != "*":
            if aa not in lst_by_12.keys():
                lst_by_12[aa] = 0
            lst_by_12[aa] += row["Freq"]

        if count == 12:
            #if lst_by_12 == [0.0, 0.0, 0.0, 0.0]:
            #    lst_by_12 = {}
            #    continue
            sum_freqs = round(sum(lst_by_12.values()))
            if sum_freqs != 3:
                lst_by_12 = {}
                count = 0

            shannon = calc_shannon_entropy(lst_by_12.values(), bases_space_size=23)
            shannon_df = shannon_df.append({"Pos": pos, "segment": segment,
                                            "sample": sample,
                                            "experiment": experiment,
                                            "shanon": shannon}, ignore_index=True)
            print(lst_by_12)
            print(shannon)
            lst_by_12 = {}
            count = 0

    return shannon_df


def calc_shannon_entropy(counts, bases_space_size=4):
    """Return the Shannon Entropy index for the given counts.

    :param counts: a list of size bases_space_size
    :param bases_space_size: represents the possible size of the counts list.
                             in our code defaults to 4 as the bases of the DNA,
                             but this can be extended to proteins, etc. (default: 4)
    :return: the shannon entropy of the provided count

    """
    entropy = 0
    count_sum = sum(counts)
    freqs = [x / count_sum for x in counts]
    for p_x in freqs:
        if p_x > 0:
            entropy += - p_x * math.log(p_x, bases_space_size)
    return entropy
