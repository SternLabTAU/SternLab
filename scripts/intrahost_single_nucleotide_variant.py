#! /usr/local/python_anaconda/bin/python3.4

import glob
from file_utilities import check_dirname, check_filename
import pandas as pd
import re
import math
pd.options.mode.chained_assignment = None

#positions to remove in TiLV - primer positions
TILV_FILTER = {1:[0, 1, 2, 3, 4, 5, 6, 7, 8, 9, 10, 11, 12, 13, 14, 15, 16, 17, 18, 19, 20, 21, 22, 23, 24, 25, 26, 27,
                  28, 29, 30, 31, 32, 33, 34, 35, 36, 37, 38, 39, 40, 41, 42, 43, 44, 45, 46, 47, 48, 49, 50, 51, 52, 53,
                  54, 55, 689, 690, 691, 692, 693, 694, 695, 696, 697, 698, 699, 700, 701, 702, 703, 704, 705, 706, 707,
                  708, 709, 878, 879, 880, 881, 882, 883, 884, 885, 886, 887, 888, 889, 890, 891, 892, 893, 894, 895,
                  896, 897, 898, 1466, 1467, 1468, 1469, 1470, 1471, 1472, 1473, 1474, 1475, 1476, 1477, 1478, 1479,
                  1480, 1481, 1482, 1483, 1484, 1485, 1486, 1487, 1488, 1489, 1490, 1491, 1492, 1493, 1494, 1495, 1496,
                  1497, 1498, 1499, 1500, 1501, 1502, 1503],
                2:[15, 16, 17, 18, 19, 20, 21, 22, 23, 24, 25, 26, 27, 28, 29, 30, 31, 32, 33, 34, 35, 36, 37, 38, 39,
                   40, 41, 42, 43, 44, 644, 645, 646, 647, 648, 649, 650, 651, 652, 653, 654, 655, 656, 657, 658, 659,
                   660, 661, 662, 663, 664, 800, 801, 802, 803, 804, 805, 806, 807, 808, 809, 810, 811, 812, 813, 814,
                   815, 816, 817, 818, 819, 820, 1330, 1331, 1332, 1333, 1334, 1335, 1336, 1337, 1338, 1339, 1340, 1341,
                   1342, 1343, 1344, 1345, 1346, 1347, 1348, 1349],
                3:[1, 2, 3, 4, 5, 6, 7, 8, 9, 10, 11, 12, 13, 14, 15, 16, 17, 18, 19, 20, 21, 22, 23, 24, 25, 26, 27,
                   28, 29, 30, 31, 503, 504, 505, 506, 507, 508, 509, 510, 511, 512, 513, 514, 515, 516, 517, 518, 519,
                   520, 521, 522, 698, 699, 700, 701, 702, 703, 704, 705, 706, 707, 708, 709, 710, 711, 712, 713, 714,
                   715, 716, 717, 1230, 1231, 1232, 1233, 1234, 1235, 1236, 1237, 1238, 1239, 1240, 1241, 1242, 1243,
                   1244, 1245, 1246, 1247],
                4:[1, 2, 3, 4, 5, 6, 7, 8, 9, 10, 11, 12, 13, 14, 15, 16, 17, 18, 19, 20, 21, 22, 23, 24, 25, 26, 27,
                   28, 29, 30, 31, 32, 33, 34, 35, 36, 37, 38],
                5:[1, 2, 3, 4, 5, 6, 7, 8, 9, 10, 11, 12, 13, 14, 15, 16, 17, 18, 19, 20, 21, 22, 23, 24, 25, 26, 27,
                   28, 29, 30, 31, 32, 33, 34, 35, 36, 37, 38, 39, 40, 41, 42, 43, 1005, 1006, 1007, 1008, 1009, 1010,
                   1011, 1012, 1013, 1014, 1015, 1016, 1017, 1018, 1019, 1020, 1021, 1022],
                6:[1, 2, 3, 4, 5, 6, 7, 8, 9, 10, 11, 12, 13, 14, 15, 16, 17, 18, 19, 20, 21, 22, 23, 24, 25, 26, 27,
                   28, 29, 30, 31, 32, 33, 34, 35, 36, 37, 38, 39, 40, 41, 42, 43, 44, 45, 46, 47, 48, 49, 50, 51, 52,
                   53, 54, 55],
                7:[1, 2, 3, 4, 5, 6, 7, 8, 9, 10, 11, 12, 13, 14, 15, 16, 17, 18, 19, 20, 21, 22, 23, 24, 25, 26, 27,
                   28, 29, 30, 31, 32, 33, 34, 35, 36, 37, 38, 39, 40, 41, 42],
                8:[], 9:[], 10:[1, 2, 3, 4, 5, 6, 7, 8, 9, 10, 11]}



def calculate_isnv_for_multiple_freqs_files(freqs_files, output_file, TiLV=False, freq_coutoff=0, read_count_cutoff=1000):
    """
    gets list of freqs file and calculated iSNV for each freqs file
    :param freqs_files: list of freqs files
    :param output_file: path of output file
    :param TiLV: if True - remove the position that are in TILV_FILTER
    :param freq_coutoff: freqs cut off
    :param read_count_cutoff: read_count cut off
    :return: saves the result df to output_file
    """
    output_file =check_filename(output_file, Truefile=False)
    result = pd.DataFrame(columns = ["segment", "sample", "experiment", "iSNV", "syn_iSNV", "non_syn_iSNV"])
    for f in freqs_files:
        df = pd.read_csv(f)
        sample = list(set(list(df["sample"])))[0]
        segment = list(set(list(df["segment"])))[0]
        experiment = list(set(list(df["experiment"])))[0]

        if TiLV:
            new_df = pd.DataFrame()
            new_df = new_df.append(df.loc[~df["Pos"].isin(TILV_FILTER[segment])])
            df = new_df

            if sample == "2":
                new_df = new_df.loc[new_df["Freq"]<0.01]
                df = new_df


        iSNV = calc_iSNV(df, read_count_cutoff=read_count_cutoff, freq_coutoff=freq_coutoff)
        syn_iSNV = calc_synonymous_iSNV(df, read_count_cutoff=read_count_cutoff, freq_coutoff=freq_coutoff)
        non_syn_iSNV = calc_non_synonymous_iSNV(df, read_count_cutoff=read_count_cutoff, freq_coutoff=freq_coutoff)

        result = result.append({"segment": segment,
                                        "sample": sample,
                                        "experiment": experiment,
                                        "iSNV": iSNV, "syn_iSNV":syn_iSNV, "non_syn_iSNV":non_syn_iSNV}, ignore_index=True)

    print(result)
    result.to_csv(output_file)

def calc_iSNV(df, read_count_cutoff=1000, freq_coutoff = 0):
    #calculate iSNV
    df = df.loc[df["Read_count"] >= read_count_cutoff]
    df = df.loc[df["Freq"] >= freq_coutoff]
    total_count = df["Read_count"].loc[df["Rank"]==0].sum()
    df.loc[:,"mutation_count"] = df["Freq"] * df["Read_count"]
    total_var = df["mutation_count"].loc[df["Rank"] !=0].sum()
    iSNV = total_var / total_count
    return iSNV


def calc_synonymous_iSNV(df, read_count_cutoff=1000, freq_coutoff = 0):
    #calculate synonymous iSNV
    df = df.loc[df["Read_count"] >= read_count_cutoff]
    df = df.loc[df["Freq"] >= freq_coutoff]
    total_count = df["Read_count"].loc[df["Rank"]==0].sum()
    df.loc[:,"mutation_count"] = df["Freq"] * df["Read_count"]
    synonymous = df.loc[df["Mutation_type"] == "synonymous"]
    total_var = synonymous["mutation_count"].loc[synonymous["Rank"] !=0].sum()
    iSNV = total_var / total_count
    return iSNV

def calc_non_synonymous_iSNV(df, read_count_cutoff=1000, freq_coutoff = 0):
    #calculate non-synonymous iSNV
    df = df.loc[df["Read_count"] >= read_count_cutoff]
    df = df.loc[df["Freq"] >= freq_coutoff]
    total_count = df["Read_count"].loc[df["Rank"]==0].sum()
    df.loc[:,"mutation_count"] = df["Freq"] * df["Read_count"]
    non_synonymous = df.loc[df["Mutation_type"] == "missense"]
    total_var = non_synonymous["mutation_count"].loc[non_synonymous["Rank"] !=0].sum()
    iSNV = total_var / total_count
    return iSNV