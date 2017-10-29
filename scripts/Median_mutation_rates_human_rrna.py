"""
@Author: odedkushnir

"""

import pandas as pd
import numpy as np
import matplotlib
matplotlib.use('Agg')
import matplotlib.pyplot as plt
import seaborn as sns
import pathlib


def main():
    # For Local Run
    suffix = "CVB3-p2.freqs"
    dir_path = "/Volumes/STERNADILABTEMP$/volume1/okushnir/Cirseq/CV/20170802_q30r2_merged_by_id/"
    plots_path = dir_path + "/plots/"
    freqs_path = dir_path + suffix
    pathlib.Path(plots_path).mkdir(parents=True, exist_ok=True)

    make_boxplot_mutation_median(freqs_path, plots_path, suffix[:-6])


def make_boxplot_mutation_median(freqs_file, out_dir, samplename):
    """
    :param freqs_file: freqs file from the cirseq pipeline
    :param out_dir: The output directory
    :param samplename: The objective name
    :return:
    """
    data = pd.read_table(freqs_file)
    data.reset_index(drop=True, inplace=True)
    flag = '-' in data.Base.values
    if flag is True:
        data = data[data.Ref != '-']
        data = data[data.Base != '-']
        data.reset_index(drop=True, inplace=True)
        # data['abs_counts'] = data['Freq'] * data["Read_count"]  # .apply(lambda x: abs(math.log(x,10))/3.45)
        # data['Frequency'] = data['abs_counts'] / data["Read_count"]
        # raise Exception("This script does not support freqs file with deletions, for support please contact Maoz ;)")
    data['Base'].replace('T', 'U', inplace=True)
    data['Ref'].replace('T', 'U', inplace=True)
    min_read_count = 100000
    data = data[data['Pos'] == np.round(data['Pos'])]  # remove insertions
    data['Pos'] = data[['Pos']].apply(pd.to_numeric)
    data = data[data['Read_count'] > min_read_count]
    data['mutation_type'] = data['Ref'] + data['Base']
    data = data[data['Ref'] != data['Base']]
    data["Mutation"] = data["Ref"] + "->" + data["Base"]
    sns.set_palette(sns.color_palette("Paired", 12))
    """ To Sample from the data"""
    # data = data[data['Pos'] > 1000]
    # data = data[data['Pos'] < 2800]

    g1 = sns.boxplot(x="Mutation", y="Freq", data=data,
                 order=["A->C", "A->G", "A->U", "C->A", "C->G", "C->U", "G->A", "G->C", "G->U", "U->A",
                        "U->C", "U->G"], color="DarkOrchid")
    grouped_df = data.groupby(["Mutation"])["Freq"]
    len_col = grouped_df.size()
    medians = grouped_df.median()
    print(medians)
    median_labels = [str(np.round(s, 6)) for s in medians]
    pos = range(len(medians))
    for tick, label in zip(pos, g1.get_xticklabels()):
        g1.text(pos[tick], medians[tick]*1.01, median_labels[tick],
                horizontalalignment='center', size='x-small', color='w', weight='semibold', fontsize=5)
    for tick, label in zip(pos, g1.get_xticklabels()):
        g1.text(pos[tick], medians[tick]+0.5, len_col[tick],
                horizontalalignment='center', size='x-small', color='navy', weight='semibold', fontsize=5)
    sns.set_style("darkgrid")
    g1.set(yscale="log")
    plt.legend(bbox_to_anchor=(1.0, 1), loc=2, borderaxespad=0., fontsize=7)
    g1.set_ylim(10 ** -6, 1)
    g1.tick_params(labelsize=7)
    plt.title(samplename + ' Mutations Frequencies', fontsize=19)
    plt.savefig(out_dir + samplename + "_median_%s.png" % str(min_read_count), dpi=300)
    plt.close("all")
    print("The Plot is ready in the folder")


if __name__ == "__main__":
    main()










