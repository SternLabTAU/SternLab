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
    suffix = "DAVIDp6_rrna_q30r3.freqs"
    dir_path = "//Users/odedkushnir/Google Drive/Studies/PhD/Projects/HRVB14/Maoz/"
    plots_path = dir_path + "/plots/"
    freqs_path = dir_path + suffix
    mutation_file = freqs_path
    mutation_rates = freqs_to_dataframe(mutation_file)
    pathlib.Path(plots_path).mkdir(parents=True, exist_ok=True)
    make_boxplot_mutation_median(mutation_rates, plots_path, suffix)


def freqs_to_dataframe(freqs_file):
    """
    This function returns arranged DataFrame without deletions
    :param freqs_file: The path of the relevant freqs file
    :return: DataFrame without deletions
    """
    data = pd.read_table(freqs_file)
    data = data[data.Ref != '-']
    data = data[data.Base != '-']
    data.reset_index(drop=True, inplace=True)
    return data


def make_boxplot_mutation_median(data, out_dir, filename):
        """
            :param data: pandas DataFrame after find_mutation_type function
            :return: pandas DataFrame ready for plotting
            """
        data['Mutation Type'] = 'Non-Synonymous'
        data['Base'].replace('T', 'U', inplace=True)
        data['Ref'].replace('T', 'U', inplace=True)
        min_read_count = 100000
        data = data[data['Pos'] == np.round(data['Pos'])]  # remove insertions
        data['Pos'] = data[['Pos']].apply(pd.to_numeric)
        data = data[data['Read_count'] > min_read_count]
        data['mutation_type'] = data['Ref'] + data['Base']
        data = data[data['Ref'] != data['Base']]
        data = data[data["Base"] != "-"]
        data["Mutation"] = data["Ref"] + "->" + data["Base"]
        sns.set_palette(sns.color_palette("Paired", 12))
        non_syn = data[data['Mutation Type'] == 'Non-Synonymous']
        # non_syn = non_syn[non_syn['Pos'] > 1000]
        # non_syn = non_syn[non_syn['Pos'] < 2800]
        g1 = sns.boxplot(x="Mutation", y="Freq", data=non_syn,
                         order=["A->C", "A->G", "A->U", "C->A", "C->G", "C->U", "G->A", "G->C", "G->U", "U->A",
                                "U->C", "U->G"], color="DarkOrchid")
        sns.set_style("darkgrid")
        grouped_df = non_syn.groupby(["Mutation"])["Freq"]
        for key, item in grouped_df:
            print(grouped_df.get_group(key), "\n\n")

        medians = grouped_df.median()
        print(medians)
        median_labels = [str(np.round(s, 6)) for s in medians]
        pos = range(len(medians))
        for tick, label in zip(pos, g1.get_xticklabels()):
            g1.text(pos[tick], medians[tick]*1.01, median_labels[tick],
                    horizontalalignment='center', size='x-small', color='w', weight='semibold', fontsize=5)

        len_col = grouped_df.size()
        for tick, label in zip(pos, g1.get_xticklabels()):
            g1.text(pos[tick], medians[tick]+0.5, len_col[tick],
                    horizontalalignment='center', size='x-small', color='navy', weight='semibold', fontsize=5)
        g1.set(yscale="log")
        plt.legend(bbox_to_anchor=(1.0, 1), loc=2, borderaxespad=0., fontsize=7)
        g1.set_ylim(10 ** -6, 1)
        g1.tick_params(labelsize=7)
        plt.title(filename + 'Mutations Frequencies', fontsize=19)
        plt.savefig(out_dir + filename[:-6] + "_median_%s.png" % str(min_read_count), dpi=300)
        plt.close("all")
        print("The Plot is ready in the folder")


if __name__ == "__main__":
    main()










