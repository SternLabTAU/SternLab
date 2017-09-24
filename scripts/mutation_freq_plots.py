
import matplotlib
matplotlib.use('Agg')
import matplotlib.pyplot as plt
import seaborn as sns
import os.path
import pathlib
from optparse import OptionParser
from cirseq_utilities import *



def main():
    # For Local Run
    suffix = "CVB3-p2.freqs"
    path = "/Volumes/STERNADILABTEMP$/volume1/okushnir/Cirseq/CV/20170906_q20r3_blastn/" + suffix
    virus = "CVB3"

    # file_name = "PV-p3.1036617.freqs"
    # virus = file_name.split(sep='-')[0]
    # virus += file_name.split(sep='-')[1].split(sep='.')[0]

    # For Cluster Run
    # parser = OptionParser("usage: %prog [options]")
    # parser.add_option("-f", "--freqs_file_path", dest="freqs_file_path", help="path of the freqs file")
    # parser.add_option("-v", "--virus", dest="virus", help="Virus name: CVB3 for CV; RVB14 for RV; PV for PV")
    # (options, args) = parser.parse_args()
    # path = options.freqs_file_path
    # virus = options.virus

    path_list = path.split(sep="-")[0].split(sep="/")
    plots_path = ''
    for i in range(0, len(path_list) - 1, 1):
        plots_path += path_list[i] + "/"



    if virus == "CVB3":
        ncbi_id ="M16572"
    if virus == "RVB14":
        ncbi_id = "NC_001490"
    if virus == "PV":
        ncbi_id ="V01149"

    if not os.path.isfile(path + ".with.mutation.type.func2.freqs"):
        append_mutation = find_mutation_type(path, ncbi_id)


    mutation_file = path + ".with.mutation.type.func2.freqs"
    mutation_rates = freqs_to_dataframe(mutation_file)

    pathlib.Path(plots_path + 'plots/').mkdir(parents=True, exist_ok=True)

    df = make_boxplot_mutation(mutation_rates, plots_path, suffix[:-6])
    make_boxplot_mutation_median(df, plots_path, suffix[:-6], "Premature Stop Codon")


# from Maoz
def make_boxplot_mutation(data, out_dir, samplename):
    """
    :param data: pandas DataFrame after find_mutation_type function
    :param out_dir:  The output directory
    :param samplename: The objective name
    :return: pandas DataFrame ready for plotting
    """

    data['Base'].replace('T', 'U', inplace=True)
    data['Ref'].replace('T', 'U', inplace=True)
    min_read_count = 100000
    data = data[data['Pos'] == np.round(data['Pos'])]  # remove insertions
    # data['Pos'] = data[['Pos']].apply(pd.to_numeric)
    data = data[data['Read_count'] > min_read_count]
    data['mutation_type'] = data['Ref'] + data['Base']
    data = data[data['Ref'] != data['Base']]
    data = data[data["Base"] != "-"]
    data['abs_counts'] = np.round(data['Freq'] * data["Read_count"])
    data['Frequency'] = data['abs_counts'].apply(lambda x: 1 if x == 0 else x) / data["Read_count"]
    data["Mutation"] = data["Ref"] + "->" + data["Base"]
    sns.set_palette(sns.color_palette("Paired", 12))
    g1 = sns.boxplot(x="Mutation Type", y="Frequency", hue="Mutation", data=data,
                     hue_order=["C->U", "U->C", "G->A", "A->G", "C->A", "G->U", "U->G", "U->A", "G->C", "A->C",
                                "A->U", "C->G"], order=["Synonymous", "Non-Synonymous", "Premature Stop Codon"])
    g1.set(yscale="log")
    sns.set_style("darkgrid")
    plt.legend(bbox_to_anchor=(1.0, 1), loc=2, borderaxespad=0., fontsize=7)
    g1.set_ylim(10 ** -6, 1)
    g1.tick_params(labelsize=7)
    plt.title(samplename + ' Mutations Frequencies', fontsize=22)
    plt.savefig(out_dir + "plots/" + samplename + "_freqs_type_%s_with_correction.png" % str(min_read_count), dpi=300)
    plt.close("all")
    print("The Plot is ready in the folder")
    return data


def make_boxplot_mutation_median(data, out_dir, samplename, mutation_type):
        """
            :param data: pandas DataFrame after find_mutation_type function
            :param out_dir:
            :param virus:
            :param samplename:
            :param mutation_type: "Synonymous"/"Non-Synonymous"/"Premature Stop Codon"
            :return: pandas DataFrame ready for plotting
            """
        min_read_count = 100000
        mutation_df = data[data['Mutation Type'] == mutation_type]
        # mutation_df = mutation_df[mutation_df['Pos'] > 1000]
        # mutation_df = mutation_df[mutation_df['Pos'] < 2800]
        mutation_df['abs_counts'] = np.round(mutation_df['Freq'] * mutation_df["Read_count"])  # .apply(lambda x: abs(math.log(x,10))/3.45)
        mutation_df['Frequency'] = mutation_df['abs_counts'].apply(lambda x: 1 if x == 0 else x) / mutation_df["Read_count"]
        g1 = sns.boxplot(x="Mutation", y="Frequency", data=mutation_df,
                         order=["A->U", "C->A", "C->G", "C->U", "G->A", "G->U", "U->A", "U->G"], color="DarkOrchid")
        sns.set_style("darkgrid")
        grouped_df = mutation_df.groupby(["Mutation"])["Freq"]
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
        plt.title(samplename + mutation_type +' Mutations Frequencies', fontsize=19)
        plt.savefig(out_dir + "plots/" + samplename + "_" + mutation_type + "_median_%s_with_correction.png" % str(min_read_count), dpi=300)
        plt.close("all")
        print("The Plot is ready in the folder")


if __name__ == "__main__":
    main()










