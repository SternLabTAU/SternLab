
import matplotlib
matplotlib.use('Agg')
import matplotlib.pyplot as plt
import seaborn as sns
import matplotlib.gridspec as gridspec
import os.path
import pathlib
from file_utilities import *
from optparse import OptionParser
from cirseq_utilities import *
sns.set_context("talk")
start_time = time.time()



def main():
    # for Cluster
    parser = OptionParser("usage: %prog [options]")
    parser.add_option("-f", "--freqs_file_path", dest="freqs_file_path", help="path of the freqs file")
    parser.add_option("-v", "--virus", dest="virus", help="Virus name: CVB3 for CV; RVB14 for RVB14, RVB6 for RVB6 ; PV for PV")
    parser.add_option("-c", "--sample", dest="sample_name", help="The name of the sample ie RV-p11")
    (options, args) = parser.parse_args()

    freqs_file = options.freqs_file_path
    virus = options.virus
    freqs_file = check_filename(freqs_file)
    sample = options.sample_name
    suffix = "%s.freqs" % sample

    #for Local
    # sample = "RV-p11"
    # suffix = "%s.freqs" % sample
    # freqs_file = "/sternadi/home/volume3/okushnir/AccuNGS/180503_OST_FINAL_03052018/merged/%s/q30_3UTR_new/%s" % \
    #              (sample, suffix)
    # virus = "RVB14"
    # seq_meth = "AccuNGS"



    # 1. Get freqs file and CirSeq running directory.
    path = freqs_file.split('/')[0:-1]
    out_dir = '' #the freqs directory
    for i in path:
        out_dir += str(i + '/')
    tmp_cirseq_dir = out_dir + 'tmp/'
    pathlib.Path(out_dir + 'plots/').mkdir(parents=True, exist_ok=True)
    out_plots_dir = out_dir + 'plots/'

    if virus == "CVB3":
        ncbi_id ="M16572"
    if virus == "RVB14":
        ncbi_id = "NC_001490"
    if virus == "RVB6.1":
        ncbi_id = "JQ747748"
    if virus == "RVB6.2":
        ncbi_id = "JN614996"
    if virus == "RVB6.3":
        ncbi_id = "JN562723"
    if virus == "RVB6.4":
        ncbi_id = "JF781502"
    if virus == "HIV":
        ncbi_id = "K03455"
    if virus == "PV":
        ncbi_id ="V01149"


    """ 2. Adding mutation types to the freqs file"""
    if not os.path.isfile(freqs_file[0:-5] + "with.mutation.type.freqs"):
         append_mutation = find_mutation_type(freqs_file, ncbi_id)
    freqs_file_mutations = freqs_file[0:-5] + "with.mutation.type.freqs"
    fig1 = plt.figure(figsize=(16, 9))
    ax = plt.subplot()
    make_boxplot_mutation(freqs_file_mutations, ax, out_plots_dir)
    plt.savefig(out_plots_dir + suffix.split(sep='.')[0] + '_All_Mutations.png', dpi=300)
    plt.close("all")
    print("The All Mutations Plot is ready in the folder")

    fig2 = plt.figure(figsize=(16, 9))
    ax = plt.subplot()
    make_boxplot_transition_mutation(freqs_file_mutations, ax, out_dir)
    plt.savefig(out_plots_dir + suffix.split(sep='.')[0] + '_Transitions_Report.png', dpi=300)
    plt.close("all")
    print("The Transition Plot is ready in the folder")


"""Graphs"""



#6. Mutation Rates
def make_boxplot_mutation(freqs_file, ax, output_dir):
    """
    Plots the mutation frequencies boxplot
    :param freqs_file: pandas DataFrame after find_mutation_type function
    :param ax: ax location
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
        # data['Frequency'] = data['abs_counts'].apply(lambda x: 1 if x == 0 else x) / data["Read_count"]
        # raise Exception("This script does not support freqs file with deletions, for support please contact Maoz ;)"
    data['Base'].replace('T', 'U', inplace=True)
    data['Ref'].replace('T', 'U', inplace=True)
    min_read_count = 1000
    data = data[data['Pos'] == np.round(data['Pos'])]  # remove insertions
    data['Pos'] = data[['Pos']].apply(pd.to_numeric)
    data = data[data['Read_count'] > min_read_count]
    data['mutation_type'] = data['Ref'] + data['Base']
    data = data[data['Ref'] != data['Base']]
    data = data[data["Base"] != "-"]
    data['abs_counts'] = data['Freq'] * data["Read_count"]  # .apply(lambda x: abs(math.log(x,10))/3.45)
    data['Frequency'] = data['abs_counts'].apply(lambda x: 1 if x == 0 else x) / data["Read_count"]
    data["Mutation"] = data["Ref"] + "->" + data["Base"]
    data.to_csv(output_dir + "data_all_mutation.csv", sep=',', encoding='utf-8')
    sns.set_palette(sns.color_palette("Paired", 12))
    g1 = sns.boxplot(x="Mutation Type", y="Frequency", hue="Mutation", data=data,
                     hue_order=["C->U", "U->C", "G->A", "A->G", "C->A", "G->U", "U->G", "U->A", "G->C", "A->C",
                                "A->U", "C->G"], order=["Synonymous", "Non-Synonymous", "Premature Stop Codon"], ax=ax)
    g1.set(yscale="log")
    plt.legend(bbox_to_anchor=(1.0, 1), loc=2, borderaxespad=0., fontsize=6)
    g1.set_ylim(10 ** -6, 1)
    g1.tick_params(labelsize=7)

def make_boxplot_transition_mutation(freqs_file,ax, output_dir):
    """
    Plots the mutation frequencies boxplot
    :param freqs_file: pandas DataFrame after find_mutation_type function
    :param ax: ax location
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
        # data['Frequency'] = data['abs_counts'].apply(lambda x: 1 if x == 0 else x) / data["Read_count"]
        # raise Exception("This script does not support freqs file with deletions, for support please contact Maoz ;)"
    data['Base'].replace('T', 'U', inplace=True)
    data['Ref'].replace('T', 'U', inplace=True)
    min_read_count = 1000
    data = data[data['Pos'] == np.round(data['Pos'])]  # remove insertions
    data['Pos'] = data[['Pos']].apply(pd.to_numeric)
    data = data[data['Read_count'] > min_read_count]
    data['mutation_type'] = data['Ref'] + data['Base']
    data = data[data['Ref'] != data['Base']]
    data = data[data["Base"] != "-"]
    data['abs_counts'] = data['Freq'] * data["Read_count"]  # .apply(lambda x: abs(math.log(x,10))/3.45)
    data['Frequency'] = data['abs_counts'].apply(lambda x: 1 if x == 0 else x) / data["Read_count"]
    data["Mutation"] = data["Ref"] + "->" + data["Base"]

    data.to_csv(output_dir + "data_transitions_mutation.csv", sep=',', encoding='utf-8')

    sns.set_palette(sns.color_palette("Paired", 12))

    g1 = sns.factorplot(x="Mutation Type", y="Frequency", data=data, col="Mutation",
                     col_order=["C->U", "U->C", "G->A", "A->G"], kind="box")
    g1.set_xticklabels(["Synonymous", "Non\nSynonymous", "Premature\nStop Codon"], fontsize=10)
    g1.set_xlabels('')
    g1.set(yscale="log", ylim=(0.000001, 0.01))


if __name__ == "__main__":
    main()
