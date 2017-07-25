
import matplotlib.pyplot as plt
import pandas as pd
import seaborn as sns
import numpy as np
import re
import glob
from Bio.Seq import Seq
import matplotlib.gridspec as gridspec
import time
import os.path
import pathlib
from Bio import Entrez
from Bio import SeqIO
from optparse import OptionParser
sns.set_context("talk")
start_time = time.time()



# ### Pipeline:
#           1. Get freqs file and CirSeq running directory.
#           2. Analyze those file and directory to get the number of tandem repeats of the cirseq, repeats length and the
#                 amount of reads per repeat
#           3. Adding mutation types to the freqs file
#           4. Run bowtie2 for human rRNA, mRNA and the virus
#           5. Subplots all relevant graph in subplots
#                 5.1. Distribution graph (=bowtie2 results)
#                 5.2. Multiple tandem repeat graph (repeat len)
#                 5.3. Reads length
#                 5.4. Amount of reads per repeat (Reads VS. Repeats Graph)
#                 5.5. Coverage
#                 5.6 Plot virus mutation frequencies(rates)

#Functions
def get_freqs_file(freqs_file):
    """
    :param freqs_file: The freqs file
    :return: pandas DataFrame of the freqs file
    """
    with open(freqs_file, 'r') as freqs_file:
        df = pd.read_csv(freqs_file, sep='\t')
        return df


def get_repeats_num(tmp_cirseq_dir, out_dir):
    files = glob.glob(tmp_cirseq_dir + "/*.fasta.blast.freqs.stats")
    repeat_summery = {}
    for file in files:
        pattern = re.compile("(\d+\t{1}\d+\n{1})", re.MULTILINE)
        text = open(file, "r").read()
        reads = pattern.findall(text)
        for r in reads:
            key = int(r.split("\t")[0])
            value = int(r.split("\t")[1].split("\n")[0])
            if key in repeat_summery:
                repeat_summery[key] += value
            else:
                repeat_summery[key] = value
    np.save(out_dir + 'repeat_summery.npy', repeat_summery)
    print("repeat_summery.npy is in your directory")
    return repeat_summery


def get_read_and_repeat_length(in_dir, out_dir):
    files = glob.glob(in_dir + "/*.fasta.blast")
    arr_names = ["sseqid", "qstart", "qend", "sstart", "send", "sstrand", "length", "btop"]
    wdf = pd.DataFrame()
    for file in files:
        data = pd.read_csv(file, delimiter="\t",  header=None, names=arr_names)
        grouped_and_filtered = data.groupby("sseqid").filter(lambda x: min(x["qend"]) - max(x["qstart"]) > 0).groupby("sseqid")["length"].agg([np.count_nonzero, np.sum, np.max])
    wdf = pd.DataFrame.append(wdf, grouped_and_filtered)
    wdf.to_csv(out_dir + '/length_df.csv', sep='\t', encoding='utf-8')
    print("length_df.csv is in your folder")
    return wdf


def parse_reads(freqs):
    """
    this method returns a vector of reads corresponding to genome positions.
    :param freqs: freqs file
    :return: an integer vector containing for each position in the genome it's num of reads.
    """

    path = freqs
    df = pd.read_csv(path, sep='\t')

    # remove all duplicates from Pos except the first occurrence
    # remove all x.number duplicates
    df[["Pos"]] = df[["Pos"]].astype(int)
    df = df.drop_duplicates("Pos")

    pos = df["Pos"]  # a vector of all positions
    reads = df["Read_count"]
    return pos, reads


def find_mutation_type(freqs_file, ncbi_id):
    """
    This function adds Mutation type to the freqs file
    :param freqs_file:  The path of the relevant freqs file
    :return:DataFrame of the frqs file with mutation type column, save it in txt file
    """
    file_name = freqs_file
    data = freqs_to_dataframe(freqs_file)
    start_pos, end_pos = find_coding_region(ncbi_id)
    data = data.loc[data['Pos'] >= start_pos]
    check_12mer(data)
    data['Codon'] = ""
    data['Ref_Protein'] = ""
    data['Potential_Protein'] = ""
    data['Mutation Type'] = ""
    data['Pos'] = data['Pos'].astype(int)
    data.reset_index(drop=True, inplace=True)
    for kmer in range(data.first_valid_index(), (len(data)), 12):
        print ("going over kmer: " + str(kmer) + "/" + str((len(data))))
        kmer_df = data[:].loc[kmer:kmer+11]
        print("finding codons....")
        find_codon(kmer_df)
        print("translating codons....")
        translate_codon(kmer_df)
        print("Sets the Mutation type....")
        kmer_df['Mutation Type'] = kmer_df[['Ref_Protein', 'Potential_Protein']].apply(
            lambda protein: check_mutation_type(protein[0], protein[1]), axis=1)
    print("After a long for loop")
    file_name += ".with.mutation.type.func2.freqs"
    data.to_csv(file_name, sep='\t', encoding='utf-8')
    print("The File is ready in the folder")
    print("--- %s sec ---" % (time.time() - start_time))
    return data


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


def find_coding_region(ncbi_id):
    """
    :param ncbi_id:
    :return: start and end Pos of CDS
    """
    try:
        Entrez.email = "A.N.Other@example.com"
        handle = Entrez.efetch(db="nucleotide", rettype="gb", retmode="text", id=ncbi_id)
        ncbi_gb = SeqIO.read(handle, "gb")
        handle.close()
        location = list(ncbi_gb.features[1].location)
        start_pos = location[0]
        end_pos = location[-4]
        return start_pos, end_pos
    except:
        print('Failed to fetch record.')


def check_12mer(data):
    """
      :param data: pandas DataFrame of the freqs file
      :return: None - checks that DataFrame is in 12 kmer
      """
    if (len(data)) % 12 != 0:
        print('The data is not in 12 mer. arrange the data')
    else:
        print('The data is in 12 mer. All OK.')
        print('The length of data is:', (len(data)))


def find_codon(data):
    """
    :param data: pandas DataFrame of the freqs file
    :return: pandas DataFrame with codon column
    """
    data['Codon'] = ""
    x = 4
    first_index_kmer = data.first_valid_index()
    for i in range(first_index_kmer, first_index_kmer+x, 1):
        data['Codon'][i] = data["Base"].loc[i] + data["Ref"].loc[first_index_kmer+x] + data["Ref"].loc[(first_index_kmer+(x*2))]
    for i in range((first_index_kmer+x), (first_index_kmer+(x*2)), 1):
            data['Codon'][i] = data["Ref"].loc[first_index_kmer] + data["Base"].loc[i] + data["Ref"].loc[first_index_kmer+(x*2)]
    for i in range((first_index_kmer+(x*2)), (first_index_kmer+(x * 3)), 1):
            data['Codon'][i] = data["Ref"].loc[first_index_kmer] + data["Ref"].loc[first_index_kmer+x] + data["Base"].loc[i]
    return data


def translate(seq):
    """
    :param seq: string of ATGC
    :return: protein sequence of seq
    """
    seq = Seq(str(seq))
    protein = seq.translate()
    return protein


def translate_codon(data):
    """
    :param data: pandas DataFrame of the freqs file
    :return: pandas DataFrame with Reference protein column and Potential protein column
    """
    data['Ref_Protein'] = ""
    data['Potential_Protein'] = ""
    protein_ref = Seq(data["Codon"][data.first_valid_index()]).translate()
    data['Ref_Protein'].loc[data.first_valid_index():data.first_valid_index()+11] = protein_ref
    data['Ref_Protein'].replace('(', '', inplace=True)
    data['Ref_Protein'].replace(')', '', inplace=True)

    data["Potential_Protein"] = data["Codon"].apply(translate)
    data["Potential_Protein"].replace('(', '', inplace=True)
    data["Potential_Protein"].replace(')', '', inplace=True)
    return data


def check_mutation_type(protein1, protein2):
    """
    :param protein1: amino acid 1
    :param protein2: amino acid 2
    :return: The mutation type
    """
    Mutation_Type = ""
    if protein1 == protein2:
        Mutation_Type = "Synonymous"
    elif protein1 != protein2:
        Mutation_Type = "Non-Synonymous"
    if protein2 == '*':
        Mutation_Type = "Premature Stop Codon"
    return Mutation_Type

"""Graphs"""
#1.Distribution graph (=bowtie2 results)


def distribution_graph(val, ax, virus):
    """
    makes distribution graph
    :param val: values from bowtie2 alignment
    :param ax:
    :return:
    """
    column_number = 3

    ind = np.arange(column_number)  # the x locations for the groups
    width = 0.35  # the width of the bars
    rects1 = ax.bar(ind, val, width, color='DarkOrchid')
    ax.set_ylabel('% of Reads')
    ax.set_xticks(ind)
    ax.set_xticklabels((virus, 'mRNA', 'rRNA'))
    labelsy = np.arange(0, 50, 10)
    ax.set_yticklabels(labelsy)
    sns.set_style("darkgrid")
    ax.set_xlim(-0.5, 3)
    return ax


#2. Multiple tandem repeat graph (repeat len)
def repeat_len_graphs(dataframe, ax):
    """
     makes a repeat length graph
    :param results:
    :param ax:
    :return:
    """
    sns.boxplot(x="count_nonzero", y="amax", data=dataframe, ax=ax)
    ax.set_xlabel("Number of Repeats")
    ax.set_ylabel("Repeat Length [bp]")
    labelsy = np.arange(0, 350, 50)
    labelsx = np.arange(1, 11, 1)
    ax.set_yticklabels(labelsy)
    ax.set_xticklabels(labelsx)
    sns.set_style("darkgrid")
    return ax

#3. Reads length
def read_len_graphs(dataframe, ax):

    graph = sns.boxplot(x="count_nonzero", y="sum", data=dataframe, ax=ax)
    ax.set_xlabel("Number of Repeats")
    ax.set_ylabel("Read Length [bp]")
    sns.set_style("darkgrid")


#4. Amount of reads per repeat (Reads VS. Repeats Graph)
def read_repeat_graph(repeat_summery, ax):
    """
    makes read VS. repeat graph
    :param repeat_summery:
    :param ax:
    :return:
    """
    keys = []
    values = []
    for key, val in repeat_summery.items():
        keys.append(key)
        values.append(val)
    graph = ax.bar(keys, values, color='DarkOrchid')
    ax.set_xlabel("Number of Repeats")
    ax.set_ylabel("Number of Reads")
    ax.set_xticklabels(list(range(1, 11)))
    ax.set_xlim(min(keys)-0.5, max(keys)+0.5)
    ax.yaxis.get_major_formatter().set_powerlimits((0, 1))
    sns.set_style("darkgrid")
    return ax


#5. Coverage
def coverage_graph(freqs, ax):
    data = parse_reads(freqs)
    pos = data[0]
    reads = data[1]
    ax.plot(pos, reads)
    ax.set_xlabel("Position In The Genome [bp]")
    ax.set_ylabel("Number Of Reads")
    sns.set_style("darkgrid")
    ax.set_xlim(0, (len(pos)+10))
    ax.set_ylim(1000, 1000000)
    ax.set_yscale("log")


#6. Mutation Rates
def make_boxplot_mutation(data, ax):
    """
    Plots the mutation frequencies boxplot
    :param data: pandas DataFrame after find_mutation_type function
    :param ax: which ax to plot
    :return:
    """
    data['Base'].replace('T', 'U', inplace=True)
    data['Ref'].replace('T', 'U', inplace=True)
    min_read_count = 100000
    data = data[data['Pos'] == np.round(data['Pos'])]  # remove insertions
    data['Pos'] = data[['Pos']].apply(pd.to_numeric)
    data = data[data['Read_count'] > min_read_count]
    data['mutation_type'] = data['Ref'] + data['Base']
    data = data[data['Ref'] != data['Base']]
    data = data[data["Base"] != "-"]
    data['abs_counts'] = data['Freq'] * data["Read_count"]  # .apply(lambda x: abs(math.log(x,10))/3.45)
    data['Frequency'] = data['abs_counts'].apply(lambda x: 1 if x == 0 else x) / data["Read_count"]
    data["Mutation"] = data["Ref"] + "->" + data["Base"]
    sns.set_palette(sns.color_palette("Paired", 12))
    g1 = sns.boxplot(x="Mutation Type", y="Frequency", hue="Mutation", data=data,
                     hue_order=["C->U", "U->C", "G->A", "A->G", "C->A", "G->U", "U->G", "U->A", "G->C", "A->C",
                                "A->U", "C->G"], order=["Synonymous", "Non-Synonymous", "Premature Stop Codon"], ax=ax)
    g1.set(yscale="log")
    plt.legend(bbox_to_anchor=(1.0, 1), loc=2, borderaxespad=0.)
    g1.set_ylim(10 ** -6, 1)
    g1.tick_params(labelsize=7)


def main():
    # for Cluster
    parser = OptionParser("usage: %prog [options]")
    parser.add_option("-f", "--freqs_file_path", dest="freqs_file_path", help="path of the freqs file")
    parser.add_option("-v", "--virus", dest="virus", help="Virus name: CVB3 for CV; RVB14 for RV; PV for PV")
    (options, args) = parser.parse_args()

    freqs_file = options.freqs_file_path
    virus = options.virus

    #for Local

    # freqs_file = 'C:/Users/Oded/Google Drive/Studies/PhD/test/CVB3-p2.freqs'
    # virus = "CVB3"



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
    if virus == "PV":
        ncbi_id ="V01149"
    file_name = freqs_file.split('/')[-1]
    virus_name = file_name.split('-')[0]
    virus_name += file_name.split(sep='-')[1].split(sep='.')[0]



    """2. Analyze those file and directory to get the number of tandem repeats of the cirseq,
    repeats length and the amount of reads per repeat"""

    if os.path.isfile(out_dir + '/repeat_summery.npy'):
            repeats_dict = np.load(out_dir + '/repeat_summery.npy').item() #dict for read vs. repeat  , encoding='latin1'
    else:
            print("Counting the repeats number")
            repeats_dict = get_repeats_num(tmp_cirseq_dir, out_dir) #dict for read vs. repeat

    if os.path.isfile(out_dir + '/length_df.csv'):
            read_and_repeat_length = pd.DataFrame.from_csv(out_dir + '/length_df.csv', sep='\t', encoding='utf-8') #pandas_df for read and repeat length
    else:
            print("Getting the read and repeat length")
            read_and_repeat_length = get_read_and_repeat_length(tmp_cirseq_dir, out_dir) #pandas_df for read and repeat length


    """ 3. Adding mutation types to the freqs file"""
    if not os.path.isfile(freqs_file + ".with.mutation.type.func2.freqs"):
         append_mutation = find_mutation_type(freqs_file, ncbi_id)

    mutation_file = freqs_file + ".with.mutation.type.func2.freqs"
    mutation_rates = freqs_to_dataframe(mutation_file)


    """4. Run bowtie2 for human rRNA, mRNA and the virus"""


    # RV
    # values = (19.32, 40.52, 45.17)
    # CV
    values = (96.09, 1.66, 1.18)

    #   5. Subplots all relevant graph in subplots
    #      5.1. Distribution graph (=bowtie2 results)
    #      5.2. Multiple tandem repeat graph (repeat len)
    #      5.3. Reads length
    #      5.4. Amount of reads per repeat (Reads VS. Repeats Graph)
    #      5.5. Coverage
    #      5.6 Plot virus mutation frequencies(rates)

    plt.close('all')
    fig = plt.figure()
    gs = gridspec.GridSpec(3, 3)
    ax0 = plt.subplot(gs[0, 0])
    ax1 = plt.subplot(gs[0, 1])
    ax2 = plt.subplot(gs[0, 2])
    ax3 = plt.subplot(gs[1, 0])
    ax4 = plt.subplot(gs[1, 1:])
    ax5 = plt.subplot(gs[2:, :])
    gs.tight_layout(fig)

    fig.suptitle(virus_name + ' Analysis', fontsize=20)
    distribution_graph(values, ax0, virus)
    # ax0.set_title('Reads Distribution')

    repeat_len_graphs(read_and_repeat_length, ax1)
    # ax1.set_title('Multiple Tandem Repeat Degree')

    read_len_graphs(read_and_repeat_length, ax2)
    # ax2.set_title('Read length')

    read_repeat_graph(repeats_dict, ax3)
    # ax3.set_title('Amount of reads per repeat')

    coverage_graph(freqs_file, ax4)
    # ax4.set_title('Coverage')

    make_boxplot_mutation(mutation_rates, ax5)
    # ax5.set_title(virus_name + ' Mutation Rates')

    # plt.show()
    plt.savefig(out_plots_dir + 'all.png', dpi=300)
    plt.close("all")
    print("The Plot is ready in the folder")

if __name__ == "__main__":
    main()
