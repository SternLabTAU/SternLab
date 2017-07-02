

# get_ipython().magic('matplotlib inline')
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
import Final_project
sns.set_context("talk")
start_time = time.time()



# ### Pipeline:
#           1. Get freqs file and CirSeq running directory.
#           2. Analyze those file and directory to get the number of tandem repeats of the cirseq, repeats length and the
#                 amount of reads per repeat
#           3. Get the coverage of the CirSeq run
#           4. Adding mutation types to the freqs file
#           5. Run bowtie2 for human rRNA, mRNA and the virus
#           6. Subplots all relevant graph in subplots
#                 6.1. Distribution graph (=bowtie2 results)
#                 6.2. Multiple tandem repeat graph (repeat len)
#                 6.3. Reads length
#                 6.4. Amount of reads per repeat (Reads VS. Repeats Graph)
#                 6.5. Coverage
#                 6.6 Plot virus mutation frequencies(rates)

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
    np.save(out_dir + '/repeat_summery.npy', repeat_summery)
    return repeat_summery


def getoverlap(a, b):
    return max(0, min(a[1], b[1]) - max(a[0], b[0]))


def get_read_and_repeat_length(tmp_cirseq_dir, out_dir):
    results = {}
    files = glob.glob(tmp_cirseq_dir + "/*.fasta.blast")
    for file in files:

        data = pd.read_csv(file, delimiter="\t",  header=None, names=["sseqid", "qstart", "qend", "sstart", "send", "sstrand", "length", "btop"])
        unique_sseqid = data.sseqid.unique() #get all unique reads

        for sseqid in unique_sseqid:
            sseqid_data = data.loc[data.sseqid == sseqid] #get the repeats of this read
            repeat_num = len(sseqid_data)
            # read with only one repeat
            if repeat_num == 1:
                read_length = sseqid_data["length"].tolist()[0]
                if not repeat_num in results:
                    results[repeat_num] = {}
                    results[repeat_num]["read_length"] = []
                    results[repeat_num]["repeat_length"] = []
                results[repeat_num]["read_length"].append(read_length)
                results[repeat_num]["repeat_length"].append(read_length)
                continue
            #check if read is good
            start = sseqid_data["qstart"].tolist()
            end = sseqid_data["qend"].tolist()
            overlap = True
            for i in range(0, len(start)-1):
                overlap = getoverlap([start[i], end[i]], [start[i+1], end[i+1]])
                if overlap == 0:
                    overlap = False

            #calculate read length and repeat lengths
            repeat_lengths = sseqid_data["length"].tolist()
            read_length = sum(repeat_lengths)

            #deal with a bad read
            if not overlap: # overlap==False
                if not "bad" in results:
                    results["bad"] = {}
                    results["bad"] = {}
                if not repeat_num in results["bad"]:
                    results["bad"][repeat_num] = {}
                    results["bad"][repeat_num]["read_length"] = []
                    results["bad"][repeat_num]["repeat_length"] = []
                results["bad"][repeat_num]["read_length"].append(read_length)
                results["bad"][repeat_num]["repeat_length"] += repeat_lengths
            else:
                if not repeat_num in results:
                    results[repeat_num] = {}
                    results[repeat_num]["read_length"] = []
                    results[repeat_num]["repeat_length"] = []
                results[repeat_num]["read_length"].append(read_length)
                results[repeat_num]["repeat_length"] += repeat_lengths
    print("After a long for-loop")
    np.save(out_dir + '/results.npy', results)
    return results


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


def find_mutation_type(freqs_file):
    file_name = freqs_file
    data = freqs_to_dataframe(freqs_file)
    check_12mer(data)
    data['Codon'] = ""
    data['Ref_Protein'] = ""
    data['Potential_Protein'] = ""
    data['Mutation Type'] = ""

    for kmer in range(data.first_valid_index(), len(data), 12):
        print ("going over kmer: " + str(kmer) + "/" + str(data.__len__()))
        kmer_df = data[:].loc[kmer:kmer+11]
        print("finding codons....")
        find_codon(kmer_df)
        print("translating codons....")
        translate_codon(kmer_df)
        print("Sets the Mutation type....")
        kmer_df['Mutation Type'] = kmer_df[['Ref_Protein', 'Potential_Protein']].apply(lambda protein: check_mutation_type(protein[0], protein[1]), axis=1)
    print("After a long for loop")
    file_name += ".with.mutation.type.txt"
    data.to_csv(file_name, sep='\t', encoding='utf-8')
    print("The File is ready in the folder")
    print("--- %s sec ---" % (time.time() - start_time))
    return data


def freqs_to_dataframe(freqs_file):
    data = pd.read_table(freqs_file)
    data = data[data.Ref != '-']
    data = data[data.Base != '-']
    data.reset_index(drop=True, inplace=True)
    return data


def check_12mer(data):
    if data.__len__() % 12 != 0:
        print('The data is not in 12 mer. arrange the data')
    else:
        print('The data is in 12 mer. All OK.')
        print('The length of data is:', data.__len__())


def find_codon(data):
    data['Codon'] = ""
    x = 4
    first_index_kmer = data.first_valid_index()
    for row in data.itertuples(index=True):
        if row.Index >= first_index_kmer and row.Index <= first_index_kmer+3:
            data['Codon'][row.Index] = data["Base"].loc[row.Index] + data["Ref"].loc[first_index_kmer + x] + data["Ref"].loc[(first_index_kmer + (x * 2))]
        elif row.Index >= first_index_kmer+4 and row.Index <= first_index_kmer+7:
            data['Codon'][row.Index] = data["Ref"].loc[first_index_kmer] + data["Base"].loc[row.Index] + data["Ref"].loc[first_index_kmer + (x * 2)]
        elif row.Index >= first_index_kmer+8 and row.Index <= first_index_kmer+11:
            data['Codon'][row.Index] = data["Ref"].loc[first_index_kmer] + data["Ref"].loc[first_index_kmer + x] + data["Base"].loc[row.Index]
    return data


def translate(seq):
    seq = Seq(seq)
    protein = seq.translate()
    return protein


def translate_codon(data):
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
    Mutation_Type = ""
    if protein1 == protein2:
        Mutation_Type = "Synonymous"
    elif protein1 != protein2:
        Mutation_Type = "Non-Synonymous"
    if protein2 == '*':
        Mutation_Type = "Premature Stop Codon"
    return Mutation_Type

#Functions
"""Graphs"""
#1.Distribution graph (=bowtie2 results)


def distribution_graph(val, ax, virus):
    """
    makes distribution graph
    :param val:
    :param ax:
    :return:
    """
    column_number = 3

    ind = np.arange(column_number)  # the x locations for the groups
    width = 0.35  # the width of the bars
    rects1 = ax.bar(ind, val, width, color='DarkOrchid')
    # add some text for labels, title and axes ticks
    ax.set_ylabel('% of Reads')
    # ax.set_yticklabels()
    # ax.set_title('Overall Alignment Rate')
    ax.set_xticks(ind)
    ax.set_xticklabels((virus, 'mRNA', 'rRNA'))
    labelsy = np.arange(0, 50, 10)
    ax.set_yticklabels(labelsy)
    sns.set_style("darkgrid")
    ax.set_xlim(-0.5, 3)
    return ax


#2. Multiple tandem repeat graph (repeat len)
def repeat_len_graphs(results, ax):
    """
     makes a repeat length graph
    :param results:
    :param ax:
    :return:
    """
    results.pop("bad", None)
    x_values = results.keys()
    y_repeat_len = [results[x]["repeat_length"] for x in x_values]

    medianprops = {'color': 'black', 'linewidth': 2}
    whiskerprops = {'color': 'black', 'linestyle': '-'}

    repeat_plt = ax.boxplot(y_repeat_len, patch_artist=True, medianprops=medianprops, whiskerprops=whiskerprops)
    for box in repeat_plt['boxes']:
        # change outline color
        box.set(color='DarkOrchid', linewidth=2)
        # change fill color
        box.set(facecolor='DarkOrchid')
    ax.set_xlabel("Number of Repeats")
    ax.set_ylabel("Fragment Size [bp]")
    labelsy = np.arange(0, 350, 50)
    labelsx = np.arange(1, 11, 1)
    ax.set_yticklabels(labelsy)
    ax.set_xticklabels(labelsx)
    sns.set_style("darkgrid")
    return ax

#3. Reads length
def read_len_graphs(results, ax):
    results.pop("bad", None)
    x_values = results.keys()
    y_read_len = [results[x]["read_length"] for x in x_values]
    medianprops = {'color': 'black', 'linewidth': 2}
    whiskerprops = {'color': 'black', 'linestyle': '-'}
    read_plt = ax.boxplot(y_read_len, patch_artist=True, medianprops=medianprops, whiskerprops=whiskerprops)
#     ax = plt.boxplot(y_read_len, patch_artist=True, medianprops=medianprops, whiskerprops=whiskerprops)
    for box in read_plt['boxes']:
#     for box in ax['boxes']:
        # change outline color
        box.set(color='DarkOrchid', linewidth=2)
        # change fill color
        box.set(facecolor='DarkOrchid')
    ax.set_xlabel("Number of Repeats")
    ax.set_ylabel("Read Length (bp)")


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
    graph = ax.plot(pos, reads, color="DarkOrchid")
    ax.set_xlabel("Position In The Genome [bp]")
    ax.set_ylabel("Number Of Reads")
    sns.set_style("darkgrid")
    ax.set_xlim(0, (len(pos)+10))
    ax.set_ylim(1000, 1000000)
    ax.set_yscale("log")


#6. Mutation Rates
def make_boxplot_mutation(data, ax):
    data['Base'].replace('T', 'U', inplace=True)
    data['Ref'].replace('T', 'U', inplace=True)
    min_read_count = 100000

    consensus_data = data[['Pos', 'Base']][data['Rank'] == 0]
    consensus_data = consensus_data[consensus_data['Pos'] == np.round(consensus_data['Pos'])]
    consensus_data['next'] = consensus_data['Base'] + consensus_data.shift(-1)['Base']
    consensus_data['prev'] = consensus_data.shift(1)['Base'] + consensus_data['Base']
    consensus_data['Pos'] = consensus_data[['Pos']].apply(pd.to_numeric)

    data['counts_for_position'] = np.round(data['Read_count'] * data['Freq'])
    data = data[data['Pos'] == np.round(data['Pos'])]  # remove insertions
    data['Pos'] = data[['Pos']].apply(pd.to_numeric)
    data = data[data['Read_count'] > min_read_count]
    data['mutation_type'] = data['Ref'] + data['Base']
    data['mutation_class'] = np.where(data["Rank"] == 0, "self", np.where(data['mutation_type'] == 'GA', 'transition', np.where(data['mutation_type'] == 'AG', 'transition', np.where(data['mutation_type'] == 'CU', 'transition', np.where(data['mutation_type'] == 'UC', 'transition','transversion')))))

    sns.set_palette(sns.color_palette("Paired", 12))
    data["complement"] = 1 - data["Freq"]
    data = data[data['Ref'] != data['Base']]
    data = pd.merge(data, consensus_data[['Pos', 'next', 'prev']], on="Pos")

    data = data[data["Base"] != "-"]
    data['abs_counts'] = data['Freq'] * data["Read_count"]
    data['Frequency'] = data['abs_counts'].apply(lambda x: 1 if x == 0 else x) / data["Read_count"]
    data["Mutation"] = data["Ref"] + "->" + data["Base"]
    g = sns.boxplot(x="Mutation Type", y="Frequency", hue="Mutation", data=data[data["Base"] != "-"],
                    hue_order=["C->U", "U->C", "G->A", "A->G", "C->A", "G->U", "U->G", "U->A", "G->C", "A->C", "A->U",
                               "C->G"], order=["Synonymous", "Non-Synonymous", "Premature Stop Codon"])
    g.set(yscale="log")
    plt.legend(bbox_to_anchor=(1.0, 1), loc=2, borderaxespad=0.)
    ax.set_ylim(10 ** -6, 1)


def main():
    # 1. Get freqs file and CirSeq running directory.
    out_dir = 'F:/Cirseq/RV/20170322_output_all_23_qscore'
    tmp_cirseq_dir = out_dir + '/tmp'
    freqs_file = out_dir + '/RVB14p2.freqs'
    out_plots_dir = out_dir + '/plots/'


    """2. Analyze those file and directory to get the number of tandem repeats of the cirseq,
    repeats length and the amount of reads per repeat"""

    if os.path.isfile(out_dir + '/repeat_summery.npy'):
            repeats_dict = np.load(out_dir + '/repeat_summery.npy', encoding='latin1').item() #dict for read vs. repeat
    else:
            repeats_dict = get_repeats_num(tmp_cirseq_dir, out_dir) #dict for read vs. repeat

    if os.path.isfile(out_dir + '/results.npy'):
            read_and_repeat_length = np.load(out_dir + '/results.npy', encoding='latin1').item() #dict for read and repeat length
    else:
            read_and_repeat_length = get_read_and_repeat_length(tmp_cirseq_dir, out_dir) #dict for read and repeat length

    """3. Get the coverage of the CirSeq run"""

    """ 4. Adding mutation types to the freqs file"""
    if not os.path.isfile(freqs_file + ".with.mutation.type.txt"):
         append_mutation = find_mutation_type(freqs_file)

    mutation_file = freqs_file + ".with.mutation.type.txt"
    mutation_rates = freqs_to_dataframe(mutation_file)


    """5. Run bowtie2 for human rRNA, mRNA and the virus"""

    #RV
    values = (19.32, 40.52, 45.17)


    # 6. Subplots all relevant graph in subplots
    #     6.1. Distribution graph (=bowtie2 results)
    #     6.2. Multiple tandem repeat graph (repeat len)
    #     6.3. Reads length
    #     6.4. Amount of reads per repeat (Reads VS. Repeats Graph)
    #     6.5. Coverage

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

    distribution_graph(values, ax0, 'RV')
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
    # ax5.set_title('Mutation Rates')

    plt.show()
    #     plt.savefig('C:/Users/Oded/Desktop/all.png', dpi=300)
    #     plt.close("all")
    #     print("The Plot is ready in the folder")

if __name__ == "__main__":
    main()