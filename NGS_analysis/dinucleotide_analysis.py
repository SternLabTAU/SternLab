import pandas as pd
import os
import matplotlib.pyplot as plt
import seaborn as sns
import argparse
import pylab
sns.set_context("poster")
sns.set_style('white')

''' constants for MS2  '''
primers = list(range(1, 20)) + list(range(1291, 1304)) + list(range(1179, 1200)) + list(range(2270, 2288)) +\
                list(range(2167, 2188)) + list(range(3548, 3570))
problematic = [17, 18, 19, 20, 21, 22, 23, 183, 188, 189, 190, 196, 274, 317, 364, 452, 762, 2719, 3117, 3133, 3139,
                   3143, 3146, 3150, 3401, 3539, 3542]

proteins = {"A_protein":list(range(130,1311)), "CP_protein":list(range(1335,1727)),\
			 "lysis_protein":list(range(1678,1905)),"replicase_protein":list(range(1761,3398))}

THRESHOLD = 0.0004
aa_palette = ['#e6194b', '#3cb44b', '#ffe119', '#0082c8	', '#f58231', '#911eb4', '#46f0f0', '#f032e6',
              '#d2f53c', '#fabebe', '#008080', '#e6beff', '#aa6e28', '#fffac8', '#800000', '#aaffc3	',
              '#808000', '#ffd8b1', '#000080', '#808080']

#pal = sns.set_palette(aa_palette)
#pal = sns.palplot(sns.color_palette("cubehelix", 20))


def main(args):
    df = pd.read_csv(args.in_file, sep='\t')
    location = args.location

    with open(args.reference, 'r') as o:
        reference = o.read()

    if location == 'next':
        df = add_next_nucleotide_2_mutation_file(df, reference)
        #df = filter_for_minor_allele_only(df_next)
    elif location == 'prev':
        df = add_previous_nucleotide_2_mutation_file(df, reference)
        #df = filter_for_minor_allele_only(df_prev)
    else:
        if location != 'both':
            raise Exception('Invalid location input, should be "both"\n')

        df_next = add_next_nucleotide_2_mutation_file(df, reference)
        df_prev = add_previous_nucleotide_2_mutation_file(df, reference)
        df = pd.concat([df_next, df_prev])



    df.reset_index(drop=True, inplace=True)

    if args.out != None:
        base = os.path.basename(args.in_file) + "_dinucleotide_{}".format(location)
        df.to_csv(os.path.join(args.out, base), sep='\t', index=False)

    # hard coded for MS2
    pos_2_remove = primers + problematic
    filter_positions(df, pos_2_remove=pos_2_remove)

    # plot
    time=13
    df = df[df['Time'] == time]
    plot_boxplot_dinucleotides_by_type(df, time)

    df = df[df['Freq'] >= THRESHOLD]

    #plot_count_dinucleotide_by_mutation(df, time)
    #plot_count_dinucleotide_by_type(df, time)
    #plot_count_dinucleotides_by_degree(df, time)
    plot_count_dinucleotide_by_mutatnt_aa(df, time)



def add_previous_nucleotide_2_mutation_file(df, reference):
    """
    This method adds to each position in df the previous nucleotide under Prev_nucleotide column
    :param df: a data frame containing mutations
    :param reference: a reference genome
    :return: a new data frame with Prev_nucleotide column
    """

    previous = []
    for pos in df['Pos'].values:
        if pos == 1:
            previous.append(None)   # no previous nucleotide for position 1
        else:
            previous.append(reference[pos -2])  # reference starts with zero

    df["Prev_nucleotide"] = previous
    df["Di_nucleotide"] = df['Prev_nucleotide'] + df["Ref"]
    df.drop("Prev_nucleotide", axis=1, inplace=True)
    return df


def add_next_nucleotide_2_mutation_file(df, reference):
    """
    This method adds to each position in df the next nucleotide under Next_nucleotide column
    :param df: a data frame containing mutations
    :param reference: a reference genome
    :return: a new data frame with Next_nucleotide column
    """

    next = []
    for pos in df['Pos'].values:
        if pos == len(reference):   # last position in reference has no next
            next.append(None)
        else:
            next.append(reference[pos])  # reference starts with zero

    df["Next_nucleotide"] = next
    df["Di_nucleotide"] = df['Ref'] + df["Next_nucleotide"]
    df.drop("Next_nucleotide", axis=1, inplace=True)
    return df


def filter_positions(df, pos_2_keep=None, pos_2_remove=None):
    """
    This method filter positions from a data frame
    :param df: data frame containing 'Pos' column
    :return: a filtered data frame
    """
    if pos_2_keep != None:
        df = df[df['Pos'].isin(pos_2_keep)]
    if pos_2_remove != None:
        df = df[~df['Pos'].isin(pos_2_remove)]

    return df


def filter_for_minor_allele_only(df):
    """
    This method filter all mutation which are not the minor alleles in each position in the genome
    :param df: a data frame of mutations
    :return: a data frame containing only minor allele mutations
    """
    df = df[(df['Rank'] == 0) | (df['Rank'] == 1)]  # get minor allele mutations. take rank 0 for
                                                    # fixed mutations
    return df

def plot_count_dinucleotides_by_degree(df, time):

    grouped = df.groupby(["Degree","Di_nucleotide", "Replica"]).size().reset_index(name="Count")
    sns.set(style="ticks")
    nuc_order = list(set(grouped['Di_nucleotide']))
    grid = sns.FacetGrid(grouped, col="Degree", row="Replica", row_order=['A', 'B'], hue="Degree",
                         palette="YlGnBu")  # palette="RdYlGn"
    grid.map(sns.barplot, "Di_nucleotide", "Count", order=nuc_order)
    grid.fig.suptitle("Count by degree Passage {}".format(time), fontsize=18)
    pylab.get_current_fig_manager().window.showMaximized()
    plt.show()




def plot_count_dinucleotide_by_type(df, time):

    grouped = df.groupby(["Degree", "Di_nucleotide", "Replica", "Type"]).size().reset_index(name="Count")
    sns.set(style="ticks")
    nuc_order = list(set(grouped['Di_nucleotide']))
    grid = sns.FacetGrid(grouped, col="Degree", row="Replica", row_order=['A', 'B'], hue="Type",
                         palette="Set2")
    grid.map(sns.barplot, "Di_nucleotide", "Count", order=nuc_order)
    sns.plt.legend(loc='upper left', bbox_to_anchor=(-1.35, 1.5))
    pylab.get_current_fig_manager().window.showMaximized()
    grid.fig.suptitle("Count by mutation type Passage {}".format(time), fontsize=18)
    plt.show()


def plot_count_by_freq_and_coverage_degree(df, time):
    """
    plots the count of each di nucleotide as a factor of freq * coverage
    :param df: di nucleotide data
    :param time: a passage to filter by
    :return: plot
    """

    grouped = df[df['Time'] == time]
    grouped['Total_count'] = df['Read_count'] * df['Freq']
    nuc_order = list(set(grouped['Di_nucleotide']))
    grid = sns.FacetGrid(grouped, col="Degree", row="Replica", row_order=['A', 'B'], hue="Degree",
                         palette="YlGnBu")  # palette="RdYlGn"

    grid.map(sns.barplot, "Di_nucleotide", "Total_count", order=nuc_order)
    pylab.get_current_fig_manager().window.showMaximized()
    plt.show()
    

def plot_count_dinucleotide_by_mutation(df, time):

    g = df.copy()
    g['Mutation'] = g['Mutation'].str[0] + '->' + g['Mutation'].str[1]
    grouped = g.groupby(["Degree", "Di_nucleotide", "Replica", "Mutation"]).size().reset_index(name="Count")
    sns.set(style="ticks")
    nuc_order = list(set(grouped['Di_nucleotide']))
    grid = sns.FacetGrid(grouped, col="Degree", row="Replica", row_order=['A', 'B'], hue="Mutation",
                         palette="Set3")
    grid.map(sns.barplot, "Di_nucleotide", "Count", order=nuc_order)
    sns.plt.legend(loc='upper left', bbox_to_anchor=(-1.35, 1.5))
    pylab.get_current_fig_manager().window.showMaximized()
    grid.fig.suptitle("Count by mutation identity Passage {}".format(time), fontsize=18)
    plt.show()


def plot_count_dinucleotide_by_mutatnt_aa(df, time):


    grouped = df.groupby(["Degree", "Di_nucleotide", "Replica", "mutAA"]).size().reset_index(name="Count")
    sns.set(style="ticks")
    pal = sns.palplot(sns.color_palette("Set2", 24))
    nuc_order = list(set(grouped['Di_nucleotide']))
    grid = sns.FacetGrid(grouped, col="Degree", row="Replica", row_order=['A', 'B'], hue="mutAA",
                         palette=pal)
    grid.map(sns.barplot, "Di_nucleotide", "Count", order=nuc_order, alpha=0.5)
    sns.plt.legend(loc='upper left', bbox_to_anchor=(-1.35, 1.5))
    pylab.get_current_fig_manager().window.showMaximized()
    grid.fig.suptitle("Count by animo acid Passage {}".format(time), fontsize=18)
    plt.show()




def plot_boxplot_dinucleotides_by_type(df, time):

    grouped = df.groupby(["Degree", "Di_nucleotide", "Replica", "Type"]).apply(lambda x : x)
    sns.boxplot(x="Type", y="Freq", hue="Di_nucleotide", data=grouped)
    plt.yscale("log")
    plt.ylim(0.0000001, 1)
    plt.title("Mutations frequencies as a function of genomic context Passage{}".format(time))
    sns.plt.legend(loc='upper left',bbox_to_anchor=(1.0, 0.7))
    pylab.get_current_fig_manager().window.showMaximized()
    plt.show()



if __name__ == "__main__":
    parser = argparse.ArgumentParser()
    parser.add_argument("-i", "--in_file", type=str, help="a path to a mutation file", required=True)
    parser.add_argument("-r", "--reference", type=str, help="a path to a reference file", required=True)
    parser.add_argument("-o", "--out", type=str, help="a path to an output directory in which the results will be saved", required=False)
    parser.add_argument("-l", "--location", type=str,
                        help="an indicator for the nucleotide to take; can be next, prev or both", default='next')
    args = parser.parse_args()
    main(args)