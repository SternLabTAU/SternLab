import pandas as pd
import os
import matplotlib.pyplot as plt
import seaborn as sns
import argparse
import pylab
sns.set_context("poster")

''' constants for MS2  '''


primers = list(range(1, 20)) + list(range(1291, 1304)) + list(range(1179, 1200)) + list(range(2270, 2288)) +\
                list(range(2167, 2188)) + list(range(3548, 3570))
problematic = [17, 18, 19, 20, 21, 22, 23, 183, 188, 189, 190, 196, 274, 317, 364, 452, 762, 2719, 3117, 3133, 3139,
                   3143, 3146, 3150, 3401, 3539, 3542]

proteins = {"A_protein":list(range(130,1311)), "CP_protein":list(range(1335,1727)),\
			 "lysis_protein":list(range(1678,1905)),"replicase_protein":list(range(1761,3398))}


def main(args):
    df = pd.read_csv(args.in_file, sep='\t')

    with open(args.reference, 'r') as o:
        reference = o.read()

    df_next = add_next_nucleotide_2_mutation_file(df, reference)
    df_next = filter_for_minor_allele_only(df_next)
    df_prev = add_previous_nucleotide_2_mutation_file(df, reference)
    df_prev = filter_for_minor_allele_only(df_prev)
    df = pd.concat([df_prev, df_next])

    if args.out != None:
        base = os.path.basename(args.in_file) + "_dinucleotide"
        df.to_csv(os.path.join(args.out, base), sep='\t', index=False)



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

def plot_count_dinucleotides_by_degree(df):

    grouped = df.groupby(["Degree","Di_nucleotide", "Replica"]).size().reset_index(name="Count")
    sns.set(style="ticks")
    grid = sns.FacetGrid(grouped, col="Degree", row="Replica", row_order=['A', 'B'], hue="Degree",
                         palette="YlGnBu")  # palette="RdYlGn"
    grid.map(sns.barplot, "Di_nucleotide", "Count")
    plt.show()

def plot_count_dinucleotide_by_type(df):

    grouped = df.groupby(["Degree", "Di_nucleotide", "Replica", "Type"]).size().reset_index(name="Count")
    sns.set(style="ticks")
    grid = sns.FacetGrid(grouped, col="Degree", row="Replica", row_order=['A', 'B'], hue="Degree",
                         palette="RdYlGn").map(sns.barplot, "Di_nucleotide", "Count").add_legend()
    sns.plt.legend(loc='upper left', bbox_to_anchor=(-0.25, 0.94))
    pylab.get_current_fig_manager().window.showMaximized()
    plt.show()

def plot_count_dinucleotide_by_mutation(df):

    grouped = df.groupby(["Degree", "Di_nucleotide", "Replica", "Mutation"]).size().reset_index(name="Count")
    sns.set(style="ticks")
    grid = sns.FacetGrid(grouped, col="Degree", row="Replica", row_order=['A', 'B'], hue="Degree",
                         palette="Set3").map(sns.barplot, "Di_nucleotide", "Count").add_legend()
    sns.plt.legend(loc='upper left', bbox_to_anchor=(-0.25, 0.94))
    pylab.get_current_fig_manager().window.showMaximized()
    plt.show()


def plot_boxplot_dinucleotides_by_type(df):

    grouped = df.groupby(["Degree", "Di_nucleotide", "Replica", "Type"]).apply(lambda x : x).reset_index()
    sns.boxplot(x="Type", y="Freq", hue="Di_nucleotides", data=grouped)
    plt.yscale("log")
    plt.show()






if __name__ == "__main__":
    parser = argparse.ArgumentParser()
    parser.add_argument("-i", "--in_file", type=str, help="a path to a mutation file", required=True)
    parser.add_argument("-r", "--reference", type=str, help="a path to a reference file", required=True)
    parser.add_argument("-o", "--out", type=str, help="a path to an output directory in which the results will be saved", required=False)
    args = parser.parse_args()
    main(args)