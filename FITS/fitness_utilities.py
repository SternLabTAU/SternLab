import pandas as pd
import matplotlib.pyplot as plt
import pylab
import seaborn as sns
import numpy as np
from scipy import stats


''' Auxiliary functions for analysing fitness as received from fits'''

def remove_sign(df):
    """
    This method removes the '?' sign from a fitnes category label
    :param df: a data frame containing a "Category" field with
    :return: data frame without the '?' sign
    """
    catg = df.Category.values
    for i, word in enumerate(catg):
        if "?" in word:
            catg[i] = word[1:]
    df["Category"] = catg
    return df


def plot_correlation(df):
    """
    This method plots the correlation between two degrees of MS2
    :param df: a data frame containing results from FITS runs
    :return: plot the correlation between both degrees
    """
    d37 = df[df.Degree == 37].Fitness_median.values
    d41 = df[df.Degree == 41].Fitness_median.values

    # test randomization
    # np.random.shuffle(d37)
    # np.random.shuffle(d41)

    hexplot = sns.jointplot(d37, d41, kind="hex", color="#4CB399")
    sns.plt.subplots_adjust(left=0.2, right=0.8, top=0.8, bottom=0.2)
    cax = hexplot.fig.add_axes([.85, .25, .01, .4])
    sns.plt.colorbar(cax=cax)
    plt.show()


def plot_degree_fitness_diff(df):
    """
    This method plots the difference between degrees in the MS2 experiment
    :param df: a data frame containing "Degree", "Pos" and "Fitness_median" field
    :return: a plot of fitness values separated by degree
    """
    sns.set_style("darkgrid")
    plt.hist(x="Pos", y = "Fitness_median", hue="Degree", capsize=0.2, data=df, palette="Reds")
    plt.title("Fitness valuse distribution separated by degree")
    plt.ylabel("Fitness values")
    plt.xlabel("Position in the genome (bp)")
    plt.show()


def plot_degree_category_diff(df):
    """
    This method plots the differences in categories by degree
    :param df: a data frame containing "Degree", "Pos" and "Category" field
    :return: a plot of fitness category separated by degree
    """
    sns.set_style("darkgrid")
    sns.barplot(x="Pos", y ="Category", hue="Degree", data=df, palette="Paired")
    plt.title("Fitness categories distribution separated by degree")
    plt.ylabel("Fitness category")
    plt.xlabel("Position in the genome (bp)")
    plt.show()

def plot_degree_replica_diff(df, degree=None):
    """
    This method plots the differences in fitness values between two replicas
    :param df: a data frame with "Replica", "Pos" and "Fitness_median" fields
    :param degree: optional. default = None. used for MS2
    :return:
    """
    if degree != None:
        df = df[df.Degree == degree]

    sns.set_style("darkgrid")
    sns.barplot(x="Pos", y="Fitness_median", hue="Replica", data=df)
    plt.title("Fitness value distribution separated by replica")
    plt.ylabel("Fitness value")
    plt.xlabel("Position in the genome (bp)")

    if degree != None:
        plt.title("Fitness value distribution separated by replica {} degrees".format(degree))
    plt.show()

def plot_mutations_diff(df, degree=None, replica=None):
    """
    This method plots a scatter of mutations fitness values
    :param df: a data frame with "Pos", "Fitnes_median" and "Mutation" fields
    :param degree: optional, default of None. used in MS2
    :param replica: optional. biological repeat (a string)
    :return: a scatter of fitness value separated by mutation type
    """
    if degree != None:
        df = df[df.Degree == degree]
    if replica != None:
        df = df[df.Replica == replica]

    # Make a custom sequential palette using the cubehelix system
    pal = sns.cubehelix_palette(10, 2, 2, light=.6, dark=.2)

    sns.set_style("darkgrid")
    sns.lmplot(x="Pos", y="Fitness_median", hue="Mutation", data=df, fit_reg=False, palette="bright")
    plt.title("Fitness values distribution by mutation type")
    if degree != None and replica != None:
        plt.title("Fitness values distribution by mutation type {} {}".format(degree, replica))
    elif degree != None and replica == None:
        plt.title("Fitness values distribution by mutation type {} degrees".format(degree))
    else:
        plt.title("Fitness values distribution by mutation type replica {}".format(replica))
    plt.ylabel("Fitness value")
    plt.xlabel("Position in the genome (bp)")
    plt.show()

def plot_mutation_types_diff(df, degree=None, replica=None):
    """
    This method plots a scatter of mutations fitness values by mutation type - synonymous, non-synonymous
    or stop
    :param df: a data frame with mutation type information - "Type" field
    :param degree: optional. default = None
    :param replica: optional. default = None
    :return: a scatter of fitness values according to type
    """
    if degree != None:
        df = df[df['Degree'] == degree]
    if replica != None:
        df = df[df['Replica'] == replica]

    # remove un coding regions.
    df = df[df['Type'] != 'non-coding']

    # Make a custom sequential palette using the cubehelix system
    pal = sns.cubehelix_palette(10, 2, 2, light=.6, dark=.2)
    current_palette_7 = sns.color_palette("hls", 5)
    x = sns.set_palette(current_palette_7)

    sns.set_style("darkgrid")
    sns.lmplot(x="Pos", y="Fitness_median", hue="Type", data=df, fit_reg=False, palette=x, legend=False)
    plt.title("Fitness values distribution by type")
    if degree != None and replica != None:
        plt.title("Fitness values distribution by type {} {}".format(degree, replica))
    elif degree != None and replica == None:
        plt.title("Fitness values distribution by type {} degrees".format(degree))
    else:
        plt.title("Fitness values distribution by type replica {}".format(replica))
    plt.ylabel("Fitness value")
    plt.xlabel("Position in the genome (bp)")
    plt.ylim(-0.06, 1.8)
    sns.plt.legend(loc='upper left', bbox_to_anchor=(-0.25, 0.94))
    pylab.get_current_fig_manager().window.showMaximized()
    plt.show()

def plot_degree_barplot(df):
    """
    This method plots the count of mutations in each fitness category by degree.
    :param df: data frame with "Degree" field
    :return: a barplot of the count of mutations in each fitness category separated by degree
    """
    grouped = df.groupby(["Degree", "Category"])["Category"].count()
    labels = grouped.index.get_level_values(0).tolist()
    deg = grouped.index.get_level_values(1).tolist()
    amount = grouped.values.tolist()

    new_df = pd.DataFrame({"Category" :labels, "Degree" :deg, "Amount" :amount})

    sns.barplot(x="Degree", y="Amount", hue="Category", data=new_df)
    plt.title("Classification count of fitness values separated according degree")
    plt.ylabel("Count")
    plt.show()

def plot_degree_boxplot(df):
    """
    This method plots a boxplot of fitness value separated by degrees and replicas. This analysis is
     suitable for MS2 - divide by degree and replica
    :param df: a data frame with "Degree" and "Replica" fields. replicas should be 'A' and 'B'.
                not generalized.
    :return: boxplot as described above
    """
    sns.boxplot(x="Degree", y="Fitness_median", hue="Replica", hue_order=['A', 'B'], data=df ,palette="Spectral")
    sns.despine(offset=10, trim=True)
    plt.title("Fitness values separated by degree and replica")
    plt.ylabel("Fitness values")
    plt.show()


def plot_heatmap(df, replica=None, degree=False):
    """
    This method plots a heatmap of fitness values throughout the genome. unfix BUG of x labels
     which are not showing well. until fix - use the manual option of defining x limits in the IPython
     window.
    :param df: data frame with fitness results
    :param replica: optional - biological repeat. default None
    :param degree: optional. default is False. if True pivot the heatmap according degree.
                    used in MS2
    :return: a heatmap of fitness values
    """

    if replica != None:
        df = df[df.Replica == replica]

    if degree:
        heat = df.pivot("Degree", "Pos", "Fitness_median")
    else:
        heat = df.pivot("Pos", "Fitness_median")
    sns.heatmap(heat, cmap="RdBu")
    if degree:
        plt.title("Fitness values by degree")
    else:
        plt.title("Fitness values throughout the genome")
    pylab.get_current_fig_manager().window.showMaximized()
    plt.show()


def plot_dfe(df, alias=None):
    """
    This method plots an DFE - distribution of fitness effects. This method is suitable for MS2
    :param df: data frame containing fitness values
    :return: a plot of a DFE separated by replicas and degree
    """

    sns.set(style="ticks")
    bins = np.arange(0, 2, 0.054)
    grid = sns.FacetGrid(df, col="Degree", row="Replica", row_order=['A', 'B'], hue="Degree", palette="YlGnBu")  # palette="RdYlGn"
    grid.map(plt.hist, "Fitness_median", bins=bins)
    if alias != None:
        grid.fig.suptitle("Distribution of Fitness Effects {}".format(alias), fontsize=18)
    else:
        grid.fig.suptitle("Distribution of Fitness Effects", fontsize=18)
    pylab.get_current_fig_manager().window.showMaximized()
    plt.show()


def test_kolmogorov_smirnov(dist1, dist2):
    """
    apply KS test on two distributions
    :param dist1: a vector of values
    :param dist2: a vector of value
    :return: p value of the KS test
    """

    res = stats.ks_2samp(dist1, dist2)
    return res[1]

def filter_replicas_ambiguity(df):
    """
    This method removes samples which are un-matched in terms of fitness category results
    :param df: a data frame of fitness outputs
    :return: a data frame filtered from position which are ambiguous
    """

    grouped = df.groupby(["Pos", "Degree"])
    grouped = grouped["Category"].apply(lambda x: len(set(x)) != 1).reset_index(name="isAmbiguous")
    pos_2_remove = grouped[grouped["isAmbiguous"] == True]["Pos"].values
    df = df[~df["Pos"].isin(pos_2_remove)]
    return df

def discretize_fitness_dfe(df, alias):
    """
    This method displays the fitness dfe according to defined bins
    :param df: a data frame of fitness values
    :param alias: an alias which will be added to the graph title
    :return: plots a histogram of fitness values - binned
    """

    df37 = df[df.Degree == 37]
    df41 = df[df.Degree == 41]
    # define bins length
    bins = [0.0, 0.2, 0.5, 0.8, 1.1, 2.0]
    g1 = pd.cut(df41["Fitness_median"], bins)
    count1 = g1.value_counts().reindex(g1.cat.categories)  # sort the bin in increasing prder
    count1.plot(kind='bar', color='brown', alpha=0.5, label="41")

    g2 = pd.cut(df37["Fitness_median"], bins)
    count2 = g2.value_counts().reindex(g2.cat.categories)
    count2.plot(kind='bar', color='blue', alpha=0.3, label="37")
    plt.title("DFE discretization {}".format(alias))
    plt.ylabel("Count")
    plt.xlabel("Fitness value")
    plt.legend()
    plt.show()




