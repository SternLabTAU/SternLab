from FITS.create_fits_input import remove_del_ins
import matplotlib.pyplot as plt
import seaborn as sns
import pandas as pd
import argparse


primers = list(range(1, 20)) + list(range(1291, 1304)) + list(range(1179, 1200)) + list(range(2270, 2288)) +\
                list(range(2167, 2188)) + list(range(3548, 3570))
problematic = [17, 18, 19, 20, 21, 22, 23, 183, 188, 189, 190, 196, 274, 317, 364, 452, 762, 2719, 3117, 3133, 3139,
                   3143, 3146, 3150, 3401, 3539, 3542]



def test_selective_swipe(df, passage, degree, replica):
    df = df[(df['Degree'] == degree) & (df['Time'] == passage) & (df['Replica'] == replica)]
    df = remove_del_ins(df)
    df["Mutation"] = df['Ref'] + df['Base']
    df = df[(df['Mutation'] == 'AG') | (df['Mutation'] == 'GA') | (df['Mutation'] == 'CT') | (df['Mutation'] == 'TC')]
    df = df[df["Freq"] >= 0.0004]

    assert(len(set(df['Pos'].values)) == df.shape[0])

    plt.scatter(x=df['Pos'].values, y=df['Freq'].values, color='green', alpha=0.5)
    plt.title("Minor transition allele frequency as a function of position in the genome\n p{}{}{}".format(passage, degree, replica))
    plt.yscale("log")
    plt.ylim(0.0001, 1)
    plt.show()



def plot_boxplot(df, x, y, alias, xlabel, ylabel,hue=None, yscale=None, save=False, out=None):
    """
    plots a boxplot out of df according to x and y axis
    :param df: a source data frame
    :param x: a name of a df column to be on the x axis
    :param y: a name of a df column to be on the y axis
    :param alias: graph title
    :param xlabel: string. label for x axis
    :param ylabel: string. label for y axis
    :param hue: a separator for each x element. defualt None
    :param yscale: the scale of the y axis. defualt None
    :param svae: boolean. if False don't save, else, save to out.
    return: a boxplot
    """


    sns.boxplot(x=x,y=y,data=df,hue=hue)
    plt.title(alias, fontsize=18)
    plt.xlabel(xlabel, fontsize=16)
    plt.ylabel(ylabel, fontsize=16)

    if yscale != None:
        plt.yscale(yscale)

    if save:
        if out == None:
            raise Exception("Error: need to include out parameter whether save = true")
        plt.savefig(out,bbox_inches='tight',dpi=400)
    plt.show()


def plot_barplot(df, x, y, alias, xlabel, ylabel,hue=None, yscale=None, save=False, out=None):
    """
    plots a barplot out of df according to x and y axis
    :param df: a source data frame
    :param x: a name of a df column to be on the x axis
    :param y: a name of a df column to be on the y axis
    :param alias: graph title
    :param xlabel: string. label for x axis
    :param ylabel: string. label for y axis
    :param hue: a separator for each x element. defualt None
    :param yscale: the scale of the y axis. defualt None
    :param svae: boolean. if False don't save, else, save to out.
    return: a barplot
    """


    sns.barplot(x=x,y=y,data=df,hue=hue)
    plt.title(alias, fontsize=18)
    plt.xlabel(xlabel, fontsize=16)
    plt.ylabel(ylabel, fontsize=16)

    if yscale != None:
        plt.yscale(yscale)

    if save:
        if out == None:
            raise Exception("Error: need to include out parameter whether save = true")
        plt.savefig(out,bbox_inches='tight',dpi=400)
    plt.show()






