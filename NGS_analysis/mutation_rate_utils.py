import numpy as np
import matplotlib.pyplot as plt
import pandas as pd
import seaborn as sns
import statsmodels.api as sm
import os
import pickle
from datetime import datetime

TIME_0_RATE = 0.00001	# 10^-5

def get_interquartile_range(df):
    """
    This method filter df to include only mutations who's frequencies are in the interquartile range (25%-75%)
    :param df: a mutation data frame
    :return: data frame filtered as described
    """

    passages = set(df.Time.values)
    all_filtered = []
    cnt = 0

    for passage in passages:
        lower_boundary = df['Freq'][df.Time == passage].quantile(0.25)
        upper_boundary = df['Freq'][df.Time == passage].quantile(0.75)

        filtered = df[(df.Time == passage) & (df.Freq >= lower_boundary) & (df.Freq <= upper_boundary)]
        all_filtered.append(filtered)
        cnt += filtered.shape[0]

    result = pd.concat(all_filtered)
    if cnt != result.shape[0]:
        raise Exception('Error: Dimensions of concatenated data frames do not fit\n')

    return result



def filter_freqs_2_regression(df, pos_2_remove=None, pos_2_keep=None, threshold=0.0):
    """
    This method filters a mutation data frame to be used as the input of the linear regression
    :param df: mutation data frame. should be in the mutation files format
    :param pos_2_remove: position to remove from calculation
    :param pos_2_keep: positions to keep - most times will be the coding regions
    :param threshold: frequency threshold
    :return: a data frame of all mutations to consider in the linear regression
    """

    # filter data
    df = df[(df['Pos'].isin(pos_2_keep)) & (~df['Pos'].isin(pos_2_remove)) & (df['Type'] == 'synonymous') & (df['Freq'] >= threshold)]

    # take mutations in 25% - 75% interval
    df = get_interquartile_range(df)

    # add a unique ID to filter positions which appear less then X times
    df['ID'] = df['Pos'].astype(str) + '_' + df['Mutation']
    filtered = df.groupby('ID').filter(lambda x: len(x) >= 3)

    return filtered


def fit_regression(x,y, points=2):
    """
    Fit a linear regression model
    :param x: passages
    :param y: correspond frequencies
    :param points: minimum number of points for regression calculation 
    :return: a tuple (slope, P value)
    """

    assert((len(x) == len(y)) and len(x) > points)
    slope, intercept, r_value, p_value, std_err = stats.linregress(x, y)

    return (slope, p_value)

def add_zero_time_point(x,y, rate):
    """
    This method adds zero time point to a given vector. x will be added '0' element and y will be added 'rate'
    :param x: a vector of passages
    :param y: a vector of frequencies
    :return: a tuple of (x,y) updated
    """
    x = np.append(x, 0)
    y = np.append(y, rate)

    return (x,y)


def calculate_regression_slopes(df, slopes, p_values, mutation_type=None, limit=3):
    """
    This method calculates the regression slopes for all mutations or a given mutation type
    :param df: containing more then one time point for each mutation.
    :param mutation_type: a list of The mutation identity e,g 'AG'. for calculating mutation rate of specific mutation type. 
                        for example - calculating only transitions mutation rate
    :param slopes: a dictionary to contain all slopes
    :param p_values: a dictionary to contain all P values
    :return: void - updated dictionaries
    """

    if mutation_type != None:
        df = df[df['Mutation'].isin(mutation_type)]
    # for each position get x any y vectors and fit regression
    for pos in np.unique(df.Pos.values):
        pos_df = df[df['Pos'] == pos]

        if pos_df.shape[0] < limit:
            raise Exception('Not enough time points for position {}'.format(pos))

        x = pos_df['Time'].values
        y = pos_df['Freq'].values

        current_mutation_type = pos_df['Mutation'].values[0]

        #with_zero_pnt = add_zero_time_point(x,y, TIME_0_RATE)
        #x = with_zero_pnt[0]
        #y = with_zero_pnt[1]

        assert (len(x) == len(y))

        reg_result = fit_regression(x,y)
        slope = reg_result[0]
        p_value = reg_result[1]

        slopes[current_mutation_type].append(slope)
        p_vals[current_mutation_type].append(p_value)


def get_median_slope(slopes, mutation_type=None):
    """
    calculate median slope values
    :param slopes: a dictionary
    :param mutation_type: mutation type to calculate according to
    :return: median slope value
    """
    if mutation_type != None:
        return np.median(slopes[mutation_type])
    return np.median(np.asarrsy(slopes.values()))


def plot_regression(x,y,slope,intercept):
    plt.scatter(x,y,color='dark orange')
    plt.plot(x, slope*x +intercept, '-')




def get_conserved_positions(position_pickle):
    """
    This method opens a pickle file containing a list of positions which are 100% conserved (the same nt in the MSA)
    :param position_pickle: a full path to a pickle file
    :return: a list of conserved positions
    """

    if position_pickle == None:
        raise Exception("Invalid input for get conserved positions, received null file path\n")
    if not os.path.exists(position_pickle):
        raise Exception("Conserved position file do not exists\n")

    with open (position_pickle, 'rb') as o:
        conserved_positions = pickle.load(o)

    return conserved_positions



#TODO - add the plots labels
def plot_all_slopes(slopes, out, max_passage):
    """
    This method plots all mutations median slopes
    :param slopes: a dictionary containing for each mutation type its slopes
    :param out: an output directory in which files will be saved
    :param max_passage: the last time point of the experiment
    :return: plot
    """
    order = ['AA', 'AC', 'AG', 'AT', 'CA', 'CC', 'CG', 'CT', 'GA', 'GC', 'GG', 'GT', 'TA', 'TC', 'TG', 'TT']
    mutations = ['AG', 'AC', 'AT', 'GA', 'GC', 'GT', 'CA', 'CG', 'CT', 'TG', 'TA', 'TC']
    median_slopes = [np.median(slopes[mutation]) for mutation in mutations]

    df1 = pd.DataFrame({"Mutation": mutations, "Slope": median_slopes})
    df2 = df1.copy()

    # in each data frame create two columns (x,y) to represent a point in the mutation regression line
    df1["x"] = 0
    df2["x"] = max_passage

    df1["y"] = TIME_0_RATE
    df2["y"] = df2.x * df2.Slope + TIME_0_RATE

    # create a third df of non mutations AA CC GG TT
    df3 = pd.DataFrame({"Mutation": ['AA', 'CC', 'GG', 'TT'], "Slope": [0, 0, 0, 0], "x": [0, 0, 0, 0], "y": [0, 0, 0, 0]})

    # unite all data frames and sort according to mutation type
    df = pd.concat([df1, df2, df3])
    df = df.sort_values(by="Mutation")
    df.reset_index(drop=True, inplace=True)

    # create a slopes vector corresponds to mutations in order list
    #all_slopes = [0] * len(order)
    #for i, mut in enumerate(order):
     #   all_slopes[i] = slop

    # create the plot
    sns.set(style="ticks")
    grid = sns.FacetGrid(df, col="mut", hue="all_slopes", col_order=order, col_wrap=4, size=2)
    grid.set(xticks=np.arange(21), yticks=[0.000000001, 0.01], xlim=(0, 21), ylim=(0.000000001, 0.01))
    grid.map(plt.plot, "x", "y", marker="o", ms=4)
    # grid.fig.suptitle("Mutation rate divided by mutation identity " + str(degree) + str(replica))
    plt.yscale("log")

    # set x and y labels
    grid.axes[0, 0].set_xlabel('A')
    grid.axes[0, 1].set_xlabel('C')
    grid.axes[0, 2].set_xlabel('G')
    grid.axes[0, 3].set_xlabel('T')
    grid.axes[0, 0].set_xlabel('A')
    grid.axes[1, 0].set_xlabel('C')
    grid.axes[2, 0].set_xlabel('G')
    grid.axes[3, 0].set_xlabel('T')


    # set titles
    for ax, title in zip(grid.axes.flat, slopes):
        ax.set_title("slope = " + str(title))

    plt.savefig(os.path.join(out, 'Mutation_rates.png'), format='png', dpi=600)
    sns.plt.show()

