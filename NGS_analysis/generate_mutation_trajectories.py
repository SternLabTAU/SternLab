import pandas as pd
import argparse
import matplotlib.pyplot as plt
import seaborn as sns
import pylab


# MS2  positions to discard
problematic = [17,18,19,20,21,22,23,183,188,189,190,196,274,317,364,452,762,2719,3117,3133,3139,3143,3146,3150,3401,3539,3542]
coding_regions = list(range(130,1312)) + list(range(1335,3399))


def main(args):
    df = pd.read_csv(args.in_file, sep='\t')
    #df = df[df["Pos"].isin(list(range(130,1312)))]
    m_type = args.mutation_type

    if args.line != None:
        line= args.line
    else:
        line = None
    if args.transitions_only != None:
         if args.transitions_only == 'Y':
             transitions_only = True
    else:
        transitions_only=False

    # filter data frame and plot by position
    df = filter_df(df, m_type, transitions_only, line)
    plot_trajectory(df, m_type, line)


def get_trajectory(df, pos):
    """
    This method generates a trajectory of a specific mutation
    :param df: a data frame with mutations time series data. should include 'Mutation' and 'Type' field. assuming NO
                duplicates (one line).
    :param pos: a position in the genome
    :param mutation_type: the type of the mutation: synonymous, non-synonymous or stop
    :param transitions: take only transitions
    :return: a data frame of the frequencies of the specific mutation over time
    """

    df = df[(df["Pos"] == pos)].sort_values(by=['Time'])
    return df




def filter_df(df, mutation_type, transitions=False, line=None):
    """
    This method filters a data frame of mutations.
    :param df: a data frame of mutations frequencies and positions
    :param mutation_type: synonymous, non-synonymous, stop
    :param transitions: include only transitions, default is False
    :param line: the line to filter by suold be in the following format: DEGREE REPLICA connected with no spaces.default None
    :return: a data frame after filtration
    """

    if df is None or mutation_type == None:
        raise Exception("Invalid input, received a object")
    if len(line) != 3:
        raise Exception("Invalid line format, should contain degree and replica only")

    # parse line if received
    if line != None:
        degree = int(line[:2])
        replica = line[-1]
        df = df[(df['Type'] == mutation_type) & (df['Replica'] == replica) & (df['Degree'] == degree)]
    else:
        df = df[df['Mutation'] == mutation_type]

    if transitions:
        df = df[(df['Mutation'] == 'AG') | (df['Mutation'] == 'GA') | (df['Mutation'] == 'CT') | (df['Mutation'] == 'TC')]

    df = df[(df["Pos"].isin(coding_regions)) & (~df["Pos"].isin(problematic))]

    return df

def plot_trajectory(df, mutation_type, line=None):
    """
    plots the frequency of a trajectory over time
    :param df: data frame containing time points and their corresponding frequencies
    :param mutation_type: synonymous, non-synonymous or stop. will be used for the plot title
    :param line: default None, will be used for the plots title
    :return: plot
    """
    ax = sns.pointplot(x="Time", y="Freq", hue="Pos", data=df, legend=False)
    plt.title("{} trajectories".format(mutation_type))
    ax.legend_.remove()

    if line != None:
        plt.title("{} trajectories at {}".format(mutation_type, line))

    plt.ylabel("Frequency")
    plt.xlabel("Passage")
    pylab.get_current_fig_manager().window.showMaximized()
    plt.show()


if __name__ == "__main__":
    parser = argparse.ArgumentParser()
    parser.add_argument("-i", "--in_file", type=str, help="a path to an input mutation file", required=True)
    parser.add_argument("-m", "--mutation_type", type=str, help="type of mutations", required=False)
    parser.add_argument("-l", "--line", type=str, help="biological line, constructed from degree and replica concatenated", required=False)
    parser.add_argument("-t", "--transitions_only", type=str, help="transition boolean, Y to include, N o.w", required=False)
    args = parser.parse_args()
    main(args)
