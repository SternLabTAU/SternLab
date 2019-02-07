import numpy as np
import pandas as pd
import random
from tqdm import tqdm
import statsmodels.stats.multitest as multi
import argparse
from Bio import SeqIO
import os

def get_critical_points(profile):
    """
    this method calculates the gradient of each one of the points in profile and returns a list of tuples containing all
    critical points (x,y)
    :param profile: a vector of numbers
    :return: a vector of the gradients of each data point
    """

    gradients = np.gradient(profile)
    critical_pts = [(x+1, profile[x]) for x, y in enumerate(gradients) if y ==0]    # add 1 to x to fit the relevant position (indices starts from 0!
    return critical_pts



def codon_scrambler(seq):
    """
    scramble the codons in seq. this method assumes that the given seq is of coding region and contains triplets.
    :param seq: a cds sequence
    :return: the scrambled codon sequence
    """

    if seq %3 != 0:
        print('seq should be a multiple of 3!')
        return

    codons = [seq[i:i+3] for i in range(0, len(seq), 3)]
    random.shuffle(codons)
    return ''.join(codons)


def stretchFinder(profile, l, m):
    """
    implementation of strechFinder as described in : "Synonymous site conservation in the HIV-1 genome"
    :param profile: a vector of entropy values
    :param l: the window size
    :param m: number of permutations
    :return:
    """
    start_index = []
    p_values = []

    for i in tqdm(range(0,len(profile) - l, l)):
        # get the current window and its average value
        w = profile[i:i+l]
        avg = np.mean(w)

        # permutation test
        avgs = []
        for j in range(m):
            new_profile = profile
            avgs.append(np.mean(new_profile[np.random.choice(len(new_profile), size=l, replace=False)]))

        # sort average in order to get the p value
        avgs.sort()
        idx = np.searchsorted(avgs, avg)
        p_value = idx/m
        p_values.append(p_value)
        start_index.append(i)

    data =  pd.DataFrame({'start':start_index, 'p_value':p_values, 'l':l})

    # correct for multiple tests
    data['corrected_pvalue'] = multi.fdrcorrection(data['p_value'])[1]

    return data


def main(args):

    df = pd.read_csv(args.data)

    c = 'seq_{}'.format(args.index - 1)
    seq = np.array(df[c].dropna())
    res = stretchFinder(seq, args.window, args.m)
    out= r'/sternadi/home/volume1/daniellem1/Entropy/stats/{}_stats.csv'.format(c)
    res.to_csv(out, index=False)





if __name__ == "__main__":
    parser = argparse.ArgumentParser()
    parser.add_argument("-i", "--index", type=int, help="array index")
    parser.add_argument("-w", "--window", type=int, help="window size", default=100)
    parser.add_argument("-m", "--m", type=int, help="total number of iterations", default=10**5)
    parser.add_argument("-d", "--data", type=str, help="file path to an input profile", required=True)

    args = parser.parse_args()

    main(args)