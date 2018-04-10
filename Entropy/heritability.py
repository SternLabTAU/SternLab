import pandas as pd
import numpy as np
import os
from collections import Counter
from tqdm import tqdm
import re
from Bio import SeqIO, Phylo
import math
from Bio.Seq import Seq
import random
from ete3 import Tree, TextFace, TreeStyle
from itertools import combinations
from scipy import stats
import matplotlib.pyplot as plt
import seaborn as sns
from utils import *


'''
This code implements the statistics behind Blomberg paper regarding testing for phylogenetic signals.

* randomization test
* k statistic
'''



def contrasts(tree, mapping):
    ''' calculates the variance of all cijs in a given tree'''

    all_contrasts = []

    #  map leaf to entropy value
    term_names = [term.name.split('.')[0] for term in tree.get_terminals()]
    entropy_values = [(name, mapping[mapping['refseq_id'] == name]['entropy_5'].values[0]) for name in term_names]

    # calculate all pairs
    pairs =  list(combinations(entropy_values, 2))
    for pair in pairs:
        sqrt_distance = np.sqrt(tree.distance(pair[0][0], pair[1][0]))
        if sqrt_distance == 0:
            sqrt_distance = 10**-5  # consider ignoring those?
        entropy_diff = np.abs(pair[0][1] - pair[1][1])

        c = entropy_diff / sqrt_distance
        all_contrasts.append(c)

    # have all contrasts, now can calculate the variance
    variance = np.var(all_contrasts)
    return variance


def randomization_test(tree, mapping, n):
    ''' shuffle the  tree n times and calculate variance'''

    vars = []
    # shuffle tree names along the tree n times, and calculate the contrasts variance
    for i in tqdm(range(n)):
        new_tree = tree
        terms = [term.name.split('.')[0] for term in tree.get_terminals()]   # poping from the list, need a new list each iterations
        for node in new_tree.get_terminals():
            node.name = terms.pop(random.randrange(len(terms)))

        # now we have a new tree and we can calculate variance of the contrasts
        var = contrasts(tree, mapping)
        vars.append(var)

    return vars

def check_signal_significant(true_tree_variance, shuffled_variance, alias, out=None):
    ''' test if the variance is significantly low. shuffled variance is a list of shuffled trees contrasts variances'''

    # need to test whether 95% of the simulations are greater then the true tree variance. look for the 5 percentile
    # and check if the true value is lower or equal

    significance = False
    shuffled_variance.append(true_tree_variance)

    cutoff = np.percentile(shuffled_variance, 5)
    if true_tree_variance <= cutoff:
        significance = True

    if out:
        sns.kdeplot(np.asarray(shuffled_variance), shade=True, color='#378F9E')
        plt.axvline(cutoff, color='r', linestyle='--', alpha=.5)
        plt.axvline(true_tree_variance, color='g', linestyle='-', alpha=.75)
        plt.xlabel('Contrasts variance')
        plt.title('{}\nRandomization test distribution result, tree significant = {}'.format(alias,significance))
        plt.savefig(os.path.join(out, 'randomization_test.png'), format='png')

    else:
        sns.kdeplot(np.asarray(shuffled_variance), shade=True, color='#378F9E')
        plt.axvline(cutoff, color='r', linestyle='--', alpha=.5)
        plt.axvline(true_tree_variance, color='g', linestyle='-', alpha=.75)
        plt.xlabel('Contrasts variance')
        plt.title('{}\nRandomization test distribution result, tree significant = {}'.format(alias, significance))
        plt.show()

    return significance


def run_randomization_test(super_folder, n, mapping):

    all_trees = []
    results = []

    for root, dirs, files in tqdm(os.walk(super_folder)):
        tree = [f for f in files if 'phyml_tree' in f]
        if tree != []:
            all_trees.append(os.path.join(root, tree[0]))

    for tree in tqdm(all_trees):
        if tree_2_string(tree) == '':
            continue
        alias = os.path.basename(tree).split('.')[0].strip()
        print(alias)
        t = Phylo.read(tree, 'newick')
        num_leafs = len([term.name.split('.')[0] for term in t.get_terminals()])
        t_var = contrasts(t, mapping)
        random_vars = randomization_test(t, mapping, n)
        out = os.path.dirname(tree)
        sgnf = check_signal_significant(t_var, random_vars, alias, out)
        df = pd.DataFrame({'group': alias, 'significance':sgnf, 'num_leafs':num_leafs})
        results.append(df)

    concat = pd.concat(results)
    concat['estimation'] = concat['num_leafs'].apply(lambda x : 'underestimation' if x > 200 else 'estimated', axis=1)
    concat.to_csv(os.path.join(super_folder, 'significance_by_family.csv'), index=False)

    return


def main():

    main_dir = r'/Volumes/STERNADILABHOME$/volume1/daniellem1/Entropy/data/Phylogeny/family'
    virus_entropy = r'/Users/daniellemiller/Google Drive/Msc Bioinformatics/Projects/entropy/virus_host_entropy/virus_host_entropy.csv'

    mapping = pd.read_csv(virus_entropy)

    n = 1000

    run_randomization_test(main_dir, n=n , mapping=mapping)
    print('Done!')


if __name__ == '__main__':
    main()

