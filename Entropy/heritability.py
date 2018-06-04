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
import matplotlib
matplotlib.use('Agg')
import matplotlib.pyplot as plt
import seaborn as sns
from utils import *
import argparse


'''
This code implements the statistics behind Blomberg paper regarding testing for phylogenetic signals.

* randomization test
* k statistic
'''


def contrasts(tree, feature, mapping):
    ''' calculates the variance of all cijs in a given tree'''

    all_contrasts = []

    #  map leaf to entropy value
    term_names = [term.name.split('.')[0] for term in tree.get_terminals()]
    values = [(name, mapping[mapping['refseq_id'] == name][feature].values[0]) for name in term_names]

    # calculate all pairs
    pairs =  list(combinations(values, 2))
    for pair in pairs:
        sqrt_distance = np.sqrt(tree.distance(pair[0][0], pair[1][0]))
        if sqrt_distance == 0:
            sqrt_distance = 10**-5  # consider ignoring those?
        diff = np.abs(pair[0][1] - pair[1][1])

        c = diff / sqrt_distance
        all_contrasts.append(c)
    return all_contrasts


def randomization_test(tree, mapping, n, feature):
    ''' shuffle the  tree n times and calculate variance'''

    vars = []
    # shuffle tree names along the tree n times, and calculate the contrasts variance
    for i in tqdm(range(n)):
        new_tree = tree
        terms = [term.name.split('.')[0] for term in tree.get_terminals()]   # poping from the list, need a new list each iterations
        for node in new_tree.get_terminals():
            node.name = terms.pop(random.randrange(len(terms)))

        # now we have a new tree and we can calculate variance of the contrasts
        var = np.var(contrasts(tree, feature, mapping))
        vars.append(var)

    return vars

def check_signal_significant(true_tree_variance, shuffled_variance, alias, feature, out=None):
    ''' test if the variance is significantly low. shuffled variance is a list of shuffled trees contrasts variances'''

    # need to test whether 95% of the simulations are greater then the true tree variance. look for the 5 percentile
    # and check if the true value is lower or equal

    significance = False
    shuffled_variance.append(true_tree_variance)

    cutoff = np.percentile(shuffled_variance, 5)
    print(cutoff, true_tree_variance)
    if true_tree_variance <= cutoff:
        significance = True

    if out:
        sns.kdeplot(np.asarray(shuffled_variance), shade=True, color='#378F9E')
        plt.axvline(cutoff, color='r', linestyle='--', alpha=.5)
        plt.axvline(true_tree_variance, color='g', linestyle='-', alpha=.75)
        plt.xlabel('Contrasts variance')
        plt.title('{}\nRandomization test distribution result, tree significant = {}'.format(alias,significance))
        plt.savefig(os.path.join(out, 'randomization_test_{}.pdf'.format(feature)), format='pdf', bbox_inches='tight', dpi=400)
        plt.gcf().clear()

    else:
        sns.kdeplot(np.asarray(shuffled_variance), shade=True, color='#378F9E')
        plt.axvline(cutoff, color='r', linestyle='--', alpha=.5)
        plt.axvline(true_tree_variance, color='g', linestyle='-', alpha=.75)
        plt.xlabel('Contrasts variance')
        plt.title('{}\nRandomization test distribution result, tree significant = {}'.format(alias, significance))
        plt.show()
        plt.clf()
        plt.cla()
        plt.close()

    return significance


def run_randomization_test(super_folder, n, mapping, feature):

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
        t_var = np.var(contrasts(t, feature, mapping))
        random_vars = randomization_test(t, mapping, n, feature)
        out = os.path.dirname(tree)
        sgnf = check_signal_significant(t_var, random_vars, alias, feature, out)
        df = pd.DataFrame({'group': alias, 'significance':sgnf, 'num_leafs':num_leafs}, index=[0])
        results.append(df)

    concat = pd.concat(results)
    #concat['estimation'] = concat['num_leafs'].apply(lambda x : 'underestimation' if x > 200 else 'estimated', axis=1)
    concat.to_csv(os.path.join(super_folder, 'significance_by_family_{}.csv'.format(feature)), index=False)

    return



def trait_correlations(feature1, feature2, mapping, tree):
    ''' simple linear regression for correlation testing between two features'''

    f1_contrasts = contrasts(tree, feature1, mapping)
    f2_contrasts = contrasts(tree, feature2, mapping)

    slope, intercept, r_value, p_value, std_err = stats.linregress(f1_contrasts, f2_contrasts)

    return r_value, p_value

def run_trait_correlations(featurs, mapping, super_folder, out):
    ''' run correlations between '''

    pairs = list(combinations(featurs, 2))
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

        for pair in pairs:
            r, p = trait_correlations(pair[0], pair[1], mapping, t)
            df = pd.DataFrame({'family': alias, 'attr1':pair[0], 'attr2':pair[1], 'p_value':p, 'r_sqrd':r,
                               'num_leafs':num_leafs}, index=[0])
            results.append(df)

    concat = pd.concat(results)
    concat.to_csv(os.path.join(out, 'contrasts_correlations.csv'), index=False)
    return concat


def main(args):

    main_dir = r'/sternadi/home/volume1/daniellem1/Entropy/data/Phylogeny/family'
    entropies = r'/sternadi/home/volume1/daniellem1/Entropy/data/entropies.csv'

    out = r'/sternadi/home/volume1/daniellem1/Entropy/data'

    mapping = pd.read_csv(entropies)

    n = 1000
    if args.type == 1:
        feature = 'k5'
        run_randomization_test(main_dir, n=n , mapping=mapping, feature=feature)
        print('Done!')
    else:
        features = [f for f in mapping.columns if f not in ['family', 'refseq_id', 'virus_name']]
        run_trait_correlations(features, mapping=mapping, super_folder=main_dir, out=out)
        print('Done!')



if __name__ == "__main__":
    parser = argparse.ArgumentParser()
    parser.add_argument("-t", "--type", type=int,
                        help="type of running - randomization test = 1, correlation test = 2", default=1)

    args = parser.parse_args()

    main(args)
