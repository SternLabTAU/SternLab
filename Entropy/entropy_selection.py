import matplotlib
matplotlib.use('Agg')
import pandas as pd
import numpy as np
import os
from utils import *
from randomization_utils import *
import re
from tqdm import tqdm
import math
from functools import reduce
import string
from Bio import SeqIO, Phylo
from itertools import combinations
from scipy import stats
import matplotlib.pyplot as plt
import seaborn as sns
import heritability
import pickle
from randomization_utils import *

# with open(r'/sternadi/home/volume1/daniellem1/Entropy/data/refseq_2_cds_positions.pickle', 'rb') as o:
#     cds_mapping = pickle.load(o)

def remove_punctuation(s):
    splitted = string.punctuation.split('_')
    marks = splitted[0]+splitted[1]

    regex = re.compile('[%s]' % re.escape(marks))
    return regex.sub('', s).strip()


def map_family_2_refseq(x, mapping):
     result = [k for k, v in mapping.items() if x in v]
     if result == []:
        return None
     else:
        return result[0]



def genome_2_entropy(fasta, k, scramble=False, rc_joint=False, out=""):
    """
    create a mapping of the complete genome to its entropy by kmer
    :param fasta: a file containing all fasta files
    :param k: the k-mer size
    :return: a data frame of ref-seq id and entropy measurement
    """

    k_entropy = []
    ref_seq_id = []

    sequences = re.split(">", open(fasta, "r").read().replace('\n', ''))[1:]
    for seq in tqdm(sequences):
        if '.' not in seq:
            print('no dot in sequence name\n')
            continue
        if ('complete genome' not in seq) and ('complete sequence' not in seq) and ('complete v' not in seq) \
                and ('genome' not in seq):
            print('not complete genome\n')
            continue

        # get identifier and genomic sequence
        splitted = seq.split('.')
        identifier = remove_punctuation(splitted[0].split()[1])
        #identifier = splitted[0].split()[0]
        genome = splitted[-1]

        if scramble:
            genome = scramble_all_sequence(genome, mapping=cds_mapping, refseq_id=identifier, how=3)
            if genome == None:
                continue
        # calculate entropy

        if rc_joint:
            rc_genome = get_reverse_complement(genome)
            entropy = joint_entropy(genome, rc_genome, k)

        else:
            entropy = entropy_by_kmer(genome, k)


        k_entropy.append(entropy)
        ref_seq_id.append(identifier)

    df = pd.DataFrame({'refseq_id':ref_seq_id, 'k{}'.format(k):k_entropy})
    df.to_csv(os.path.join(out, 'k{}.csv'.format(k)), index=False)
    return df

def string_by_codon_position(seq, position):
    """
    extract only the specified position of a codon triplet
    :param seq: a string
    :param position: the position in the codon can be 0,1,2
    :return: the sting containing only the codon at position in each triplet
    """

    triplets = [seq[i:i+3] for i in range(0, len(seq), 3)]
    filtered = ''
    filtered = filtered.join([tripl[position] for tripl in triplets if len(tripl) == 3])
    return filtered

def refseq_2_cds(cds):
    """
    creates a mapping of the cds to refseq
    :param cds: a fasta files
    :return:
    """


    # create a mapping between each id all the coding sequences it has.
    mapping = {}
    sequences = re.split(">", open(cds, "r").read().replace('\n', ''))[1:]
    for seq in sequences:
        splitted = seq.split('|')
        # if splitted[2] != '':
            # print(splitted[0], splitted[2], splitted[3])
        refseq_id = remove_punctuation(splitted[3])
        if 'geneid' in refseq_id.lower():
            refseq_id = remove_punctuation(splitted[4])
        # refseq_id = splitted[3]
        if '}' in seq:
            coding = splitted[-1].split('}')[-1]
        else:
            coding = splitted[-1].split(')')[-1]

        # add the sequence to the dictionary
        if refseq_id in mapping:
            mapping[refseq_id].append(coding)
        else:
            mapping[refseq_id] = [coding]

    return mapping

    # now we have for each id the cds seqs.

def entropy_calculator_cds(key, values, jump=False, position=0, scramble=True, rc_joint=False):
    '''
    calculate the entropy according to coding regions only for a given set of . here we have two options - by reading frame or by
    some codon position
    :param key: refseq id of the sequences in values
    :param values: the cds sequences of that entry
    :param jump: reading frame indication
    :param position: the codon position
    :return: entropy
    '''

    # for each sequence get the entropy by reading frame, by k=5, and by k=1 according to each position.

    kmers = {}
    entropy = 0
    total_len = 0

    if not rc_joint:
        for seq in values:
            if scramble:
                seq = scramble_synonymous_codons(seq)
            total_len += len(seq)
            if jump:

                assert(position in [0,1,2])
                for i in range(position, len(seq) - 3, 3):
                    kmer = seq[i:i + 3]
                    if kmer in kmers:
                        kmers[kmer] += 1
                    else:
                        kmers[kmer] = 1
            else:
                assert(position != None)

                codon_trimmed = string_by_codon_position(seq, position)
                for i in range(len(codon_trimmed)):
                    kmer = codon_trimmed[i]
                    if kmer in kmers:
                        kmers[kmer] += 1
                    else:
                        kmers[kmer] = 1

        # calculate entropy for all the cds of this refseq id
        total_kmers = sum(kmers.values())
        for kmer in kmers:
            p = kmers[kmer] / total_kmers
            entropy += -(p * math.log2(p))


    if rc_joint:
        kmers1 = {}
        kmers2 = {}

        for seq in values:
            if scramble:
                seq = scramble_synonymous_codons(seq)
            total_len += len(seq)
            if jump:
                assert (position in [0, 1, 2])

                rc_seq = get_reverse_complement(seq)

                # calculate k-mers for each sequence, and combine the distributions
                # use k=3 for reading frames
                for i in range(position, len(seq) - 3, 3):
                    kmer = seq[i:i + 3]
                    if kmer in kmers1:
                        kmers1[kmer] += 1
                    else:
                        kmers1[kmer] = 1

                for i in range(position, len(rc_seq) - 3, 3):
                    kmer = rc_seq[i:i + 3]
                    if kmer in kmers2:
                        kmers2[kmer] += 1
                    else:
                        kmers2[kmer] = 1


            else:
                assert (position != None)

                # use k=1 in this case by looking at a specific codon position in the cdc
                codon_trimmed = string_by_codon_position(seq, position)

                rc_seq = get_reverse_complement(seq)
                rc_codon_trimmed = string_by_codon_position(rc_seq, position)
                for i in range(len(codon_trimmed)):
                    kmer = codon_trimmed[i]
                    if kmer in kmers1:
                        kmers1[kmer] += 1
                    else:
                        kmers1[kmer] = 1

                for i in range(len(rc_codon_trimmed)):
                    kmer = codon_trimmed[i]
                    if kmer in kmers2:
                        kmers2[kmer] += 1
                    else:
                        kmers2[kmer] = 1

        # calculate entropy for all the cds of this refseq id
        total_kmers_1 = sum(kmers1.values())
        total_kmers_2 = sum(kmers2.values())

        total = total_kmers_1 + total_kmers_2

        # compare the kmers space to be equal at both
        for kmer in kmers1:
            if kmer not in kmers2:
                kmers2[kmer] = 0

        for kmer in kmers2:
            if kmer not in kmers1:
                kmers2[kmer] = 0

        joint_entropy = 0
        for kmer1 in kmers1:
            for kmer2 in kmers2:
                p_xy = (kmers1[kmer1] + kmers2[kmer2]) / total

                joint_entropy += -(p_xy * math.log2(p_xy))


        entropy = joint_entropy

    return entropy


def entropy_by_cds(mapping, out=""):
    """
    create a data frame with all entropy values for a given refseq id
    :param mapping: dictionary with refseq id and cds sequences
    :return: a data frame with all information.
    """
    result = []

    for key in mapping.keys():
        rf_entropy = entropy_calculator_cds(key, mapping[key], jump=True, position=0, rc_joint=True)
        rf2_entropy = entropy_calculator_cds(key, mapping[key], jump=True, position=1, rc_joint=True)
        rf3_entropy = entropy_calculator_cds(key, mapping[key], jump=True, position=2, rc_joint=True)
        codon1_entropy = entropy_calculator_cds(key, mapping[key], jump=False, position=0, rc_joint=True)
        codon2_entropy = entropy_calculator_cds(key, mapping[key], jump=False, position=1, rc_joint=True)
        codon3_entropy = entropy_calculator_cds(key, mapping[key], jump=False, position=2, rc_joint=True)
        df = pd.DataFrame({'refseq_id':key, 'reading_frame':rf_entropy, 'reading_frame_shift_1':rf2_entropy,
                           'reading_frame_shift_2':rf3_entropy, 'codon_position_1':codon1_entropy,
                           'codon_position_2': codon2_entropy, 'codon_position_3':codon3_entropy}, index=[0])
        result.append(df)

    all_values = pd.concat(result)
    all_values.to_csv(os.path.join(out, 'cds_entropies.csv'), index=False)
    return all_values



def tree_slope_calculator(tree, mapping, feature, alias):

    value_2_tree_distance = []

    term_names = [term.name.split('.')[0] for term in tree.get_terminals()]
    values = [(name, mapping[mapping['refseq_id'] == name][feature].values[0]) for name in term_names]

    # calculate all pairs
    pairs =  list(combinations(values, 2))
    for pair in pairs:
        if np.isnan(pair[0][1]) or np.isnan(pair[1][1]):
            # print("we have None for feature {}".format(feature))
            continue
        distance = np.sqrt(tree.distance(pair[0][0], pair[1][0]))
        diff = np.abs(pair[0][1] - pair[1][1])

        value_2_tree_distance.append((diff, distance))

    # linear regression for the coefficients
    #x = [np.log(val[1]) for val in value_2_tree_distance]
    x = [val[1] for val in value_2_tree_distance]
    y = stats.boxcox(np.asarray([val[0] + 0.00001 for val in value_2_tree_distance]))[0]
    slope, intercept, r_value, p_value, std_err = stats.linregress(x, y)

    n = len(x)

    # get the estimator for the slope

    t = stats.t.ppf(0.975, n-2)

    estimator = t * std_err
    lower_bound = slope - estimator
    upper_bound = slope + estimator


    df = pd.DataFrame({'feature':feature, 'family':alias, 'slope':slope, 'intercept':intercept,
                       'r_value':r_value, 'p_value':p_value, 'lower_CI':lower_bound, 'upper_CI':upper_bound},
                      index=[0])
    return df

def get_family_slopes_and_graphs(tree, mapping, features, alias, out=None):
    sns.set_style('white')
    all_df = []
    for feature in features:
        df = tree_slope_calculator(tree, mapping, feature, alias)
        all_df.append(df)

    features_df = pd.concat(all_df)


    # create the plot
    x = np.linspace(0,1,50)
    #colors = ['#E9C025', '#63AB2B', '#2BA5AB', '#9C59DC', '#E9ADE6', '#961232', '#D1357C', '#36DCAA', '#DC6B36']
    colors = list(sns.color_palette('Set2', len(features)).as_hex())
    for i, feature in enumerate(features):
        feature_info = features_df[features_df['feature'] == feature]
        slope = feature_info['slope'][0]    # no need in values because its only one line
        intercept = feature_info['intercept'][0]
        y = [slope * i for i in x]
        plt.plot(x, y, '--', label=feature, color=colors[i])

    plt.xlabel('Evolutionary time', fontsize = 16)
    plt.ylabel('Entropy value difference', fontsize = 16)
    plt.title('Entropy difference as a factor of evolutionary time for {} family'.format(alias), fontsize=22)
    plt.legend()
    sns.despine()
    if out != None:
        plt.savefig(os.path.join(out, 'selection_test.png'), format='png', dpi=400, bbox_inches='tight')
        plt.gcf().clear()

    return features_df

def run_entropy_selection_test(super_folder, mapping, out):

    all_trees = []
    results = []

    #features = [f for f in mapping.columns if f not in ['family', 'refseq_id', 'virus_name']]
    # features = ['k5','reading_frame', 'codon_position_1', 'codon_position_2', 'codon_position_3']
    features = ['k5','reading_frame', 'first codon position', 'second codon position', 'third codon position']

    do_not_consider = ['Avsunviroidae', 'Pospiviroidae', 'Deltasatellite']

    for root, dirs, files in tqdm(os.walk(super_folder)):
        tree = [f for f in files if 'phyml_tree' in f]
        if tree != []:
            all_trees.append(os.path.join(root, tree[0]))

    for t in tqdm(all_trees):
        if heritability.tree_2_string(t) == '':
            continue
        alias = os.path.basename(t).split('.')[0].strip()
        if alias in do_not_consider:
            continue
        print(alias)
        tree = Phylo.read(t, 'newick')
        df = get_family_slopes_and_graphs(tree, mapping, features, alias, os.path.dirname(t))
        results.append(df)

    final = pd.concat(results)
    final.to_csv(os.path.join(out, 'entropy_selection_test.csv'), index=False)


def cds_statistics_by_family(tree, family, mapping, features):
    """
    calculate the number of cds sequences per each leaf in the tree
    :param tree: phylo tree object
    :param family: a string indicating the family name
    :param mapping: a data frame
    :return: a data frame with statistics for each leaf
    """

    features_dict = {}
    term_names = [term.name.split('.')[0] for term in tree.get_terminals()]
    num_leafs = len(term_names)
    for feature in features:
        print(family, feature)
        print(term_names)
        values = [(name, mapping[mapping['refseq_id'] == name][feature].values[0]) for name in term_names]
        num_not_none_leafs = len([val2 for val1, val2, in values if not np.isnan(val2)])
        features_dict[feature] = num_not_none_leafs/num_leafs

    df = pd.DataFrame({'family':family, 'k5':features_dict['k5'], 'codon_position_1':features_dict['codon1_entropy']
                          , 'codon_position_2':features_dict['codon2_entropy'], 'codon_position_3':features_dict['codon3_entropy'],
                       'reading_frame':features_dict['rf_entropy'], 'num_leafs':num_leafs}, index=[0])

    return df

def run_cds_statistics_by_family(super_folder, mapping, out):

    all_trees = []
    results = []

    features = [f for f in mapping.columns if f not in ['family', 'refseq_id', 'virus_name']]


    for root, dirs, files in tqdm(os.walk(super_folder)):
        tree = [f for f in files if 'phyml_tree' in f]
        if tree != []:
            all_trees.append(os.path.join(root, tree[0]))

    for t in tqdm(all_trees):
        if heritability.tree_2_string(t) == '':
            continue
        alias = os.path.basename(t).split('.')[0].strip()
        print(alias)
        tree = Phylo.read(t, 'newick')
        df = cds_statistics_by_family(tree, alias, mapping, features)
        results.append(df)

    final = pd.concat(results)
    final.to_csv(os.path.join(out, 'genomic_cds_entropy_stats.csv'), index=False)

def test_selection_validity(slopes,feature1,feature2, out=None):
    """
    tests the selection of each family by codon3 entropy vs. entropy k=5 using ci and hypothesis testing
    :param slopes: a dataframe with slopes
    :param feature1: the wanted upper slope
    :param feature2: the wanted lower slope
    :param families:
    :return:
    """

    mapping = {}
    features = [feature1, feature2]
    slopes = slopes[slopes['feature'].isin(features)]
    for family in tqdm(set(slopes['family'])):
        df = slopes[slopes['family'] == family]
        is_significant = df[df['feature'] == feature1]['p_value'].values[0] < 0.05 and \
                         df[df['feature'] == feature2]['p_value'].values[0] < 0.05

        differ = df[df['feature'] == feature1]['lower_CI'].values[0] > \
                 df[df['feature'] == feature2]['upper_CI'].values[0]

        mapping[family] = (is_significant, differ)

    slopes['pval_significant'] = slopes['family'].apply(lambda x: mapping[x][0])
    slopes['slope_difference'] = slopes['family'].apply(lambda x: mapping[x][1])
    if out != None:
        slopes.to_csv(os.path.join(out, 'selection_test_stats_{}_{}.csv'.format(feature1, feature2)), index=False)

    return slopes



#
#
# def main():
#     cds = r'/sternadi/home/volume1/daniellem1/Entropy/data/virushostdb.cds.fna'
#     genomics = r'/sternadi/home/volume1/daniellem1/Entropy/data/virushostdb.genomic.fna'
#
#
#     out = r'/sternadi/home/volume1/daniellem1/Entropy/data/'
#
#     dfs = []
#
#     for k in range(1,11):
#         basic_genome_entropy = genome_2_entropy(genomics, k, rc_joint=True)
#         dfs.append(basic_genome_entropy)
#
#     mapping = refseq_2_cds(cds)
#     cds_mapping = entropy_by_cds(mapping)
#     dfs.append(cds_mapping)
#
    # merged = reduce(lambda left, right: pd.merge(left, right, on=['refseq_id'],
    #                                              how='outer'), dfs)
    #
    # with open(r'/sternadi/home/volume1/daniellem1/Entropy/data/family_2_refseq.pickle', 'rb') as f:
    #     d = pickle.load(f)
    #
    # merged['family'] = merged['refseq_id'].apply(lambda x : map_family_2_refseq(x, d))
#
#     merged.to_csv(os.path.join(out, 'entropies.csv'), index=False)
#
#
#     # entropy_mapping = pd.read_csv(ent_map)
#     # entropy_mapping.dropna(subset=['entropy_5'], inplace=True)
#     # run_cds_statistics_by_family(super_fold, entropy_mapping, out)
#
# if __name__ == '__main__':
#     main()



