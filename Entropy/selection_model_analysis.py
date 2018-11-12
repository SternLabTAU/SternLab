import pandas as pd
import numpy as np
import matplotlib.pyplot as plt
import seaborn as sns
from utils import joint_entropy, get_reverse_complement
from entropy_selection import string_by_codon_position
from tqdm import tqdm
import re
import os
import scipy.stats as stats
import statsmodels.stats.multitest as multi
import RNA

cp2_path = r'/sternadi/home/volume1/daniellem1/Entropy/data/OU_model/simulations_significance_bm_cp2.csv'
rf_path = r'/sternadi/home/volume1/daniellem1/Entropy/data/OU_model/simulations_significance_bm_rf.csv'
k5_path = r'/sternadi/home/volume1/daniellem1/Entropy/data/OU_model/simulations_significance_bm_k5.csv'

# # read files
# cp2 = pd.read_csv(cp2_path)
# rf = pd.read_csv(rf_path)
# k5 = pd.read_csv(k5_path)
#
#
# # create basic tables and print them to screen
# # overall false bm
# print('BM simulations')
# print('codon position 2')
# print(cp2['isBM'].value_counts()/cp2.shape[0])
# print('reading frame')
# print(rf['isBM'].value_counts()/rf.shape[0])
# print('k5')
# print(k5['isBM'].value_counts()/k5.shape[0])
#
#
# # plot distributions
# # overall false bm
# print('BM simulations')
# print('codon position 2')
# print(cp2['model'].value_counts()/cp2.shape[0])
# print('reading frame')
# print(rf['model'].value_counts()/rf.shape[0])
# print('k5')
# print(k5['model'].value_counts()/k5.shape[0])


#### auxiliary functions for OU significant families ####

def get_kmers_distribution(fasta, k, out=None):
    """
    get the kmers distribution plot for each family separately
    :param fasta: fasta file
    :param k: the kmer length
    :return: saves the plot
    """

    alias = os.path.basename(fasta).split('.')[0]
    all_values = []

    sequences = re.split(">", open(fasta, "r").read().replace('\n', ''))[1:]
    for seq in sequences:
        # get identifier and genomic sequence
        splitted = seq.split('.')
        genome = splitted[-1]
        #rc_genome = get_reverse_complement(genome)

        kmers_1 = {}
        kmers_2 = {}

        if k == 5:
            # sliding window of k
            for i in range(len(genome) - k):
                kmer = genome[i:i+k]
                if kmer in kmers_1:
                    kmers_1[kmer] += 1
                else:
                    kmers_1[kmer] = 1

            # for i in range(len(rc_genome) - k):
            #     kmer = rc_genome[i:i+k]
            #     if kmer in kmers_2:
            #         kmers_2[kmer] += 1
            #     else:
            #         kmers_2[kmer] = 1

        elif k == 3:
            # reading frame
            for i in range(0, len(genome) - 3, 3):
                kmer = seq[i:i + 3]
                if kmer in kmers_1:
                    kmers_1[kmer] += 1
                else:
                    kmers_1[kmer] = 1

            # for i in range(0, len(rc_genome) - 3, 3):
            #     kmer = rc_genome[i:i + 3]
            #     if kmer in kmers_2:
            #         kmers_2[kmer] += 1
            #     else:
            #         kmers_2[kmer] = 1
        else:
            assert(k==1)
            codon_trimmed = string_by_codon_position(seq, 2)
            # rc_codon_trimmed = get_reverse_complement(seq)
            for i in range(len(codon_trimmed)):
                kmer = codon_trimmed[i]
                if kmer in kmers_1:
                    kmers_1[kmer] += 1
                else:
                    kmers_1[kmer] = 1

            # for i in range(len(rc_codon_trimmed)):
            #     kmer = rc_codon_trimmed[i]
            #     if kmer in kmers_2:
            #         kmers_2[kmer] += 1
            #     else:
            #         kmers_2[kmer] = 1


        # create one dictionary for all kmers
        all_kmers = {x: kmers_1.get(x, 0) + kmers_2.get(x, 0) for x in set(kmers_1) | set(kmers_2)}
        values = [int(x) for x in all_kmers.values()]
        all_values.append(all_kmers)
        if out != None:
            # sns.distplot(values, hist=False, kde_kws={'shade':True})
            plt.hist(values, alpha=0.8, normed=True)

    if out != None:
        plt.title('Distribution of kmers {}'.format(alias), fontsize=18)
        plt.xlabel('# kmers appearence', fontsize=18)
        plt.ylabel('Count', fontsize=18)
        sns.despine(offset=10)
        plt.savefig(os.path.join(out,'{}_kmers_distribution_hist_normed.png'.format(alias)), format='png', dpi=400, bbox_inches='tight')
        plt.gcf().clear()

    return all_values



def get_entropy_profile(fasta, w, out=None):
    """
    sliding window entropy profile of all sequences in a family
    :param fasta: a fasta file contatining viral sequences
    :param w: the window size
    :param out: optional. if != None a profile will be saved as a png
    :return: the vector of profile entropy
    """
    all_entropies = {}
    alias = os.path.basename(fasta).split('.')[0]

    sequences = re.split(">", open(fasta, "r").read().replace('\n', ''))[1:]
    for i, seq in tqdm(enumerate(sequences)):
        entropies = []
        # get identifier and genomic sequence
        splitted = seq.split('.')
        genome = splitted[-1]

        for j in range(len(genome) - w):
            sub_genome = genome[j:j+w]
            try:
                rc_sub_genome = get_reverse_complement(sub_genome)
                entropy = joint_entropy(sub_genome, rc_sub_genome, 5)
                entropies.append(entropy)
            except:
                break

        print('Done with seq {}'.format(i))
        all_entropies['seq_{}'.format(i)] = entropies
        plt.plot(entropies)
    plt.title('Entropy profile {}'.format(alias), fontsize=18)
    plt.xlabel('Genome position', fontsize=18)
    plt.ylabel('Entropy', fontsize=18)
    sns.despine(offset=10)
    plt.savefig(os.path.join(out, '{}_profile.png'.format(alias)), format='png', dpi=400,
                bbox_inches='tight')
    plt.gcf().clear()

    df = pd.DataFrame(dict([(k,pd.Series(v)) for k,v in all_entropies.items()]))
    df.to_csv(os.path.join(out, '{}_profile.csv'.format(alias)), index=False)

    return df




def construct_heatmap(filepath, out):
    """
    create a heatmap form entropy profiles
    :param filepath: a path to a csv file, containing sequences as columns and genome positions as rows
    :param out: output file path to save the plot
    :return: 
    """
    df = pd.read_csv(filepath)
    alias = os.path.basename(filepath).split('_')[0]
    ax = sns.heatmap(df, xticklabels=False, yticklabels=False)
    plt.title('{} entropy profile'.format(alias), fontsize=20)
    plt.savefig(os.path.join(out, '{}_profile.png'.format(alias)), format='png', dpi=400,
                bbox_inches='tight')
    plt.gcf().clear()



def construct_clustermap(filepath, out):
    """
    create a heatmap form entropy profiles
    :param filepath: a path to a csv file, containing sequences as columns and genome positions as rows
    :param out: output file path to save the plot
    :return:
    """
    df = pd.read_csv(filepath)
    df = df.dropna()
    alias = os.path.basename(filepath).split('_')[0]
    ax = sns.clustermap(df, row_cluster=False, yticklabels=False)
    plt.savefig(os.path.join(out, 'clustered_{}_profile.png'.format(alias)), format='png', dpi=400,
                bbox_inches='tight')
    plt.gcf().clear()


def get_kmers_statistics(fasta):
    """
    get the most representable kmers for a certain fasta file
    :param fasta: a fasta file containing dna sequences
    :return: a list of kmers and thier relative appearance in the the kmers distribution
    """

    alias = os.path.basename(fasta).split('.')[0]
    count_by_kmer = get_kmers_distribution(fasta, 5)
    total_count = sum([sum(d.values()) for d in count_by_kmer]) # all kmers counts across al sequences in a family

    results = {}


    for i, d in enumerate(count_by_kmer):
        # get the cutoff of the kmers distribution as the highest 1%
        top_1_prect = np.percentile(list(d.values()), 99)
        sorted_counts = sorted(d.items(), key = lambda t: t[1], reverse = True)
        top_kmers = [(key, value) for key, value in sorted_counts if value >= top_1_prect]

        results['seq_{}'.format(i)] = top_kmers

    # get all kmers counts in the top list
    top_kmers_2_count = {}
    L = [dict(l) for l in list(results.values())]

    for d in L:
        top_kmers_2_count.update(d)

    dfs = []
    for kmer in top_kmers_2_count:
        iskmer_isTop1 = sum([sum([v for k, v in d if k == kmer]) for d in results.values()])
        isKmer = sum([sum([v for k,v in d.items() if k == kmer]) for d in count_by_kmer])
        isKmer_notTop1 = isKmer - iskmer_isTop1
        isTop1 = sum([sum([v[1] for v in d]) for d in results.values()])
        isTop1_notKmer = isTop1 - iskmer_isTop1
        notKmer_notTop1 = total_count - iskmer_isTop1 - isTop1_notKmer - isKmer_notTop1


        # fisher test
        oddsratio, pvalue = stats.fisher_exact([[iskmer_isTop1, isTop1_notKmer],[isKmer_notTop1,notKmer_notTop1]])
        df = pd.DataFrame({'fisher':pvalue, 'family':alias, 'kmerID':kmer}, index=[0])
        dfs.append(df)

    final_df = pd.concat(dfs)
    final_df['corrected_pvalue'] = multi.fdrcorrection(final_df['fisher'])[1]

    return final_df



def get_n_moment(data, n):
    """
    returns an array with first n momrnts of a dataset
    :param data: array of integers
    :param n: number of moments
    :return: a list of moments corresponding to data
    """

    moments = []
    for i in range(1, n+1):
        moments.append(stats.moment(data,i))

    return moments


def deltaG_calculator(seq):
    """
    calculate the minimum free energy (G) for a given sequence
    :param seq: an rna sequence
    :return: minimum free energy
    """
    ss, mfe = RNA.fold(seq)
    return mfe

def deltaG_profile(fasta, w, out):

    """
    sliding window free energy profile of all sequences in a family
    :param fasta: a fasta file contatining viral sequences
    :param w: the window size
    :param out: output file
    :return: the vector of profile entropy
    """
    all_deltaG = {}
    alias = os.path.basename(fasta).split('.')[0]

    sequences = re.split(">", open(fasta, "r").read().replace('\n', ''))[1:]
    for i, seq in tqdm(enumerate(sequences)):
        values = []
        # get identifier and genomic sequence
        splitted = seq.split('.')
        genome = splitted[-1]
        rna_genome = genome.replace('t', 'u')

        for j in range(len(rna_genome) - w):
            sub_genome = rna_genome[j:j+w]
            mfe = deltaG_calculator(sub_genome)
            values.append(mfe)

        print('Done with seq {}'.format(i))
        all_deltaG['seq_{}'.format(i)] = values

    df = pd.DataFrame(dict([(k,pd.Series(v)) for k,v in all_deltaG.items()]))
    df.to_csv(os.path.join(out, '{}_deltaG_profile.csv'.format(alias)), index=False)

    return df







