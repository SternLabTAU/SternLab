import pandas as pd
import numpy as np
# import matplotlib.pyplot as plt
# import seaborn as sns
from utils import joint_entropy, get_reverse_complement
from entropy_selection import string_by_codon_position
from tqdm import tqdm
import re
import os

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

def get_kmers_distribution(fasta, k, out):
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
        # sns.distplot(values, hist=False, kde_kws={'shade':True})
        plt.hist(values, alpha=0.8)

    plt.title('Distribution of kmers {}'.format(alias), fontsize=18)
    plt.xlabel('# kmers appearances', fontsize=18)
    plt.ylabel('Count', fontsize=18)
    sns.despine(offset=10)
    plt.savefig(os.path.join(out,'{}_kmers_distribution_hist.png'.format(alias)), format='png', dpi=400, bbox_inches='tight')
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
    #     plt.plot(entropies)
    # plt.title('Entropy profile {}'.format(alias), fontsize=18)
    # plt.xlabel('Genome position', fontsize=18)
    # plt.ylabel('Entropy', fontsize=18)
    # sns.despine(offset=10)
    # plt.savefig(os.path.join(out, '{}_profile.png'.format(alias)), format='png', dpi=400,
    #             bbox_inches='tight')
    # plt.gcf().clear()

    df = pd.DataFrame(dict([(k,pd.Series(v)) for k,v in all_entropies.items()]))
    df.to_csv(os.path.join(out, '{}_profile.csv'.format(alias)), index=False)

    return df











