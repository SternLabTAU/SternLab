import pandas as pd
import numpy as np
import matplotlib.pyplot as plt
import seaborn as sns
from utils import joint_entropy, get_reverse_complement, entropy_by_kmer
from entropy_selection import string_by_codon_position
from tqdm import tqdm
import re
import os
import scipy.stats as stats
import statsmodels.stats.multitest as multi
import RNA
from Bio import SeqIO
from random import sample
import glob
from functools import reduce

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


    for rec in SeqIO.parse(fasta, "fasta"):
        # get identifier and genomic sequence
        genome = rec.seq
        rc_genome = str(get_reverse_complement(genome))

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
                kmer = genome[i:i + 3]
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
            codon_trimmed = string_by_codon_position(genome, 2)
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



def get_joint_entropy_profile(fasta, w, out=None):
    """
    sliding window entropy profile of all sequences in a family
    :param fasta: a fasta file contatining viral sequences
    :param w: the window size
    :param out: optional. if != None a profile will be saved as a png
    :return: the vector of profile entropy
    """
    all_entropies = {}
    alias = os.path.basename(fasta).split('.')[0]

    i = 0
    for rec in SeqIO.parse(fasta, "fasta"):
        entropies = []
        # get identifier and genomic sequence

        genome = str(rec.seq)

        for j in range(len(genome) - w):
            sub_genome = genome[j:j+w]
            rc_sub_genome = str(get_reverse_complement(sub_genome))
            entropy = joint_entropy(sub_genome, rc_sub_genome, 5)
            entropies.append(entropy)

        print('Done with seq {}'.format(i))
        all_entropies['seq_{}'.format(i)] = entropies
        i += 1

    df = pd.DataFrame(dict([(k,pd.Series(v)) for k,v in all_entropies.items()]))

    df.to_csv(os.path.join(out, '{}_Joint_profile.csv'.format(alias)), index=False)

    return df



def get_entropy_profile(fasta, w, out=None, type='fasta'):
    """
    sliding window entropy profile of all sequences in a family
    :param fasta: a fasta file contatining viral sequences
    :param w: the window size
    :param out: optional. if != None a profile will be saved as a png
    :return: the vector of profile entropy
    """
    all_entropies = {}
    alias = os.path.basename(fasta).split('.')[0]

    i = 0
    for rec in SeqIO.parse(fasta, type):
        entropies = []
        # get identifier and genomic sequence

        genome = str(rec.seq)

        for j in range(len(genome) - w):
            sub_genome = genome[j:j+w]
            entropy = entropy_by_kmer(sub_genome, 5)
            entropies.append(entropy)


        print('Done with seq {}'.format(i))
        all_entropies['seq_{}'.format(i)] = entropies
        i += 1

    df = pd.DataFrame(dict([(k,pd.Series(v)) for k,v in all_entropies.items()]))

    df.to_csv(os.path.join(out, '{}_profile.csv'.format(alias)), index=False)

    return df

def get_joint_entropy_profile_per_sequence(seq, w, alias, out=None):
    """
    sliding window entropy profile of all sequences in a family
    :param fasta: a fasta file contatining viral sequences
    :param w: the window size
    :param out: optional. if != None a profile will be saved as a png
    :return: the vector of profile entropy
    """
    all_entropies = {}

    entropies = []
    # get identifier and genomic sequence
    genome = seq

    for j in range(len(genome) - w):
        sub_genome = genome[j:j+w]
        try:
            rc_sub_genome = get_reverse_complement(sub_genome)
            entropy = joint_entropy(sub_genome, rc_sub_genome, 5)
            entropies.append(entropy)
        except:
            break

    df = pd.DataFrame({'{}'.format(alias):entropies})
    if out != None:
        df.to_csv(os.path.join(out, '{}_profile.csv'.format(alias)), index=False)

    return df

def get_entropy_profile_per_sequence(seq, w, alias, out=None):
    """
    sliding window entropy profile of all sequences in a family
    :param fasta: a fasta file contatining viral sequences
    :param w: the window size
    :param out: optional. if != None a profile will be saved as a png
    :return: the vector of profile entropy
    """
    all_entropies = {}

    entropies = []
    # get identifier and genomic sequence
    genome = seq

    for j in range(len(genome) - w):
        sub_genome = genome[j:j+w]
        entropy = entropy_by_kmer(sub_genome, 5)
        entropies.append(entropy)


    df = pd.DataFrame({'{}'.format(alias):entropies})
    if out != None:
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

def deltaG_profile(fasta, w, out=None, type='fasta'):

    """
    sliding window free energy profile of all sequences in a family
    :param fasta: a fasta file contatining viral sequences
    :param w: the window size
    :param out: output file
    :return: the vector of profile entropy
    """
    all_deltaG = {}
    alias = os.path.basename(fasta).split('.')[0]

    i = 0
    for rec in SeqIO.parse(fasta, type):
        values = []

        genome = str(rec.seq)

        for j in range(len(genome) - w):
            sub_genome = genome[j:j+w]
            mfe = deltaG_calculator(sub_genome)
            values.append(mfe)

        print('Done with seq {}'.format(i))
        all_deltaG['seq_{}'.format(i)] = values
        i += 1

    df = pd.DataFrame(dict([(k,pd.Series(v)) for k,v in all_deltaG.items()]))
    df.to_csv(os.path.join(out, '{}_deltaG_profile.csv'.format(alias)), index=False)

    return df


def regress_entropy2deltaG(entropy_filepath, deltaG_filepath, start=None, end=None, family=None):
    """
    regress two files of entropy and delta g and get the p value and significance
    :param entropy_filepath: a path to an entropy profile file
    :param deltaG_filepath: a path to an delta g profile file
    :return: p value list
    """

    entropy = pd.read_csv(entropy_filepath)
    deltag = pd.read_csv(deltaG_filepath)

    # remove empty columns in entropy
    entropy = entropy.dropna(how='all', axis=1)
    deltag = deltag[entropy.columns]

    if start!=None:

        entropy = entropy.iloc[list(range(start, end))]
        deltag = deltag.iloc[list(range(start, end))]
        entropy = entropy.dropna(how='all', axis=1)
        deltag = deltag[entropy.columns]

    r_values = []
    p_values = []


    for c in entropy.columns:
        x = entropy[c].values
        y = deltag[c].values

        data = pd.DataFrame({'x':x, 'y':y})
        data = data.dropna()
        slope, intercept, r_value, p_value, std_err = stats.linregress(data['x'], data['y'])
        r_values.append(r_value)
        p_values.append(p_value)


    result = pd.DataFrame({'R_value':r_values, 'p_value':p_values})
    result['corrected_pvalue'] = multi.fdrcorrection(result['p_value'])[1]
    if family != None:
        result['family'] = family
    return result


def enrichment_test_OU(df, feature, c):
    """
    test if there exists any enrichments in the OU families
    :param df:
    :return:
    """

    ou_feature = df[(df['model'] == 'OU') & (df[feature] == c)]
    ou_Notfeature = df[(df['model'] == 'OU') & (df[feature] != c)]
    bm_feature = df[(df['model'] == 'BM') & (df[feature] == c)]
    bm_Notfeature = df[(df['model'] == 'BM') & (df[feature] != c)]

    x11 = len(set(ou_feature['family']))
    x12 = len(set(bm_feature['family']))
    x21 = len(set(ou_Notfeature['family']))
    x22 = len(set(bm_Notfeature['family']))

    oddsratio, pvalue = stats.fisher_exact([[x11, x12], [x21, x22]])

    return oddsratio, pvalue


def sample_from_fasta(fasta, n, out, ref_id='NC_004162'):
    """
    This method randomly samples n sequences from a fasta file and saves it to a new file
    :param fasta: a fasta file containing < 300 sequences
    :param n: number of sequences to sample
    :param out: output file
    :param ref_id: an id of the reference sequence
    :return: saves the new fasta file to "out"
    """

    seqs = list(SeqIO.parse(fasta, 'fasta'))
    seqs_to_sample = []
    reference = ''

    for s in seqs:
        if s.id == ref_id:
            reference = s
        else:
            seqs_to_sample.append(s)

    sampled = sample(seqs_to_sample, n-1)
    sample.append(reference)
    with open(out, 'w') as o:
        SeqIO.write(sampled, o, "fasta")

    return


def parse_r4s_results(data, out=None):
    """
    parse rate 4 site output file
    :param data: a data output file
    :param out: output filepath
    :return: a data frame containing the r4s scores per position
    """

    # read the data frame and skip the 13 first rows (comments)
    df = pd.read_table(data, skiprows=13, header=None)

    pos = []
    seq = []
    score = []
    upper = []
    lower = []
    std = []

    for i in range(df.shape[0]):
        row = df.iloc[i,0].split()
        # we are at the end of the file, two comments not records
        if len(row) < 6 :
            break

        # get all floats from the row
        p = re.compile(r'\d+\.*\d*')
        floats = [float(i) for i in p.findall(df.iloc[i, 0])]    # pos, score, lower, upper, std
        pos.append(int(row[0]))
        seq.append(row[1])
        score.append(floats[1])
        lower.append(floats[2])
        upper.append(floats[3])
        std.append(floats[4])


    res = pd.DataFrame({'pos':pos, 'seq':seq, 'score':score, 'lower':lower, 'upper':upper, 'std':std})

    # calculate the confidence in each positions
    scores = np.sort(res['score'].values)
    bins = np.array([elem[0] for elem in np.split(scores, 10)]) # get all bins
    lower_bins = np.digitize(res['lower'].values, bins)
    upper_bins = np.digitize(res['upper'].values, bins)
    diff = np.abs(upper_bins - lower_bins)
    score_bins = np.digitize(res['score'], bins)

    # update data frame
    res['lower_bin'] = lower_bins
    res['upper_bin'] = upper_bins
    res['diff'] = diff
    res['score_bin'] = score_bins

    res['confidence'] = res.apply(lambda row: row['score_bin'] if row['diff'] <= 4 else 20, axis=1)


    if out != None:
        res.to_csv(os.path.join(os.path.dirname(data), os.path.basename(data)+'_parsed.csv'), index=False)
    return res


def r4s_2_heatmap(filepath, addition = None):
    """
    create a data frame with position to conservation by all alignments sizes
    :param filepath: an input directory with parsed csv files from r4s
    :return: a data frame with position and the averaged score by the alignment
    """

    all_r4s = glob.glob(os.path.join(filepath, '*.csv'))

    if addition != None:
        all_r4s = glob.glob(os.path.join(filepath, '*.csv')) + glob.glob(os.path.join(addition, '*.csv'))
    num_seqs = set([os.path.basename(x).split('_')[1] for x in all_r4s])

    all_means = []

    for n in num_seqs:
        csvs = glob.glob(os.path.join(filepath, 'n_{}_*.csv'.format(n)))
        if addition != None:
            csvs = glob.glob(os.path.join(filepath, 'n_{}_*.csv'.format(n))) + glob.glob(os.path.join(addition, 'n_{}_*.csv'.format(n)))
        dfs = []
        for c in csvs:
            df = pd.read_csv(c)
            df = df[['pos','confidence']]
            dfs.append(df)

        df_merged = reduce(lambda left, right: pd.merge(left, right, on=['pos'], how='outer'), dfs)
        df_merged = df_merged.dropna()
        df_merged.set_index(['pos'], inplace=True)

        df_merged[n] = df_merged.median(axis=1)

        df_merged = df_merged[[n]]
        all_means.append(df_merged)

    result = pd.concat(all_means, axis=1)
    return result
