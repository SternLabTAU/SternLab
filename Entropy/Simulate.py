import pandas as pd
import numpy as np
import seaborn as sns
import matplotlib.pyplot as plt
from tqdm import tqdm
from utils import joint_entropy, entropy_by_kmer, get_reverse_complement
import os
import glob
import argparse
from selection_model_analysis import get_entropy_profile, get_joint_entropy_profile
from stats_utils import find_sequential_drops, stretchFinder
import re

def simulate_genome_by_composition(p, n, size, mode):
    """
    simulate genomes of changing nucleotide compositions
    :param p: the proportions of each character
    :param n: length of simulated sequence
    :param size: number of sequences to simulate
    :param mode: the mode of simulation: 1= no structure, 2= structure
    :return: sequences and corresponding names
    """

    sequences = []
    names = []
    if mode == 1:
        # repetitive sequences
        for i in tqdm(range(size)):
            seq = ''.join(np.random.choice(['a', 'c', 'g', 't'], p=p, size=n))
            sequences.append(seq)
            names.append('mode_{}_seq_{}'.format(mode,i))

    else:
        # only structure - generate a perfect stem loop
        for i in tqdm(range(size)):
            seq = ''.join(np.random.choice(['a', 'c', 'g', 't'], p=p, size=n // 2))
            seq = seq + 'aaaaaa' + str(get_reverse_complement(seq))
            sequences.append(seq)
            names.append('mode_{}_seq_{}'.format(mode, i))

    return sequences, names

def simulate_genome_by_drops(size, w, genome_size=5000):
    """
    simulate genomes of changing nucleotide compositions
    :param n: length of simulated sequence
    :param size: number of sequences to simulate
    :param w: the drop size
    :return: sequences and corresponding names
    """

    metagenome = ''.join(np.random.choice(['a', 'c', 'g', 't'], p=[0.25,0.25,0.25,0.25], size=genome_size))

    sequences = []
    names = []

    # simulate size genomes
    for i in tqdm(range(size)):
        # drops 1 - homogeneous sequence
        drop_1 = np.random.choice(['a', 'c', 'g', 't']) * w

        # drop 2 - repetitive and structure
        letter1 = np.random.choice(['a', 'c', 'g', 't'])
        letter2 = np.random.choice([x for x in ['a', 'c', 'g', 't'] if x != letter1])   # we want a different nuc.

        drop_2 = letter1 * w//2 + letter2 * w//2

        #drop 3 - pure structure
        stem_arm = ''.join(np.random.choice(['a', 'c', 'g', 't'], p=[0.25,0.25,0.25,0.25], size=w//2 - 5))
        loop = np.random.choice(['a', 'c', 'g', 't']) * 10  # loop of 10 nuc.

        drop_3 = stem_arm + loop + str(get_reverse_complement(stem_arm))

        # drop 4 - bias in nucleotide composition
        nucs = ['a', 'c', 'g', 't']
        np.random.shuffle(nucs)
        drop_4 = ''.join(np.random.choice(nucs, p=[0.6,0.2,0.1,0.1], size=w))


        # insert genomes to metagenome and save the indices.
        simulated_genome = metagenome[:1000] + drop_1 + metagenome[1000:2000] + drop_2 +\
            metagenome[2000:3000] + drop_3 + metagenome[3000:4000] + drop_4 + metagenome[4000:]

        sequences.append(simulated_genome)
        names.append('seq_{}'.format(i))

    return sequences, names





def to_fasta_file(output_path, sequences, names):
    """
    saves a fasta file of all sequences int the list
    :param output_dir: a path to save the fasta file
    :param sequences: list of sequences
    :param names: list of corresponding names
    :return: saves a fasta file
    """

    with open(output_path, 'w') as o:
        for i in range(len(sequences)):
            o.write('>' + names[i] + '\n' + sequences[i] + '\n')


def simulation_runner():
    # run all simulations
    p_rep = [0.6,0.2,0.1,0.1]
    p_random = [0.25,0.25,0.25,0.25]
    outdir = r'/Users/daniellemiller/Google Drive/Msc Bioinformatics/Projects/Entropy/Simulations/Classifications'

    seqs, names = simulate_genome_by_composition(p_rep, 10**4, 1000, 1)
    output = os.path.join(outdir, 'repetative.fasta')
    to_fasta_file(output,seqs, names)

    seqs, names = simulate_genome_by_composition(p_rep, 10**4, 1000, 2)
    output = os.path.join(outdir, 'repetative_structure.fasta')
    to_fasta_file(output,seqs, names)

    seqs, names = simulate_genome_by_composition(p_random, 10**4, 1000, 2)
    output = os.path.join(outdir, 'structure.fasta')
    to_fasta_file(output,seqs, names)

    seqs, names = simulate_genome_by_composition(p_random, 10**4, 1000, 1)
    output = os.path.join(outdir, 'random.fasta')
    to_fasta_file(output,seqs, names)

    
def run_stretchFinder(type, index, m=10**4, w=100, joint=False):
    """
    run stretch finder on each sequence in a family data frame
    :param family: the dataset type name
    :param index: an index to be used by the PBS system
    :param w: the size of the linear window stretch (default =100)
    :return: saves the data frame into a file
    """

    df = pd.read_csv(
        r'/sternadi/home/volume1/daniellem1/Entropy/Simulations/Classifications/{}_profile.csv'.format(type))
    if joint:
        df = pd.read_csv(
            r'/sternadi/home/volume1/daniellem1/Entropy/Simulations/Classifications/{}_Joint_profile.csv'.format(type))

    col = 'seq_{}'.format(args.index - 1)
    seq = np.array(df[col].dropna())
    res = stretchFinder(seq, w, m)
    out_dir = r'/sternadi/home/volume1/daniellem1/EntropySimulations/Classifications/{}'.format(type)

    # create a directory per family if do not exist, if so, don't touch
    if not os.path.exists(out_dir):
        os.makedirs(out_dir)

    out = r'/sternadi/home/volume1/daniellem1/Entropy/Simulations/Classifications/{}/{}_{}_stats.csv'.format(type, type, col)
    if joint:
        out = r'/sternadi/home/volume1/daniellem1/Entropy/Simulations/Classifications/{}/{}_joint_{}_stats.csv'.format(type,
                                                                                                                 type,
                                                                                                                 col)
    res.to_csv(out, index=False)



def merge_stats_by_type(type, joint=False):
    """
    the method merges all stats csv files into a single file
    :param family: the family name
    :param joint: indicator for joint entropy. if true run the analysis on joint entropy
    :return: saves a csv with all p values
    """

    # define the folder from which all stats are going to be generated
    family_folder = r'/sternadi/home/volume1/daniellem1/Entropy/Simulations/Classifications/{}'.format(type)

    stats = glob.glob(os.path.join(family_folder, '*seq*stats.csv'))
    mapping = {}
    for f in stats:
        print(f)
        df = pd.read_csv(f)
        seq_id = re.findall(r'seq_\d+', f)[0]
        mapping[seq_id] = df['corrected_pvalue']

    # create a complete dataframe with all p values
    stats_by_seq = pd.DataFrame(dict([(k, pd.Series(v)) for k, v in mapping.items()]))

    output = os.path.join(family_folder, 'all_stats.csv')

    stats_by_seq.to_csv(output, index=False)
    print('Done!')
    return stats_by_seq


def analyse_simulations(data_dir, type):
    """
    filter the entropy and joint entropy measures to include only significant drops
    :param data_dir:
    :return:
    """





def main(args):
    output_dir = r'/sternadi/home/volume1/daniellem1/Entropy/Simulations/Drops'
    #files = glob.glob(output_dir)
    #fasta = files[args.index - 1]



    fasta = r'/volumes/STERNADILABHOME$/volume1/daniellem1/Entropy/Simulations/Drops/simulated.fasta'
    if args.mode == 1:
        ent = get_entropy_profile(fasta, 200, output_dir)
    else:
        joint = get_joint_entropy_profile(fasta, 200, output_dir)



if __name__ == "__main__":
    parser = argparse.ArgumentParser()
    #parser.add_argument("-i", "--index", type=int, help="array index")
    parser.add_argument("-m", "--mode", type=int, help="mode of running - entropy or join")

    args = parser.parse_args()

    main(args)