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

def simulate_genome(p, n, size, mode):
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
        for i in tqdm(range(n)):
            seq = ''.join(np.random.choice(['a', 'c', 'g', 't'], p=p, size=size // 2))
            seq = seq + 'aaaaaa' + str(get_reverse_complement(seq))
            sequences.append(seq)
            names.append('mode_{}_seq_{}'.format(mode, i))

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

    seqs, names = simulate_genome(p_rep, 10**4, 1000, 1)
    output = os.path.join(outdir, 'repetative.fasta')
    to_fasta_file(output,seqs, names)

    seqs, names = simulate_genome(p_rep, 10**4, 1000, 2)
    output = os.path.join(outdir, 'repetative_structure.fasta')
    to_fasta_file(output,seqs, names)

    seqs, names = simulate_genome(p_random, 10**4, 1000, 2)
    output = os.path.join(outdir, 'structure.fasta')
    to_fasta_file(output,seqs, names)

    seqs, names = simulate_genome(p_random, 10**4, 1000, 1)
    output = os.path.join(outdir, 'random.fasta')
    to_fasta_file(output,seqs, names)

def main(args):
    # calculate entropy values and profiles for each file
    output_dir = r'/sternadi/home/volume1/daniellem1/Entropy/Simulations/Classifications'
    files = glob.glob(output_dir)
    fasta = files[args.index - 1]

    if args.mode == 1:
        ent = get_entropy_profile(fasta, 200, output_dir)
    else:
        joint = get_joint_entropy_profile(fasta, 200, output_dir)



if __name__ == "__main__":
    parser = argparse.ArgumentParser()
    parser.add_argument("-i", "--index", type=int, help="array index")
    parser.add_argument("-m", "--mode", type=int, help="mode of running - entropy or join")

    args = parser.parse_args()

    main(args)