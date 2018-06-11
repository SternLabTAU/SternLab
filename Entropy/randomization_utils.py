import matplotlib.pyplot as plt
import seaborn as sns
from utils import *
from itertools import combinations
from scipy import stats
import random
from Bio.SeqUtils import CodonUsage
from entropy_selection import *
import pickle


def refseq_2_cds_positions(cds, out=None):
    """
    creates a mapping of refseq id to cds positions
    :param cds: cds fasta file
    :return: a dictionary containing the mapping
    """


    # create a mapping between each id all the coding sequences it has.
    mapping = {}
    sequences = re.split(">", open(cds, "r").read().replace('\n', ''))[1:]
    for seq in tqdm(sequences):
        splitted = seq.split('|')
        # if splitted[2] != '':
            # print(splitted[0], splitted[2], splitted[3])
        refseq_id = remove_punctuation(splitted[3])
        if 'geneid' in refseq_id.lower():
            refseq_id = remove_punctuation(splitted[4])
        # refseq_id = splitted[3]
        if '}' in seq:
            str_positions = splitted[-1].split('}')[0].split('|')[-1]
            lst_positions = re.findall(r'\d+', str_positions)
        else:
            str_positions = splitted[-1].split(')')[0].split('|')[-1]
            lst_positions = re.findall(r'\d+', str_positions)

        # add the sequence to the dictionary
        if refseq_id in mapping:
            mapping[refseq_id].extend(lst_positions)
        else:
            mapping[refseq_id] = lst_positions

    for key in mapping.keys():
        positions = mapping[key]
        mapping[key] = sorted([(positions[i], positions[i+1]) for i in range(0, len(positions),2)], key=lambda x: x[0])

    with open (out, 'wb') as o:
        pickle.dump(mapping, o)
    return mapping



def get_cds_by_position(seq, positions):
    """
    return a list of coding regions sequences by positions
    :param seq: a whole genome sequence
    :param positions: a list of tuples (start position, end position)
    :return: a list of coding regions sequences
    """

    return [seq[t1:t2] for t1,t2 in positions]


def scramble_codons_arbitrarily(seq):
    """

    :param seq: a cds sequence
    :param positions: a list of tuple (start positions, end position)
    :return: a shuffled sequence by codons (changes the protein completly)
    """
    all_codons = [seq[i:i+3] for i in range(0, len(seq), 3)]
    random.shuffle(all_codons)
    return ''.join(all_codons)

def replace_codon(codon, syn_dict, codon_dict):
    """
    replaces the codon given by a synonymous codon that is still available
    :param codon: a given codon
    :param syn_dict: synonymous codons dict
    :param codon_dict: codon distribution of a sequence
    :return: new synonymous codon available
    """
    aa_type = [k for k,v in syn_dict.items() if codon in v][0]
    syn_codons = list(syn_dict[aa_type])
    valid_options = [c for c in syn_codons if codon_dict[c] != 0 and c != codon]
    if valid_options == []:
        return codon
    else:
        random.shuffle(valid_options)
        return valid_options[0]


def scramble_synonymous_codons(seq):
    """
    scramble synonymous codons to keep the protein the same, with the same nucleotide compositions
    :param seq: a cds sequence
    :return: a shuffled sequence
    """

    seq = seq.upper()
    codon_dict = CodonUsage.CodonsDict
    syn_codons = CodonUsage.SynonymousCodons
    all_codons = [seq[i:i + 3] for i in range(0, len(seq), 3) if len(seq[i:i + 3]) == 3]

    for codon in all_codons:
        codon_dict[codon] += 1

    shuffled = ''
    for i in range(0, len(seq), 3):
        codon = seq[i:i+3]
        print(codon)
        print(syn_codons[codon])
        new_codon = replace_codon(codon, syn_codons, codon_dict)

        # remove that option from the available codons
        codon_dict[new_codon] -= 1
        shuffled += new_codon

    return shuffled.lower()

def scramble_all_sequence(seq, mapping, refseq_id, how=1):
    """
    scramble whole genome by how
    :param seq: a whole genome sequence
    :param mapping: a mapping between the genome and its coding regions positions
    :param refseq_id: refseq id
    :param how: how to scramble - 1: random genome scrambling (default) 2: random codons scramble 3: synonymous codons scramble
    :return: scrambled genome
    """

    if refseq_id not in mapping:
        print('No id in mapping')
        return

    if how == 1:
        return scrambler(seq)

    positions = mapping[refseq_id]
    print(positions)
    all_positions =[]
    for t1, t2 in positions:
        all_positions.append(t1)
        all_positions.append(t2)
    print(all_positions)
    scrambled = scrambler(seq[:all_positions[0]])

    for i in range(0, len(all_positions)):
        # if we are in coding region
        if i % 2 == 0:
            coding = seq[all_positions[i]: all_positions[i+1] + 1]
            if how == 2:
                coding = scramble_codons_arbitrarily(coding)
            else:
                coding = scramble_synonymous_codons(coding)

            scrambled += coding
        else:   # not coding
            if i == len(all_positions) - 1: # end of the list
                non_coding = scrambler(seq[all_positions[i]+1:])
            else:
                non_coding = scrambler(seq[all_positions[i]+1: all_positions[i+1]])

            scrambled += non_coding

    return scrambled


