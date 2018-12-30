import pandas as pd
import numpy as np
import os
from tqdm import tqdm
import re
from Bio import SeqIO, Phylo
import math
from Bio.Seq import Seq
import random
from ete3 import Tree, TextFace, TreeStyle
from gradients import *



def refseq_2_seq(refseqs, fasta):
    """
    takes a list of refseq id's and returns the corresponding sequences
    :param refseqs: a list of refseq id's
    :param fasta: a path to a fasta file containing all genomes
    :return: a dictionary of refseq id as a key and its corresponding sequence
    """

    id_2_seq = {}

    sequences = re.split(">", open(fasta, "r").read().replace('\n', ''))[1:]
    for seq in tqdm(sequences):
        if '.' not in seq:
            print('no dot in sequence name\n')
            continue
        if 'complete genome' not in seq:
            print('not complete genome\n')
            continue

        # get identifier and genomic sequence
        splitted = seq.split('.')
        identifier = splitted[0].split()[0]
        genome = splitted[-1]

        if identifier in refseqs:
            if identifier in id_2_seq:
                print('Refseq id is not unique!')

            id_2_seq[identifier] = genome

    return id_2_seq

def extract_fasta_seqs(in_fasta, out_fasta, ids):
    """
    parses a fasta file to contain only the sequences in ids
    :param in_fasta: input fasta file
    :param out_fasta: output fasta path
    :param ids: a list of identifications
    :return: saves the new_fasta file
    """

    fin = open(in_fasta, 'r')
    fout = open(out_fasta, 'w')

    for record in SeqIO.parse(fin, 'fasta'):
        for item in ids:
            if item.strip() == record.id and 'complete genome' in record.description:
                fout.write(">" + record.id + "\n")
                fout.write(str(record.seq) + "\n")

    fin.close()
    fout.close()


def entropy_by_kmer(seq, k):
    """
    calculate the entropy of a sequence according to its k-mers
    :param seq: a genome string
    :param k: size of the k-mer to use
    :return: entropy
    """

    # update kmers
    kmers = {}
    for i in range(len(seq) - k):
        kmer = seq[i:i+k]
        if kmer in kmers:
            kmers[kmer] += 1
        else:
            kmers[kmer] = 1

    # calculate entropy
    total_kmers = sum(kmers.values())
    entropy = 0
    for kmer in kmers:
        p = kmers[kmer] / total_kmers
        entropy += -(p * math.log2(p))

    return entropy

def joint_entropy (seq1, seq2, k):
    """
    calculates the joint entropy of two sequences.
    :param seq1: sequence #1
    :param seq2: sequence #2
    :param k: k-mer length
    :return: joint entropy value
    """

    kmers_1 = {}
    kmers_2 = {}

    # kmers in sequence #1
    for i in range(len(seq1) - k):
        kmer = seq1[i:i+k]
        if kmer in kmers_1:
            kmers_1[kmer] += 1
        else:
            kmers_1[kmer] = 1

    for i in range(len(seq2) - k):
        kmer = seq2[i:i+k]
        if kmer in kmers_2:
            kmers_2[kmer] += 1
        else:
            kmers_2[kmer] = 1

    # calculate joint entropy
    total_kmers_1 = sum(kmers_1.values())
    total_kmers_2 = sum(kmers_2.values())

    total = total_kmers_1 + total_kmers_2

    # compare the kmers space to be equal at both
    for kmer in kmers_1:
        if kmer not in kmers_2:
            kmers_2[kmer] = 0

    for kmer in kmers_2:
        if kmer not in kmers_1:
            kmers_2[kmer] = 0

    joint_entropy = 0
    for kmer1 in kmers_1:
        for kmer2 in kmers_2:
            p_xy = (kmers_1[kmer1] + kmers_2[kmer2]) / total

            joint_entropy += -(p_xy * math.log2(p_xy))

    return joint_entropy


def information_storage (seq1, seq2, k):
    """
    calculates the information storage of two sequences.
    :param seq1: sequence #1
    :param seq2: sequence #2
    :param k: k-mer length
    :return: information storage value
    """

    kmers_1 = {}
    kmers_2 = {}

    # kmers in sequence #1
    for i in range(len(seq1) - k):
        kmer = seq1[i:i+k]
        if kmer in kmers_1:
            kmers_1[kmer] += 1
        else:
            kmers_1[kmer] = 1

    for i in range(len(seq2) - k):
        kmer = seq2[i:i+k]
        if kmer in kmers_2:
            kmers_2[kmer] += 1
        else:
            kmers_2[kmer] = 1

    # calculate joint entropy
    total_kmers_1 = sum(kmers_1.values())
    total_kmers_2 = sum(kmers_2.values())

    total = total_kmers_1 + total_kmers_2

    # compare the kmers space to be equal at both
    for kmer in kmers_1:
        if kmer not in kmers_2:
            kmers_2[kmer] = 0

    for kmer in kmers_2:
        if kmer not in kmers_1:
            kmers_2[kmer] = 0

    inf_storage = 0
    for kmer1 in kmers_1:
        for kmer2 in kmers_2:

            p_xy = (kmers_1[kmer1] + kmers_2[kmer2]) / total
            p_x = kmers_1[kmer1] / total_kmers_1
            p_y = kmers_2[kmer2] / total_kmers_2

            if p_x == 0 or p_y == 0:
                continue
            inf_storage += p_xy * math.log2(p_xy/(p_x*p_y))

    return inf_storage


def reading_frame_entropy(seq):
    """
    get entropy by reading frame
    :param seq: a genome sequence
    :param k: kmer length
    :return: entropy value by reading frame
    """

    kmers = {}
    for i in range(0,len(seq) - 3,3):
        kmer = seq[i:i+3]
        if kmer in kmers:
            kmers[kmer] += 1
        else:
            kmers[kmer] = 1

    # calculate entropy
    total_kmers = sum(kmers.values())
    entropy = 0
    for kmer in kmers:
        p = kmers[kmer] / total_kmers
        entropy += -(p * math.log2(p))

    return entropy



def get_reverse_complement(seq):
    """
    get reverse complement genome
    :param seq: a genome sequence
    :return: a string of reverse complement
    """
    seq = Seq(seq)
    reverse_complement = seq.reverse_complement()
    return reverse_complement


def _scrambler(word):
    word_to_scramble = list(word)
    random.shuffle(word_to_scramble)
    new_word = ''.join(word_to_scramble)
    return new_word

def scrambler(word):
   new_word = _scrambler(word)
   while new_word == word and len(word) > 1:
       new_word = _scrambler(word)
   return new_word



def entropy_2_refseq(genomic, cds, out, k):
    """
    calculate for each complete genome its entropy, joint entropy with reverse complement and information storage
    :param db_path: a path to a db with viruses information
    :param genomic: a genomic fasta file
    :param cds: coding regions fasta file
    :param out: output path to save the files
    :param k: kmer length
    :return: a data frame with entropy and virus information
    """

    refseq_g = []
    entropy_g = []
    joint_g = []
    information_storage_g = []
    size_g = []
    scrambled_joint_g = []

    refseq_p = []
    frame_entropy_p = []
    entropy_p = []
    protein_p = []
    joint_p = []
    information_storage_p = []
    size_p = []
    scrambled_joint_p = []



    genomic_sequences = re.split(">", open(genomic, "r").read().replace('\n', ''))[1:]
    for seq in tqdm(genomic_sequences):
        if '.' not in seq:
            print('no dot in sequence name\n')
            continue
        if 'complete genome' not in seq:
            #print('not complete genome\n')
            continue

        # get identifier and genomic sequence
        splitted = seq.split('.')
        identifier = splitted[0].split()[0]
        genome = splitted[-1]

        rc_genome = get_reverse_complement(genome)
        rc_genome_scrambled = scrambler(rc_genome)

        joint = joint_entropy(genome, rc_genome, k)
        info_storage = information_storage(genome, rc_genome, k)
        entropy = entropy_by_kmer(genome, k)
        joint_scrambled = joint_entropy(genome, rc_genome_scrambled, k)

        refseq_g.append(identifier)
        entropy_g.append(entropy)
        joint_g.append(joint)
        information_storage_g.append(info_storage)
        size_g.append(len(genome))
        scrambled_joint_g.append(joint_scrambled)

    cds_sequences = re.split(">", open(cds, "r").read().replace('\n', ''))[1:]
    for seq in tqdm(cds_sequences):

        # get identifier and genomic sequence
        splitted = seq.split('|')
        identifier = splitted[0]
        assert(splitted[2] == '')
        refseq_id = splitted[3]
        genome = splitted[-1]

        rc_genome = get_reverse_complement(genome)
        rc_genome_scrambled = scrambler(rc_genome)

        joint = joint_entropy(genome, rc_genome, k)
        info_storage = information_storage(genome, rc_genome, k)
        entropy = entropy_by_kmer(genome, k)
        rf_entropy = reading_frame_entropy(genome)
        joint_scrambled = joint_entropy(genome, rc_genome_scrambled, k)

        refseq_p.append(refseq_id)
        entropy_p.append(entropy)
        joint_p.append(joint)
        information_storage_p.append(info_storage)
        protein_p.append(identifier)
        frame_entropy_p.append(rf_entropy)
        size_p.append(len(genome))
        scrambled_joint_p.append(joint_scrambled)

    # create a data frame for each cds\ genomic data and merge
    genomic_df = pd.DataFrame({'refseq_id': refseq_g, 'entropy_genome':entropy_g, 'genome_size':size_g,
                               'joint_entropy_genome':joint_g, 'information_storage_genome':information_storage_g,
                               'scrambled_joint_genome': scrambled_joint_g})

    genomic_df.to_csv(os.path.join(out, 'entropy_measurements_genomics_{}.csv'.format(k)), index=False)

    cds_df = pd.DataFrame({'refseq_id': refseq_p, 'protein_id':protein_p, 'entropy_cds': entropy_p, 'protein_size': size_p,
                               'joint_entropy_cds': joint_p, 'information_storage_cds': information_storage_p,
                           'rf_entropy':frame_entropy_p,'scrambled_joint_cds':scrambled_joint_p})

    cds_df.to_csv(os.path.join(out, 'entropy_measurements_proteins_{}.csv'.format(k)), index=False)

    combined = pd.concat([genomic_df, cds_df], axis=0, ignore_index='True')

    combined.to_csv(os.path.join(out, 'entropy_measurements_{}.csv'.format(k)), index=False)

    return combined


def split_fasta_by_feature(input_fasta, out_dir, db_path, feature):
    """
    splits a fasta file to sub-fasta files according some feature. creates a directory for each feature value and store
    the fasta file in it
    :param input_fasta: input fasta file of complete genomes
    :param out_dir: outputdit to create the feature values directories in
    :param db_path: a path to a db containing features and refseq ids
    :return: saves fasta files in a directory
    """

    virus_host_db = pd.read_excel(db_path)

    for elem in tqdm(set(virus_host_db[feature].dropna())):
        refseq_by_feature = virus_host_db[virus_host_db[feature] == elem]['refseq_id'].values
        alias = elem.strip()
        directory = os.path.join(out_dir, alias)
        if not os.path.exists(directory):
            os.makedirs(directory)

        out_fasta = os.path.join(directory, '{}.fasta'.format(alias))
        extract_fasta_seqs(input_fasta, out_fasta, refseq_by_feature)


def phylo_statistic_by_feature(input_dir, out):
    """
    create a statistics table with sequences information, number of sequences, sd of sequences length,
    mean sequence length,
    :param input_dir: the directory in which the folder containing fasta files are in
    :param out: an output directory path to save the result
    :return: a data frame with the information about every family
    """

    feature = os.path.basename(input_dir)

    # get all fasta files. each fasta file name contains the relevant class that will correspond to a line in the table
    all_fasta = []

    separation_class = []
    all_sd = []
    all_mean= []
    all_median = []
    all_min = []
    all_max = []
    all_num_seqs = []


    for root, dirs, files in os.walk(input_dir):
        fasta_file = [f for f in files if '.fasta' in f]
        if fasta_file != []:
            all_fasta.append(os.path.join(root, fasta_file[0]))

    for f in tqdm(all_fasta):
        sequences = re.split(">", open(f, "r").read().replace('\n', ''))[1:]
        if sequences == []:
            print('no sequences for type {}'.format(os.path.basename(f).split('.')[0]))
            continue
        splitted = [re.split(r'(\d+)', s) for s in sequences]
        only_sequence = [elem[-1] for elem in splitted]

        len_sequences = [len(seq) for seq in only_sequence]

        num_seqs = len(only_sequence)
        min_len = np.min(len_sequences)
        max_len = np.max(len_sequences)
        mean_len = np.mean(len_sequences)
        median_len = np.median(len_sequences)
        std = np.std(len_sequences)
        name = os.path.basename(f).split('.')[0]

        # add the values to the correspond list
        separation_class.append(name)
        all_sd.append(std)
        all_mean.append(mean_len)
        all_median.append(median_len)
        all_min.append(min_len)
        all_max.append(max_len)
        all_num_seqs.append(num_seqs)

    # create a data frame and save it to the output directory
    result = pd.DataFrame({'class':separation_class, 'num_seqs':all_num_seqs, 'mean':all_mean, 'median':all_median,
                           'min':all_min, 'max':all_max, 'std':all_sd})
    result.to_csv(os.path.join(out,'{}_sequences_information.csv'.format(feature)), index=False)

    return result



def value_2_color(values, all_colors):
    """

    :param values:
    :param all_colors:
    :return:
    """
    result = {}

    for i, sorted_idx in enumerate(np.argsort(values)):
       result[str(values[sorted_idx])] = all_colors[i]

    return result

def tree_2_string(tree_file):
    with open(tree_file, "r") as handle:
        return handle.read()

def add_labels_n_colors_2_tree(tree_path, refseq_2_entropy, colors, feature, out):
    """
    add entropy measurement to a tree and color it's leaves
    :param tree_path: a path to a tree
    :param refseq_2_entropy: a data base with entropy values for each refseq id
    :param out: output directory path to save the results
    :return: a tree with colors and entropy labels
    """
    print(tree_path)
    str_tree = tree_2_string(tree_path)
    if str_tree == '':
        return
    tree = Tree(str_tree)


    # get all refseq ids in the tree and all corresponding entropy values
    all_leaves = []

    for node in tree:
        if node.is_leaf():
            all_leaves.append(node.name)

    #all_values = refseq_2_entropy[refseq_2_entropy['refseq_id'].isin(all_leaves)]['entropy_5'].values
    all_values = refseq_2_entropy[feature].values

    num_leaves = len(all_leaves)

    # create a gradient color bar
    #colors = linear_gradient("#DC9221", "#33C193", n=num_leaves)

    #mapping = value_2_color(all_values, colors['hex'])
    mapping = value_2_color(all_values, colors['hex'])

    for node in tree:
        if node.is_leaf():
            value = refseq_2_entropy[refseq_2_entropy['refseq_id'] == node.name.split('.')[0]][feature].values[0]
            virus_name = refseq_2_entropy[refseq_2_entropy['refseq_id'] == node.name.split('.')[0]]['virus_name'].values[0]

            # change virus name
            node.name = virus_name
            c = mapping[str(value)]
            node.add_feature("label", value)
            label_face = TextFace(node.label, fgcolor=c, fsize=22)
            node.add_face(label_face, column=0, position='branch-right')


    family = os.path.basename(tree_path).split('.')[0]

    ts = TreeStyle()
    ts.mode = 'c'
    ts.scale = 5
    # Disable the default tip names config
    ts.show_leaf_name = True
    ts.show_scale = True
    ts.show_branch_length = True
    ts.title.add_face(TextFace("{} values for {}".format(feature, family), fsize=20), column=0)
    # Draw Tree
    tree.render(os.path.join(out, '{}_{}_tree.pdf'.format(family, feature)), dpi=300, tree_style=ts)



def phylo_entropy_construction(db, df ,out):
    """
    create a data frame with extreme entropy changed in each tree, containing the change in sequences length and the sum
    of branchs between them both extreme values
    :param db: a path to a folders of viruses folders
    :param df: mapping of tree leaf name (id) to the requested numerical value (value) and to sequence length
    :param out: a pth to save the results
    :return: the described data frame
    """
    all_df = []

    all_trees = []
    for root, dirs, files in tqdm(os.walk(db)):
        tree = [f for f in files if 'phyml_tree' in f]
        if tree != []:
            all_trees.append(os.path.join(root, tree[0]))

    for f in tqdm(all_trees):
        if tree_2_string(f) == '':
            continue
        tree = Phylo.read(f, 'newick')
        term_names = [term.name.split('.')[0] for term in tree.get_terminals()]
        entropy_values = [(name, df[df['refseq_id'] == name]['entropy_5'].values[0]) for name in term_names]

        max_tup = max(entropy_values, key=lambda x: x[1])
        min_tup = min(entropy_values, key=lambda x: x[1])

        len_max_seq =  df[df['refseq_id'] == max_tup[0]]['genome_size_x'].values[0]
        len_min_seq = df[df['refseq_id'] == min_tup[0]]['genome_size_x'].values[0]

        branch_distance = tree.distance(max_tup[0], min_tup[0])

        tree_df = pd.DataFrame({'max_refseq_id': max_tup[0], 'min_refseq_id': min_tup[0], 'max_entropy': max_tup[1],
                                'min_entropy': min_tup[1], 'max_genome_size': len_max_seq, 'min_genome_size': len_min_seq,
                                'sum_branch_length': branch_distance,'family': os.path.basename(f).split('.')[0],
                               'num_sequences_in_tree': len(term_names)}, index=[0])
        all_df.append(tree_df)


    result = pd.concat(all_df)
    result['delta_len'] = result.apply(lambda row: math.fabs(row['max_genome_size'] - row['min_genome_size']), axis=1)
    result['delta_entropy'] = result.apply(lambda row: math.fabs(row['max_entropy'] - row['min_entropy']), axis =1)

    result.to_csv(out, index=False)

    return result


def refseq_2_content(genomic, virus_host, out):
    ''' crate a data frame of ids and their correspond nuc content in genome'''

    result = []

    genomic_sequences = re.split(">", open(genomic, "r").read().replace('\n', ''))[1:]
    for seq in tqdm(genomic_sequences):
        if '.' not in seq:
            print('no dot in sequence name\n')
            continue
        if 'complete genome' not in seq:
            continue

        # get identifier and genomic sequence
        splitted = seq.split('.')
        identifier = splitted[0].split()[0]
        genome = splitted[-1]
        n = len(genome)
        a_content = (genome.count('a') / n) * 100
        c_content = (genome.count('c') / n) * 100
        g_content = (genome.count('g') / n) * 100
        t_content = (genome.count('t') / n) * 100

        df = pd.DataFrame({'refseq_id':identifier, 'a_content':a_content, 'c_content':c_content,
                           'g_content':g_content,'t_content':t_content}, index=[0])
        result.append(df)

    concat = pd.concat(result)
    merged = pd.merge(virus_host, concat, on='refseq_id')
    merged.to_csv(os.path.join(out, "refseq_2_content.csv"), index=False)

def phylo_content_construction(db, df, feature, out):
    """
    create a data frame with extreme entropy changed in each tree, containing the change in sequences length and the sum
    of branchs between them both extreme values
    :param db: a path to a folders of viruses folders
    :param df: mapping of tree leaf name (id) to the requested numerical value (value) and to sequence length
    :param out: a pth to save the results
    :return: the described data frame
    """
    all_df = []

    all_trees = []
    for root, dirs, files in tqdm(os.walk(db)):
        tree = [f for f in files if 'phyml_tree' in f]
        if tree != []:
            all_trees.append(os.path.join(root, tree[0]))

    for f in tqdm(all_trees):
        if tree_2_string(f) == '':
            continue
        tree = Phylo.read(f, 'newick')
        term_names = [term.name.split('.')[0] for term in tree.get_terminals()]
        entropy_values = [(name, df[df['refseq_id'] == name][feature].values[0]) for name in term_names]

        max_tup = max(entropy_values, key=lambda x: x[1])
        min_tup = min(entropy_values, key=lambda x: x[1])

        len_max_seq =  df[df['refseq_id'] == max_tup[0]]['genome_size_x'].values[0]
        len_min_seq = df[df['refseq_id'] == min_tup[0]]['genome_size_x'].values[0]

        branch_distance = tree.distance(max_tup[0], min_tup[0])

        tree_df = pd.DataFrame({'max_refseq_id': max_tup[0], 'min_refseq_id': min_tup[0], 'max_{}'.format(feature): max_tup[1],
                                'min_{}'.format(feature): min_tup[1], 'max_genome_size': len_max_seq, 'min_genome_size': len_min_seq,
                                'sum_branch_length': branch_distance,'family': os.path.basename(f).split('.')[0],
                               'num_sequences_in_tree': len(term_names)}, index=[0])
        all_df.append(tree_df)


    result = pd.concat(all_df)
    result['delta_len'] = result.apply(lambda row: math.fabs(row['max_genome_size'] - row['min_genome_size']), axis=1)
    result['delta_{}'.format(feature)] = result.apply(lambda row: math.fabs(row['max_{}'.format(feature)] - row['min_{}'.format(feature)]), axis =1)

    result.to_csv(out, index=False)

    return result


def remove_duplicated_from_aln(aln, alias):
    """
    removes duplicate sequences for a given alignment
    :param aln: a path to an alignment file
    :return: alignment without duplications
    """

    record_dict = SeqIO.to_dict(SeqIO.parse(aln, "fasta"))
    no_duplicates = {key: val for key, val in record_dict.items() if '.1' not in key}
    output = os.path.join(os.path.dirname(aln), alias + '.filtered.aln.fas')
    with open(output, 'w') as o:
        SeqIO.write(no_duplicates.values(), o, "fasta")




def prune_tree(aln, tree_path, alias, threshold=0.5):
    """
    removes the tips in which their branch length exceeds some threshold
    :param aln: a path to an aln file
    :param tree: a path to an existing tree files
    :return: saves a new aln file named filtered pruned, to the prune directory
    """
    if tree_2_string(tree_path) == '':
        return

    tree = Phylo.read(tree_path, 'newick')

    # get all terminal names to exlude from the aln
    need_2_remove = []
    for clade in list(tree.find_clades()):
        if clade.branch_length == None:
            continue
        if clade.branch_length > threshold:
            need_2_remove.extend([l.name for l in list(clade.get_terminals())])

    # create a new aln exluding the problematic nodes.
    record_dict = SeqIO.to_dict(SeqIO.parse(aln, "fasta"))
    no_duplicates = {key: val for key, val in record_dict.items() if '.1' not in key and key not in need_2_remove}
    output = os.path.join('/Volumes/STERNADILABHOME$/volume1/daniellem1/Entropy/data/Phylogeny/pruned_trees/', alias +'.filtered_prune.aln.fas')

    with open(output, 'w') as o:
        SeqIO.write(no_duplicates.values(), o, "fasta")
