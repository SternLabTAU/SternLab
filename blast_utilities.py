#! /usr/local/python_anaconda/bin/python3.4

import pandas as pd
import re

mutation_pattern = re.compile('[ACGT]{2}')
insertion_pattern = re.compile('-[ACGT]')
deletion_pattern = re.compile('[ACGT]-')
number_pattern = re.compile('\d+')
pattern = re.compile('[ACGT-]{2}|\d+')



def blast_to_mutations_list(blast_path, out_csv_path):
    '''
    Gets a blast file and creates a file containing all single nucleotide mutations
    (point mutations, insertions, and deletions) with the matching read name.
    :param blast_path: input file path, needs to be blast output
    :param out_csv_path: output file path
    '''
    blast = pd.read_csv(blast_path, sep="\t", header=None)
    blast.columns = ["read", "start_ref", "end_ref", "start_read", "end_read", "strand", "length", "btop"]
    blast['parsed_btop'] = blast.apply(parse_btop, axis=1)
    reads = []
    positions = []
    refs = []
    bases = []
    for index, row in blast.iterrows():
        for mutations in list(row['parsed_btop'].values()):
            for i in mutations:
                reads.append(row['read'])
                positions.append(i[0])
                refs.append(i[1][0])
                bases.append(i[1][1])
    out = pd.DataFrame()
    out['read'] = reads
    out['position'] = positions
    out['ref'] = refs
    out['base'] = bases
    out.to_csv(out_csv_path, index=False)
    return

def parse_btop(row):
    '''
    Gets a pandas dataframe row from a blast file and parses the btop field.
    Used in get_mutations_from_blast function.
    '''
    insertions = []
    deletions = []
    mutations = []
    location = float(row['start_ref'])
    btops = pattern.findall(row['btop'])
    for b in btops:
        if number_pattern.findall(b) != []:
            location = float(int(location) + int(b))
        if mutation_pattern.findall(b) != []:
            mutations.append((location, b))
            location = float(int(location) + 1)
        if insertion_pattern.findall(b) != []:
            location += 0.01
            insertions.append((location, b))
        if deletion_pattern.findall(b) != []:
            deletions.append((location, b))
            location = float(int(location) + 1) 
    return {'mutations':mutations, 'deletions':deletions, 'insertions':insertions}