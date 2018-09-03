#! /usr/local/python_anaconda/bin/python3.4

import pandas as pd
import re
import math

def blast_to_df(blast_path):
    '''
    Creates pandas dataframe from blast file
    :param blast_path: input file path, needs to be blast output
    :returns: pandas dataframe
    '''
    blast_df = pd.read_csv(blast_path, sep='\t', header = None)
    blast_df.columns = ("read", "start_ref", "end_ref", "start_read", "end_read", "strand", "length", "btop")
    return blast_df

def blast_to_mutations_list(blast_path, out_csv_path=None):
    '''
    Gets a blast file and creates a file containing all single nucleotide mutations
    (point mutations, insertions, and deletions) with the matching read name.
    :param blast_path: input file path, needs to be blast output
    :param out_csv_path: output file path
    :returns: pandas dataframe
    '''
    blast = blast_to_df(blast_path)
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
    if out_csv_path:
        out.to_csv(out_csv_path, index=False)
    return out

def parse_btop(row):
    '''
    Gets a pandas dataframe row from a blast file and parses the btop field.
    Used in blast_to_mutations_list function.
    '''
    mutation_pattern = re.compile('[ACGTN]{2}')
    insertion_pattern = re.compile('-[ACGTN]')
    deletion_pattern = re.compile('[ACGTN]-')
    number_pattern = re.compile('\d+')
    pattern = re.compile('[ACGTN-]{2}|\d+')

    insertions = []
    deletions = []
    mutations = []
    location = float(row['start_ref'])
    btops = pattern.findall(row['btop'])
    for b in btops:
        if number_pattern.findall(b) != []:
            location = float(math.floor(location) + int(b))
        if mutation_pattern.findall(b) != []:
            mutations.append((float(math.ceil(location)), b))
            location = float(math.ceil(location) + 1)
        if insertion_pattern.findall(b) != []:
            location += 0.01
            insertions.append((location, b))
        if deletion_pattern.findall(b) != []:
            deletions.append((float(math.ceil(location)), b))
            location = float(math.ceil(location) + 1) 
    return {'mutations':mutations, 'deletions':deletions, 'insertions':insertions}