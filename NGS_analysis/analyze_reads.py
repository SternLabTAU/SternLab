import pandas as pd
import os
from itertools import combinations
from tqdm import tqdm
import argparse



def main(args):

    blast_files = get_blast_files(args.in_dir)
    merged_blast = merge_blast_files(blast_files)

    mutations = ['A1664G', 'A1744G']
    read_mapping = get_reads_mutations(merged_blast, mutations)
    results = create_df_from_dict(read_mapping, mutations)

    base = os.path.basename(args.in_dir)
    results.to_csv(os.path.join(args.out, "{]_{}".format(base, 'reads_stats.csv')), index=False)



def get_blast_files(running_directory):
    """
    returns all blast files in pipeline tmp directory
    :param running_directory: a path to an output pipeline directory
    :return: a list with all blast paths
    """

    blasts = [f for f in os.listdir(os.path.join(running_directory, 'tmp')) if f.endswith('.blast')]
    return blasts



def merge_blast_files(blast_files):
    """
    merges all part blast files into a single data frame
    :param blast_files: a list of blast files path
    :return: a data frame with all blast files.
    """

    all_parts = []

    for f in blast_files:
        col_names = ['read','start', 'end', 'x','y','strand','size','read_mutation']
        blast_part = pd.read_csv(f, sep='\t', names = col_names)
        all_parts.append(blast_part)

    unite_blast = pd.concat(all_parts)
    return unite_blast

def get_reads_mutations(reads_data, mutations):
    """
    for each read check how many mutations are on it and identify them.
    :param reads_data: a data frame containing the read mutation positions on read and their start indices
    :param mutations: a vector of mutations in a form : ref pos base
    :return: a dictionary with count of mutation that come together for each possible mutation combination
    """

    # create the dictionary - keys are all mutations combinations
    mut_combinations = []
    for i in range(len(mutations)):
        mut_combinations.extend(list(combinations(mutations, i)))

    mapping = dict()
    for comb in mut_combinations:
        key = [x for x in comb]
        mapping[key] = 0

    # parse each read to get the mutations on that read.
    for i in tqdm(reads_data.shape[0]):
        start = reads_data.loc[i]['start']
        read = reads_data.loc[i]['read_mutation']

        # get the index of each mutations on the read. a mutation can start with A T C G
        mut_indices = [x for x, ltr in enumerate(read) if ltr == 'A' or ltr == 'C' or ltr == 'T' or ltr == 'G' or ltr == '-']
        mutated_positions = []

        # if the read has no mutations, break and assign to the dictionary +1 in the empty list key
        if mut_indices == []:
            mapping[mutated_positions] += 1
            break

        # iterate over the read to get the positions of mutations on it
        offset = 0
        for j in tqdm(range(0, len(mut_indices), 2)):
            idx = mut_indices[j]
            if i == 0:
                if read[:idx] == "":
                    offset +=1
                else:
                    offset += (int(read[:idx]) + 1)
            else:
                if read[mut_indices[j-1] +1: idx] == "":
                    offset += 1
                else:
                    offset += (int(read[mut_indices[j-1] +1: idx]) + 1)

            position = start + offset
            mut = "{}{}{}".format(read[mut_indices[j]], position, read[mut_indices[j + 1]])
            mutated_positions.append(mut)

        # update dictionary
        mapping[mutated_positions] += 1

    return mapping



def create_df_from_dict(mutations_mapping_dict, mutations):
    """
    create a data frame with with columns as indicators for each mutation True = the mutation is on the read False ow.
    :param mutations_mapping_dict: a dict with lists of mutation combinations and the count of the reads for each combinations
    :param mutations: the list of mutations received by the user
    :return: a data frame with mutations as columns, number of reads in which a combination share and total number of reads
    """

    mapping = dict()
    keys = mutations_mapping_dict.keys()
    vals = mutations_mapping_dict.values()

    # iterate over mutations and get a vector of true false for appearances in keys list
    for mut in mutations:
        is_in_idx = [keys.index(x) for x in keys if mut in x]
        appear = []
        for i in range(len(keys)):
            if i in is_in_idx:
                appear.append(True)
            else:
                appear.append(False)
        mapping[mut] = appear

    df = pd.DataFrame(mapping, columns=keys)
    df['Number_of_reads'] = vals
    df['Total'] = sum(vals)

    return df



if __name__ == "__main__":
    parser = argparse.ArgumentParser()
    parser.add_argument("-i", "--in_dir", type=str, help="input dir with pipeline running results. should include tmp file", required=True)
    parser.add_argument("-o", "--out", type=str, help="output directory to save the files", required=True)
    args = parser.parse_args()
    main(args)













