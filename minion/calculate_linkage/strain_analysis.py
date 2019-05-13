#! /usr/local/python_anaconda/bin/python3.4


import pandas as pd
import argparse

def count_haplotypes(args):
    '''
    This function separates the reads into strains. A strain is defined as any combination 
    of the mutations or lack of mutations that are provided to the function.
    
    Format for the csv with the mutations to check strains for: every mutation should 
    have its own row, the csv should have no header row, and the mutations should be written 
    in the following format: "A1664.0G". For example:
    "A1664.0G
    A535.0G
    T1764.0-"
    '''
    recognized_mutations = pd.read_csv(args.input_chosen_mutation, header=None)[0].tolist()
    recognized_positions = [float(p[1:-1]) for p in recognized_mutations]
    blast_df = pd.read_csv(args.input_blast_df)
    # choose only reads that were mapped only once in blast
    blast_df['read_count'] = blast_df.groupby('read')['start_ref'].transform('count')
    blast_df = blast_df[(blast_df.read_count == 1)]
    
    # choose only reads that are mapped from at least start_pos_read to end_pos_read
    blast_df = blast_df[(blast_df.start_ref < min(recognized_positions)) & (blast_df.end_ref > max(recognized_positions))]
    blast_df = blast_df[['read', 'start_ref', 'end_ref']]
    
    mutations_df = pd.read_csv(args.input_mutation_df)
    mutations_df = mutations_df[(mutations_df.ref != '-')]
    
    # drop reads containing a variation that is not the recognized mutation or the WT in the positions we are checking combinations for.
    reads_to_drop = mutations_df[(mutations_df.position.isin(recognized_positions)) & ~(mutations_df.full_mutation.isin(recognized_mutations))].read.tolist()
    blast_df = blast_df[~(blast_df.read.isin(reads_to_drop))]
    
    mutations_df = mutations_df[mutations_df.full_mutation.isin(recognized_mutations)]
    df = pd.merge(mutations_df, blast_df[['read']], how='right', on='read')
    df = df.sort_values('position')
    df['full_mutation'] = df.full_mutation.astype(str)
    df['mutations_on_read'] = df.groupby('read')['full_mutation'].transform(', '.join)
    df = df[['mutations_on_read', 'read']].drop_duplicates()
    df_counts = df.groupby('mutations_on_read').read.count().reset_index().rename(columns={'read':'read_count'}).sort_values('read_count', ascending=False)
    df_counts['read_percentage'] = df_counts.read_count * 100 / df_counts.read_count.sum()
    df_counts['mutations_on_read'] = df_counts.mutations_on_read.str.replace('nan', 'WT')
    
    # For every strain, calculate its percentage out of all the population containing the 
    # lowest appearing variant in the strain. This percentage can later be used as a cutoff
    # using the 90th percentile error frequency for example.
    recognized_dict = {}
    for i in recognized_mutations:
        recognized_dict[i] = df_counts[df_counts.mutations_on_read.str.contains(i)].read_percentage.sum()
    df_counts[['critical_variant_total_precent', 'critical_variant']] = df_counts.mutations_on_read.apply(lambda x: get_critical_variant_freq(x, recognized_dict)).apply(pd.Series)
    df_counts['percent_for_error_cutoff'] = df_counts.read_percentage / df_counts.critical_variant_total_precent
    df_counts.to_csv(args.output_file, index=False)
    return

def get_critical_variant_freq(mutations_on_read, recognized_dict):
    mutations_on_read = mutations_on_read.split(', ')
    smallest_variant = None
    smallest_freq = 100
    for m in recognized_dict:
        if m in mutations_on_read:
            if recognized_dict[m] < smallest_freq:
                smallest_freq = recognized_dict[m]
                smallest_variant = m
        else:
            if (100 - recognized_dict[m]) < smallest_freq:
                smallest_freq = (100 - recognized_dict[m])
                smallest_variant = m[:-1] + m[0]
    return (smallest_freq, smallest_variant)

    
if __name__ == "__main__":
    parser = argparse.ArgumentParser()
    parser.add_argument("-b", "--input_blast_df", type=str, help="path to blasts df csv", required=True)
    parser.add_argument('-m', '--input_mutation_df', type=str, help='path to mutations df csv', required=True)
    parser.add_argument('-p', '--input_chosen_mutation', type=str, help='path to csv file with mutations to check linkage of. Every mutation should have its own row, no header row, and the mutations should be written in the following format: "A1664.0G"', required=True)
    parser.add_argument("-o", "--output_file", type=str, help="a path to an output file", required=True)
    args = parser.parse_args()
    if not vars(args):
        parser.print_help()
        parser.exit(1)
    count_haplotypes(args)

