#! /usr/local/python_anaconda/bin/python3.4


import pandas as pd
import numpy as np
import os
import scipy.stats
import argparse


def association_test(args):
    '''
    This function runs association tests for position couples with index matching
    the job id. The function gets the blast dataframe path, the mutations dataframe
    path, the association coupe index file and a directory to save the results to.
    '''
    blast_df = pd.read_csv(args.input_blast_df)
    # choose only reads that were mapped only once in blast
    blast_df['read_count'] = blast_df.groupby('read')['start_ref'].transform('count')
    blast_df = blast_df[(blast_df.read_count == 1)]
    
    mutations_df = pd.read_csv(args.input_mutation_df)
    association_test_dir = args.output_dir + '/' + str(args.pbs_job_array_id) + '/'
    # create directory for specific job id if does not exist already
    if not os.path.exists(association_test_dir):
        os.mkdir(association_test_dir)    
    couples = pd.read_csv(args.input_position_pairs_df)
    # get position couples with the association index the same as the job id.
    position_tuple_list = list(couples[(couples.association_index == int(args.pbs_job_array_id))][['i','j']].itertuples(index=False, name=None))
    
    for (i, j) in position_tuple_list:
        if not os.path.isfile(association_test_dir + str(i) + '_' + str(j)) and not os.path.isfile(association_test_dir + str(j) + '_' + str(i)):
            print((i,j))
            temp_matrix = pd.DataFrame(np.zeros(shape=(2,2)), columns=[j,0], index=[i,0])
            b = blast_df.copy()
            m = mutations_df.copy()
            # look only at mutations in positions i and j
            m = m[(m.position.isin([i,j]))]
            # look only at alignments that contain both positions i and j
            b = b[(b.start_ref < i) & (b.end_ref > i) & (b.start_ref < j) & (b.end_ref > j)]
            b = b[['read']]
            # keep all reads that contain positions i and j (both containing mutations in i and/or j 
            # and not containing mutations there), and their appropriate mutations.
            relevant_mutations = pd.merge(b, m, on='read', how='left').fillna(0)
            
            # from reads containing both positions, count reads with mutation in i, with mutation 
            # in j, with mutations in both or with no mutations and create matrix for association test.
            reads_count = len(relevant_mutations[['read']].drop_duplicates())
            reads_with_mut_i_mut_j  = len(pd.merge(relevant_mutations[(relevant_mutations.position == i)], relevant_mutations[(relevant_mutations.position == j)], on='read', how='inner').read.drop_duplicates())
            reads_with_mut_j = len(relevant_mutations[(relevant_mutations.position == j)][['read']].drop_duplicates())
            reads_with_mut_i = len(relevant_mutations[(relevant_mutations.position == i)][['read']].drop_duplicates())
            reads_with_mut_i_wt_j = reads_with_mut_i - reads_with_mut_i_mut_j
            reads_with_wt_i_mut_j = reads_with_mut_j - reads_with_mut_i_mut_j
            reads_with_wt_i_wt_j = reads_count - reads_with_mut_i_mut_j - reads_with_mut_i_wt_j - reads_with_wt_i_mut_j
            temp_matrix.at[i,j] = reads_with_mut_i_mut_j
            temp_matrix.at[i,0] = reads_with_mut_i_wt_j
            temp_matrix.at[0,j] = reads_with_wt_i_mut_j
            temp_matrix.at[0,0] = reads_with_wt_i_wt_j
            # run association test on matrix and write results to file.
            if temp_matrix.sum(axis=0).all() > 0 and temp_matrix.sum(axis=1).all() > 0:
                a = (scipy.stats.chi2_contingency(temp_matrix))
                with open(association_test_dir + str(i) + '_' + str(j), 'w') as f:
                    f.write('\n'.join([str(i) for i in a]))
                temp_matrix.to_csv(association_test_dir + str(i) + '_' + str(j) + '.csv')
    return


if __name__ == "__main__":
    parser = argparse.ArgumentParser()
    parser.add_argument("-b", "--input_blast_df", type=str, help="path to blasts df csv", required=True)
    parser.add_argument('-m', '--input_mutation_df', type=str, help='path to mutations df csv', required=True)
    parser.add_argument('-i', '--input_position_pairs_df', type=str, help='path to position pairs index df, created by create_couples_index_file.py', required=True)
    parser.add_argument("-o", "--output_dir", type=str, help="a path to an output directory", required=True)
    parser.add_argument('-j', '--pbs_job_array_id', type=int, help='pbs job array id, used to choose pairs of positions to check.', required=True)
    args = parser.parse_args()
    if not vars(args):
        parser.print_help()
        parser.exit(1)
    association_test(args)