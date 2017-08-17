import pandas as pd
import process_freqs
import os
import argparse
from collections import Counter


''' This script parses a unite freq file into a FITS input file'''


def main(args):
    df = pd.read_csv(args.in_file, sep='\t')
    df = remove_del_ins(df)

    # remove 18, 20 from calculation - for MS2
    df = df[(df['Time'] != 18) & (df['Time'] != 20)]

    # filter according line if needed
    if args.line != None:
        degree = int(args.line[:2])
        replica = args.line[-1]
        df = df[(df['Degree'] == degree) & (df['Replica'] == replica)]

    # create FITS input files
    split_by_mutation(df, args.out, ID=args.line)




def remove_del_ins(df):
    """
    This method receives a data frame of a freq file including a Time column derived out of ID column
    and removes insertions (X.x) and deletions ('-')
    :param df: a data frame of a freq file containing information about passage (Time columns needed!)
    :return: a data frame without deletions and insertions
    """

    # remove insertions
    df = process_freqs.remove_dot_from_pos(df)
    # remove deletions
    df = df[df['Base'] != '-']

    return df



def sum_2_one(df):
    """
    This method receives a data frame and make the wt frequency to be 1 - mutant frequency
    :param df: a data frame
    :return:
    """
    df = df.sort_values(['Time','Pos', 'Base'])
    frequencies = df['Freq'].values
    fit_frequencies = []
    # take the mutant frequency which is in indices 1,3,5,7.. and reduce from 1 to get the non mutant frequency
    for i in range(0, len(frequencies), 2):
        fit_frequencies.extend([1 - frequencies[i+1], frequencies[i+1]])

    # change the frequency column to the new one after fitting
    df['Freq'] = fit_frequencies
    return df


def split_by_mutation(df, out, ID=None):
    """
    This method receives a data frame and splits it to data frames containing mutation type
    :param df: a data frame containing all frequencies
    :param out: a directory to save the results
    :param ID: an id that will be added to the filename, optional, default is none.
    :return: saves the input files to out directory
    """

    mutations = ['AA', 'AG', 'GG', 'GA', 'CC', 'CT', 'TT', 'TC']
    df['Mutation'] = df['Ref'] + df['Base']
    for i in range(0, len(mutations), 2):
        print("Starting splitting mutation type {}....\n".format(mutations[i+1]))
        # filter according to mutation and wt
        mut_df = df[(df['Mutation'] == mutations[i]) | (df['Mutation'] == mutations[i+1])]
        print(set(mut_df.Mutation.values))

        # change bases to numeric value 0 for wt and 1 for mutant
        mut_df['Base'] = 1 - (mut_df['Base'] == mut_df['Ref']).astype(int)
        mut_df['Ref'] = 0
        mut_df.reset_index(drop=True, inplace=True)

        # DEBUG
        # for p in range(0, len(mut_df['Base'].values), 2):
        #     if mut_df['Base'].values[p] != 0 or mut_df['Base'].values[p+1] != 1:
        #         print(mut_df.ix[p])

        count = Counter(mut_df['Base'].values)
        if count[0] != count[1]:
            raise Exception('Error: incompatible fits format. Base to index failed\n')

        # change the frequencies so the sum of each mutant wt will be 1
        mut_df = sum_2_one(mut_df)

        # remove unnecessary columns and change the order of the columns
        mut_df = mut_df.drop(['Prob', 'Rank', 'Replica', 'Sample', 'Degree', 'Mutation'], axis=1)
        mut_df.rename(columns={"Time": "Gen"}, inplace=True)
        mut_df = mut_df[['Gen', 'Base', 'Freq', 'Ref', 'Read_count', 'Pos']]


        # save the result into a file
        if ID != None:
            filename = 'FITS_input_file_{}_{}'.format(ID, mutations[i+1]) # mutation type is at i + 1
        else:
            filename = 'FITS_input_file_{}'.format(mutations[i + 1])  # mutation type is at i + 1

        print("Done parsing mutation type {}. Saving..\n".format(mutations[i+1]))
        mut_df.to_csv(os.path.join(out, filename), sep ='\t', encoding='utf-8', index = False)

if __name__ == "__main__":
    parser = argparse.ArgumentParser()
    parser.add_argument("-i", "--in_file", type=str, help="a path to an input frequency file with different time points", required=True)
    parser.add_argument("-o", "--out", type=str, help="output directory to save the files", required=True)
    parser.add_argument("-l", "--line", type=str, help="an ID of the experiment", required=False)
    args = parser.parse_args()
    main(args)
