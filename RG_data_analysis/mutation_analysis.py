import glob
import os

import pandas as pd
from optparse import OptionParser

from Bio.Seq import Seq

from freqs_utilities import add_mutation_to_freq_file


def create_freqs_with_muts():
    # freq_files = glob.glob('/Users/omer/PycharmProjects/SternLab/RG_data_analysis/ET86_2s/*.freqs')
    freq_files = glob.glob('/Users/omer/PycharmProjects/SternLab/RG_data_analysis/ET86_2s/504193_S37.freqs') #TODO- remove

    for file in freq_files:
        sample_id = os.path.splitext(os.path.basename(file))[0]
        print('Handling sample: ' + sample_id)

        # use global reference
        original_freq = pd.read_csv(file, sep='\t')
        global_ref = glob.glob('/Users/omer/PycharmProjects/SternLab/RG_data_analysis/C_ET86_ref.freqs')
        ref = pd.read_csv(global_ref[0], sep='\t')
        transformed_freq = inject_alternative_reference(original_freq, ref)

        output = str(file).replace('.freqs','_muts.freqs')
        add_mutation_to_freq_file(output, freqs= transformed_freq, combine_major_mutations = True)


def main():
    parser = OptionParser("usage: %prog [options]\nTry running %prog --help for more information")
    parser.add_option("-f", "--freqs", dest="freqs", help="frequency file")
    parser.add_option("-r", "--reference", dest="global_ref", help="global reference")
    parser.add_option("-o", "--output_freqs", dest="output_file", help="output freqs file")
    (options, args) = parser.parse_args()
    freq_file = options.freqs
    global_ref_file = options.global_ref
    output_freq_file = options.output_file

    add_mutation_to_freq_file_with_global_ref(output_freq_file, freq_file, global_ref_file)


def add_mutation_to_freq_file_with_global_ref(output_freq_file, freq_file, global_ref_file):
    print('Handling file:' + freq_file) # TODO- remove
    # use global reference
    original_freq = pd.read_csv(freq_file, sep='\t')
    ref = pd.read_csv(global_ref_file, sep='\t')
    transformed_freq = inject_alternative_reference(original_freq, ref)

    add_mutation_to_freq_file(output_freq_file, freqs= transformed_freq)


def inject_alternative_reference(original_freq, alternative_ref_freq):
    alternative_ref_freq = alternative_ref_freq.drop_duplicates("Pos")
    alternative_ref_freq = alternative_ref_freq[['Pos', 'Ref']]
    transformed_freq = original_freq.set_index(original_freq.Pos).join(alternative_ref_freq.set_index(alternative_ref_freq.Pos), rsuffix='_r')
    transformed_freq.Ref = transformed_freq.Ref_r
    transformed_freq = transformed_freq.loc[(transformed_freq.Pos >= 1357) & (transformed_freq.Pos <= 4587)][
        ['Pos', 'Base', 'Freq', 'Ref', 'Read_count', 'Rank', 'Prob']]

    return transformed_freq


#############################################


def create_mut_summaries_from_freq_with_muts():
    freqs_with_muts = glob.glob('/Users/omer/PycharmProjects/SternLab/RG_data_analysis/*_muts.freqs')

    for file in freqs_with_muts:
        sample_id = str(file).split("/")[7]
        print('Handling sample: ' + sample_id)

        freq_df = pd.read_csv(file, sep='\t')
        # output = str(file).replace('.freqs','_summary.csv')
        output = str(file).replace('.freqs','_summary_check.csv') #TODO- remove
        aa_mutations_summary(freq_df, output)

def create_per_patient_mut_summaries():
    summary_table = pd.read_csv('/Users/omer/PycharmProjects/SternLab/RG_data_analysis/final_ZA04.csv',sep=',')
    summary_table['sample_date'] = pd.to_datetime(summary_table['sample_date'], format='%d/%m/%Y')
    summary_table = summary_table.sort_values(by=['ind_id', 'sample_date'])

    # for patient_id in [13003]:
    for patient_id in summary_table.ind_id:
        print('Patient: ' + str(patient_id))
        patient = summary_table.loc[summary_table.ind_id == patient_id][['sample_date','sample_id']]
        patient_sum = pd.DataFrame()

        for sample_id in patient.sample_id:
            print('Handling sample: ' + sample_id)
            file = glob.glob('/Users/omer/PycharmProjects/SternLab/RG_data_analysis/ET86_2s/with_muts/{}.freqs'.format(sample_id))[0]
            freq_df = pd.read_csv(file, sep='\t')
            output = str(file).replace('.freqs','_summary.csv')
            aa_sum = aa_mutations_summary(freq_df, output)
            aa_sum = aa_sum[['AA_Pos', 'AA_Mut_summary']]
            aa_sum = aa_sum.groupby('AA_Pos')['AA_Mut_summary'].apply(lambda x: "%s" % ', '.join(x))
            # aa_sum.AA_Pos = aa_sum['AA_Pos'].astype(int) # TODO- handle

            if patient_sum.empty:
                patient_sum = aa_sum
            else:
                patient_sum = pd.merge(patient_sum,
                                        aa_sum,
                                        on='AA_Pos',
                                        how='outer',
                                        suffixes= ('','_' + sample_id))

        patient_sum.sort_values(by=['AA_Pos'])
        # print(patient_sum)
        patient_sum.to_csv('/Users/omer/PycharmProjects/SternLab/RG_data_analysis/ET86_2s/with_muts/{}.csv'.format(patient_id), index=True, sep=',')


def aa_mutations_summary(freq_with_muts, output):
    # get Pol interval
    current = freq_with_muts
    current['AA_Pos'] = ((current.Pos - 1453 - ((current.Pos - 1) % 3)) / 3)

    # filters
    current = current.loc[current.Mutation_type != 'consensus']
    current = current.loc[current.Mutation_type != 'synonymous']
    current = current.loc[(current.Pos >= 1456) & (current.Pos <= 4488)]
    current = current.loc[current.Freq > 0.1]
    current = current.loc[current.Read_count > 100]
    # current = current.loc[current.AA_Pos == 184]
    # current = current.loc[(current.AA_Pos >= 180) & (current.AA_Pos <= 190)]

    # formating
    current.AA_Pos = current.AA_Pos.astype(int)
    current.AA_Pos = current.AA_Pos.astype(str)
    current.Freq = current.Freq.round(3).apply(str)
    current = current[['AA_Pos', 'wt_aa', 'mut_aa', 'Freq', 'Rank']]
    # current['AA_Mut_summary'] = current.wt_aa + current.AA_Pos + current.mut_aa + ' /' + current.Freq # TODO- restore this line
    current['AA_Mut_summary'] = current.wt_aa + current.mut_aa + ' /' + current.Freq
    # print(current.loc[current.AA_Pos == '60'])
    # print(current)
    # print(current[['AA_Mut_summary']])

    current.to_csv(output, index=False, sep='\t')
    return current

def main2():
    # create_freqs_with_muts()
    create_mut_summaries_from_freq_with_muts()
    # TODO- compare mutation summaries

def add_mutation_to_freq_file_combine_mutations_by_rank(output, freqs_file = None, freqs = None, forced_rf_shift = 0):
    # assumes that position 1 is the beginning of the CDS
    # removes positions that at the beginning or at the end that are not part of a full codon
    if freqs_file == None and type(freqs) == "NoneType":
        raise Exception("Need to specify or freqs file path or a freqs pandas object")
    elif freqs_file != None and freqs != None:
        print(freqs_file, freqs)
        print(type(freqs))
        raise Exception("Need to specify EITHER freqs file path OR a freqs pandas object - only one!")
    elif freqs_file != None:
        freqs = pd.read_csv(freqs_file, sep="\t")
    freqs = freqs[freqs.Pos % 1 == 0] #removes insertions
    freqs = freqs[freqs.Base != "-"] #removes deletions
    freqs.reset_index(drop=True, inplace=True)

    first_pos = int(freqs.loc[1].Pos) #gets the first position in the right frameshift
    first_pos += forced_rf_shift
    if first_pos == 1:
        start_from = first_pos
    elif first_pos % 3 == 1:
        start_from = first_pos
    elif first_pos % 3 == 2:
        start_from = first_pos + 2
    elif first_pos % 3 == 0:
        start_from = first_pos + 1

    freqs["Mutation_type"] = None
    freqs["wt_aa"] = None
    freqs["mut_aa"] = None
    freqs["wt_codon"] = None
    freqs["mut_codon"] = None
    freqs["Mutation"] = None

    for pos in range(start_from, int(max(freqs.Pos)), 3): # add mutation information
        temp = freqs.loc[freqs['Pos'].isin([pos, pos+1, pos+2])]
        if len(temp) != 12: #12 - is 4 * 3 [A, C, G, T] (ordered by Rank)
            continue

        first = temp.iloc[0].Ref
        second = temp.iloc[4].Ref
        third = temp.iloc[8].Ref
        wt_codon = "".join([first, second, third])
        wt_aa = str(Seq(wt_codon).translate())

        pos = temp.iloc[0].Pos
        for n in range(8, 12):
            ref_base = temp.iloc[n].Ref
            mut_base = temp.iloc[n].Base
            mut_base_prev = temp.iloc[n-4].Base
            mut_base_prev_prev = temp.iloc[n-8].Base
            mut_codon = "".join([mut_base_prev_prev, mut_base_prev, mut_base])

            mut_aa = str(Seq(mut_codon).translate())

            if wt_codon == mut_codon:
                mutation_type = "consensus"
            elif wt_aa == mut_aa:
                mutation_type = "synonymous"
            elif wt_aa != "*" and mut_aa == "*":
                mutation_type = "stop"
            else:

                mutation_type = "missense"
            for i in range(0, 2):
                freqs.loc[(freqs["Pos"] == pos + i) & (freqs["Base"] == mut_base), "Mutation_type"] = mutation_type
                freqs.loc[(freqs["Pos"] == pos + i) & (freqs["Base"] == mut_base), "wt_aa"] = wt_aa
                freqs.loc[(freqs["Pos"] == pos + i) & (freqs["Base"] == mut_base), "mut_aa"] = mut_aa
                freqs.loc[(freqs["Pos"] == pos + i) & (freqs["Base"] == mut_base), "wt_codon"] = wt_codon
                freqs.loc[(freqs["Pos"] == pos + i) & (freqs["Base"] == mut_base), "mut_codon"] = mut_codon
                freqs.loc[(freqs["Pos"] == pos + i) & (freqs["Base"] == mut_base), "Mutation"] = ref_base + mut_base

    freqs = freqs[freqs.Mutation_type.notnull()] #removes Nones - rows at the beginning and the end
    freqs.to_csv(output, index=False, sep='\t')
    return freqs

if __name__ == "__main__":
    main2()
