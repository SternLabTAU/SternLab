import glob
import os
import pandas as pd
from optparse import OptionParser
from freqs_utilities import add_mutation_to_freq_file

# This file requires clean-up & verification

def add_mutation_to_freq_file_with_global_ref(output_freq_file, freq_file, global_ref_file):
    print('Handling file:' + freq_file)
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


def main2():
    # freq_files = glob.glob('/Users/omer/PycharmProjects/SternLab/RG_data_analysis/ET86_2s/*.freqs')
    freq_files = glob.glob('/Users/omer/PycharmProjects/SternLab/RG_data_analysis/ET86_2s/504212_S56.freqs') #TODO- remove

    for file in freq_files:
        sample_id = os.path.splitext(os.path.basename(file))[0]
        print('Handling sample: ' + sample_id)

        # use global reference
        original_freq = pd.read_csv(file, sep='\t')
        global_ref = glob.glob('/Users/omer/PycharmProjects/SternLab/RG_data_analysis/C_ET86_ref.freqs')
        ref = pd.read_csv(global_ref[0], sep='\t')
        transformed_freq = inject_alternative_reference(original_freq, ref)

        # output = str(file).replace('.freqs','_muts.freqs')
        # output = str(file).replace('.freqs','_muts_combine_check.freqs') #TODO- remove
        output = f'/Users/omer/PycharmProjects/SternLab/RG_data_analysis/ET86_2s/with_muts_combined/{sample_id}.freqs'
        add_mutation_to_freq_file(output, freqs= transformed_freq, add_combined_major_mutation= False)
        add_mutation_to_freq_file(output, freqs= transformed_freq, add_combined_major_mutation= True)


#############################################


def create_mut_summaries_from_freq_with_muts():
    # freqs_with_muts = glob.glob('/Users/omer/PycharmProjects/SternLab/RG_data_analysis/ET86_2s/*_muts.freqs')
    freqs_with_muts = glob.glob('/Users/omer/PycharmProjects/SternLab/RG_data_analysis/ET86_2s/with_muts_combined/504193_S37.freqs') #TODO- remove

    for file in freqs_with_muts:
        sample_id = os.path.splitext(os.path.basename(file))[0]
        print('Handling sample: ' + sample_id)

        freq_df = pd.read_csv(file, sep='\t')
        # output = str(file).replace('.freqs','_summary.csv')
        output = str(file).replace('.freqs','_summary_check.csv') #TODO- remove
        aa_mutations_summary(freq_df, output)


def create_per_patient_mut_summaries():
    summary_table = pd.read_csv('/Users/omer/PycharmProjects/SternLab/RG_data_analysis/final_ZA04.csv',sep=',')
    summary_table['sample_date'] = pd.to_datetime(summary_table['sample_date'], format='%d/%m/%Y')
    summary_table = summary_table.sort_values(by=['ind_id', 'sample_date'])

    # for patient_id in [26892]:
    for patient_id in summary_table.ind_id:
        print('Patient: ' + str(patient_id))
        patient = summary_table.loc[summary_table.ind_id == patient_id][['sample_date','sample_id']]
        patient_sum = pd.DataFrame()

        for sample_id in patient.sample_id:
            print('Handling sample: ' + sample_id)
            file = glob.glob('/Users/omer/PycharmProjects/SternLab/RG_data_analysis/ET86_2s/with_muts_coordinated/{}.freqs'.format(sample_id))[0]
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
        patient_sum.to_csv(f'/Users/omer/PycharmProjects/SternLab/RG_data_analysis/ET86_2s/with_muts_coordinated/{patient_id}.csv', index=True, sep=',')


def aa_mutations_summary(freq_with_muts, output):
    # TODO - improve in next use
    # get Pol interval
    current = freq_with_muts
    pol_ET86_start_adjusted = 1453
    protease_ET86_start_adjusted = 1642 # start_aa_idx-63, length- 99 aa
    RT_ET86_start_adjusted = 1939 # start_aa_idx-162, length- 560 aa
    integrase_ET86_start_adjusted = 3622 # start_aa_idx-723, length- 288 aa
    current['AA_Pos'] = ((current.Pos - RT_ET86_start_adjusted - ((current.Pos - 1) % 3)) / 3)

    # filters
    current = current.loc[current.Mutation_type != 'consensus']
    current = current.loc[current.Mutation_type != 'synonymous']
    # current = current.loc[(current.Mutation_type == 'stop') | (current.Mutation_type == 'missense') | (current.CMM_type == 'stop') | (current.CMM_type == 'missense')] # TODO- cmm
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

    # current = current[['AA_Pos', 'wt_aa', 'mut_aa', 'Freq', 'Rank', 'CMM_aa', ]] #TODO- cmm
    # current['AA_Mut_summary'] = current.wt_aa + current.mut_aa + ' /' + current.Freq + ', m: ' + current.wt_aa + current.CMM_aa #TODO- cmm

    # print(current.loc[current.AA_Pos == '60'])
    # print(current)
    # print(current[['AA_Mut_summary']])

    current.to_csv(output, index=False, sep=',')
    return current


def main3():
    # create_freqs_with_muts()
    # create_mut_summaries_from_freq_with_muts() #sample summary
    create_per_patient_mut_summaries() # patient summary (loop on all patients)
    # TODO- compare mutation summaries


#############################################


if __name__ == "__main__":
    main2()


