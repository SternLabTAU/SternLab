import glob
import os

import numpy as np
import pandas as pd

from RG_HIVC_analysis import constants


def convert_freqs_to_zn_input_format():
    freq_files = glob.glob('/Users/omer/PycharmProjects/SternLab/RG_HIVC_analysis/ET86_2s/*.freqs')
    # freq_files = glob.glob('/Users/omer/PycharmProjects/SternLab/RG_HIVC_analysis/ET86_2s/100888_S14.freqs')

    for file in freq_files:
        freq_df = pd.read_csv(file, sep='\t')
        freq_df['counts_for_position'] = np.round(freq_df['Read_count'] * freq_df['Freq'])
        freq_df = freq_df[['Pos', 'Base', 'counts_for_position']]
        freq_df = freq_df[freq_df.Pos % 1 == 0]  # removes insertions
        freq_df.reset_index(drop=True, inplace=True)

        df = freq_df.pivot(index='Pos', columns='Base', values='counts_for_position')

        # re-order columns
        cols = df.columns.tolist()
        cols = cols[1:] + cols[:1]
        df = df[cols]
        # add N column
        df['N'] = 0.0
        # print(df.iloc[-1])
        # print(df)

        # filling missing indices (with zeros)
        idx_reference = pd.DataFrame({'Pos': range(1, constants.ET86_length + 1)})
        complete_idx = idx_reference.merge(df, on='Pos', how='left', sort=False,suffixes=('', '_r'))
        complete_idx.fillna(0)
        # print(complete_idx)
        cols = complete_idx.columns.tolist()
        counts_df = complete_idx[cols[1:]]

        # convert to numpy
        counts_array = counts_df.to_numpy().astype(int).transpose()
        # print(counts_array.shape)

        output_folder = '/Users/omer/PycharmProjects/SternLab/RG_HIVC_analysis/ZN_input/single_nucleotide_variants/'
        npy_file_name = os.path.splitext(os.path.basename(file))[0] + '.npy'
        path = output_folder + npy_file_name
        print(path)
        np.save(path, counts_array)


def generate_samples_tables():
    # visibility
    pd.set_option('display.width', 600)
    pd.set_option('display.max_columns', 16)

    # fetching all patients data
    summary_table = pd.read_csv('/Users/omer/PycharmProjects/SternLab/RG_HIVC_analysis/final_ZA04.csv', sep=',')
    # ** FILTERING low quality SAMPLES **
    summary_table = summary_table[~summary_table['sample_id'].isin(constants.excluded_samples)]
    # sorting by patient + sampling date
    summary_table['sample_date'] = pd.to_datetime(summary_table['sample_date'], format='%d/%m/%Y')
    summary_table = summary_table.sort_values(by=['ind_id', 'sample_date'])
    # adding "days_since_infection"
    first_samples_dates = summary_table.groupby('ind_id').first().reset_index()
    first_samples_dates = first_samples_dates[['ind_id', 'sample_date']]
    summary_table = summary_table.merge(first_samples_dates,on='ind_id',  how='left',  sort=False, suffixes= ('','_r'))
    summary_table['days since infection'] = (summary_table['sample_date'] - summary_table['sample_date_r']) / np.timedelta64(1, 'D')
    summary_table['days since infection'] = summary_table['days since infection'].astype(int)

    # for patient_id in [26892]:
    for patient_id in set(summary_table.ind_id):
        print('Patient: ' + str(patient_id))
        samples_table = summary_table.loc[summary_table.ind_id == patient_id][['sample_id', 'ind_id', 'days since infection', 'sample_VL']]
        samples_table.rename({'sample_id': 'id', 'ind_id': 'patient', 'sample_VL': 'viral load'}, axis='columns', inplace = True)
        # print(samples_table)
        samples_table.to_csv(f'/Users/omer/PycharmProjects/SternLab/RG_HIVC_analysis/ZN_input/tables/samples_{patient_id}.tsv',index=False, sep='\t')

        # for sample_id in samples_table.sample_id:


if __name__ == "__main__":
    convert_freqs_to_zn_input_format()
    # aft = np.load('/Users/omer/PycharmProjects/SternLab/RG_HIVC_analysis/ZN_input/single_nucleotide_variants/100888_S14.npy')
    # print(len(aft.argmax(axis=0)))
    # generate_samples_tables()
