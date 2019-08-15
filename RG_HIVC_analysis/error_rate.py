import os

import pandas as pd
import numpy as np
import glob
import matplotlib.pyplot as plt
import seaborn as sns;

from RG_HIVC_analysis.constants import excluded_samples

sns.set_context("poster")

def create_unified_samples_to_patient_and_dsi():
    extension = 'tsv'
    all_filenames = [i for i in glob.glob('/Users/omer/PycharmProjects/SternLab/RG_HIVC_analysis/ZN_input/tables/samples_*.{}'.format(extension))]

    # combine all files in the list
    combined_csv = pd.concat([pd.read_csv(f) for f in all_filenames])
    # export to csv
    combined_csv.to_csv("samples_to_patient_and_dsi.csv", index=False, encoding='utf-8-sig')


def generate_error_rates_dataframe(min_read_count = 100):
    freq_files_with_muts = glob.glob('/Users/omer/PycharmProjects/SternLab/RG_HIVC_analysis/ET86_2s/with_muts_coordinated/*.freqs')
    # freq_files_with_muts = glob.glob('/Users/omer/PycharmProjects/SternLab/RG_HIVC_analysis/ET86_2s/with_muts_coordinated/504201_S45.freqs')

    freq_dfs_with_id = []
    samples_to_patient_and_dates = pd.read_csv('/Users/omer/PycharmProjects/SternLab/RG_HIVC_analysis/samples_to_patient_and_dsi.csv', sep='\t')
    samples_to_patient_and_dates.set_index('id', inplace= True)

    for file in freq_files_with_muts:
        sample_id = os.path.splitext(os.path.basename(file))[0]
        if sample_id in excluded_samples:
            print('Excluded sample: {} - Skipping'.format(sample_id))
            continue
        print('Handling sample: ' + sample_id)

        freq_df = pd.read_csv(file, sep='\t')
        # filtering freq file
        freq_df = freq_df[freq_df["Read_count"] > min_read_count]

        #adding id & dsi columns
        patient_id = samples_to_patient_and_dates.loc[f'{sample_id}', 'patient']
        days_since_infection = samples_to_patient_and_dates.loc[f'{sample_id}', 'days since infection']
        freq_df['ind_id'] = patient_id
        freq_df['years_since_infection'] = str(float(days_since_infection) / float(365))

        # print(freq_df)
        freq_dfs_with_id.append(freq_df)


    unified_freq_df_with_ids_ysi = pd.concat(freq_dfs_with_id)
    print(unified_freq_df_with_ids_ysi.shape)
    return unified_freq_df_with_ids_ysi



def error_rate_plots(unified_freq_df):
    pd.set_option('display.width', 600)  # TODO- remove
    pd.set_option('display.max_columns', 16)  # TODO- remove

    g = sns.catplot(
        x='years_since_infection',
        y='muts_rate',
        col='ind_id',
        hue='type',
        data=unified_freq_df,
        col_wrap=5,
        kind='box',
        facet_kws={'sharex': True, 'legend_out':True},
        )
    g.set(yscale="log")
    g.set_ylabels("Muts rate")
    g.set_xlabels("ET first sample (Years)")
    # g.set_xticklabels(rotation=45, fontsize=14)

    # extracting plot
    # plt.show()
    plt.savefig(fname= '/Users/omer/PycharmProjects/SternLab/RG_HIVC_analysis/figures/muts_rate_trends.png')
    # g.savefig('/Users/omer/PycharmProjects/SternLab/RG_HIVC_analysis/figures/pi_rates_ET86_2.png')



if __name__ == "__main__":
    unified_freq_df = generate_error_rates_dataframe()
    unified_freq_df.to_csv('/Users/omer/PycharmProjects/SternLab/RG_HIVC_analysis/ET86_2s/with_muts_coordinated/unified_freqs_with_ids_ysi.csv', index=False)
    # unified_freq_df = pd.read_csv('/Users/omer/PycharmProjects/SternLab/RG_HIVC_analysis/ET86_2s/with_muts_coordinated/unified_freqs_with_ids_ysi.csv', sep='\t')
    error_rate_plots(unified_freq_df)
