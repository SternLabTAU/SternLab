import os

import pandas as pd
import numpy as np
import glob
import matplotlib.pyplot as plt
import seaborn as sns;

from RG_HIVC_analysis.constants import excluded_samples

sns.set_context("poster")


def error_rates_summary(min_read_count = 100):
    freq_files_with_muts = glob.glob('/Users/omer/PycharmProjects/SternLab/RG_HIVC_analysis/ET86_2s/with_muts_coordinated/*.freqs')
    error_rates = pd.DataFrame(
        columns=['sample_id', 'syn', 'non_syn', 'stop'])

    i = 0
    for file in freq_files_with_muts:
        sample_id = os.path.splitext(os.path.basename(file))[0]
        if sample_id in excluded_samples:
            print('Excluded sample: {} - Skipping'.format(sample_id))
            continue
        print('Handling sample: ' + sample_id)
        freq_df = pd.read_csv(file, sep='\t')
        freq_df = freq_df[freq_df["Read_count"] > min_read_count]
        freq_df = freq_df.groupby('Mutation_type')
        current_error_rates = freq_df['Freq'].agg(['sum'])

        row = [sample_id] + [float(current_error_rates.loc['synonymous'])/9031] + [float(current_error_rates.loc['missense'])/9031] + [float(current_error_rates.loc['stop'])/9031]
        # print(row)
        error_rates.loc[i] = row
        i = i + 1

    print(error_rates)
    error_rates.to_csv('/Users/omer/PycharmProjects/SternLab/RG_HIVC_analysis/error_rates.csv', index=False)
    return error_rates



def error_rate_plots():
    pd.set_option('display.width', 600)  # TODO- remove
    pd.set_option('display.max_columns', 16)  # TODO- remove

    # joining diversity values with patient & date info
    samples_to_patient_and_dates = pd.read_csv('/Users/omer/PycharmProjects/SternLab/RG_HIVC_analysis/final_ZA04.csv',
                                sep=',')[['sample_id', 'ind_id', 'sample_date']]
    samples_to_patient_and_dates = samples_to_patient_and_dates.set_index('sample_id')
    error_rates = pd.read_csv('/Users/omer/PycharmProjects/SternLab/RG_HIVC_analysis/error_rates.csv', sep=',').set_index('sample_id')

    muts_by_ind = samples_to_patient_and_dates.join(error_rates)

    # sorting by patient + sampling date
    muts_by_ind['sample_date'] = pd.to_datetime(muts_by_ind['sample_date'], format='%d/%m/%Y')
    ordered_pis_by_ind = muts_by_ind.sort_values(by=['ind_id', 'sample_date'])

    # coverting "sample_date" to "time_since_infection"
    first_samples_dates = ordered_pis_by_ind.groupby('ind_id').first().reset_index()
    first_samples_dates = first_samples_dates[['ind_id', 'sample_date']]
    # print(first_samples_dates)

    ordered_pis_by_ind = ordered_pis_by_ind.merge(first_samples_dates,
                                        on='ind_id',
                                        how='left',
                                        sort=False,
                                        suffixes= ('','_r'))
    ordered_pis_by_ind['years_since_infection'] = (ordered_pis_by_ind['sample_date'] - ordered_pis_by_ind['sample_date_r']) / np.timedelta64(1, 'Y')
    print(ordered_pis_by_ind)

    # generating plot
    ordered_pis_by_ind = ordered_pis_by_ind.melt(id_vars= ('ind_id', 'years_since_infection'),
                        value_vars= ('syn', 'non_syn', 'stop'),
                        var_name='type',  value_name='muts_rate'
                        )
    ordered_pis_by_ind = ordered_pis_by_ind.sort_values(by=['ind_id', 'years_since_infection'])
    # print(ordered_pis_by_ind[ordered_pis_by_ind['ind_id'] == 16207])

    g = sns.relplot(
        x='years_since_infection',
        y='muts_rate',
        col='ind_id',
        hue='type',
        data=ordered_pis_by_ind,
        col_wrap=5,
        kind='line',
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
    error_rates_summary()
    error_rate_plots()
