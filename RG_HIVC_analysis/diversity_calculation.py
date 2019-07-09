import os

import pandas as pd
import glob
import numpy as np
import matplotlib.pyplot as plt
import seaborn as sns;

from RG_HIVC_analysis.coordinates import gag_ET86_interval, pol_ET86_interval, env_ET86_interval, excluded_samples

sns.set_context("poster")
import sys


def get_simple_diversity_stats():
    freq_files = glob.glob('/Users/omer/PycharmProjects/SternLab/RG_data_analysis/freq_files/*/*.freqs')
    basic_muts_count_stats = pd.DataFrame(
        columns=['sample_id', 'major_subs_count'])

    i = 0
    for file in freq_files:
        sample_id = str(file).split("/")[7]
        major_subs_count = count_major_subs(file)

        row = [sample_id] + [major_subs_count]
        # print(row)
        basic_muts_count_stats.loc[i] = row
        i = i + 1

    basic_muts_count_stats = basic_muts_count_stats.sort_values(by='major_subs_count', ascending=False)
    return basic_muts_count_stats


def count_major_subs(freq_file):
    df = pd.read_csv(freq_file, sep='\t')
    subs_count = df['Pos'].loc[(df['Read_count'] > 1000) & (df['Base'] != df['Ref']) & (df['Rank'] == 0) & (df['Base'] != '-') & (df['Ref'] != '-') ].count()
    return subs_count


def pis_calc(data, pivot_cols=[], min_read_count = 0, freq_threshold = 0, interval = (0, sys.maxsize)): # start_pos=0, end_pos= sys.maxsize
    """
    Calculates PI diversity per position, than calculates mean per group according to pivot_vols. Assumes data is not indexed.
    :param data:
    :param pivot_cols:
    :param min_read_count:
    :return:
    """
    def pairwise_differences_proportion(row):
        if row["Minor"] == 0:
            return 0
        total_part = (row["Total"] * (row["Total"] - 1))
        numerator = total_part - ((row["Major"] * (row["Major"] - 1)) + (row["Minor"] * (row["Minor"] - 1)))
        denominator = total_part
        return numerator * 1.0 / denominator

    # Filters
    # TODO for each filter: place here or extract from method?
    # transitions only
    data["mutation_type"] = data['Base'] + data['Ref']
    filtered_data = data[data["mutation_type"].isin(['GA','AG','GG','AA','CT','TC','CC','TT'])]
    # remove indels
    filtered_data = filtered_data[(filtered_data["Base"] != "-") & (filtered_data["Ref"] != "-")]
    # TODO- remove insertions
    # remove low coverage
    filtered_data = filtered_data[filtered_data["Read_count"] > min_read_count]
    # set low frequency to 0
    # filtered_data = filtered_data[filtered_data["Freq"] >= freq_threshold]
    filtered_data["Freq"] = np.where(filtered_data["Freq"] >= freq_threshold, filtered_data["Freq"], 0)
    # choose interval
    filtered_data = filtered_data[(filtered_data["Pos"] >= interval[0]) & (filtered_data["Pos"] <= interval[1])]

    if filtered_data.empty:
        # TODO - change to warning
        print('No relevant data after filtering. Skipping')
        return None

    filtered_data['counts_for_position'] = np.round(filtered_data['Read_count'] * filtered_data['Freq'])
    selecting_cols = pivot_cols[:]
    selecting_cols.extend(["Pos", "Base", "counts_for_position", "Rank"])
    filtered_data = filtered_data[selecting_cols]

    filtered_data["Rank"] = np.where(filtered_data["Rank"] == 0, "Major", "Minor")

    group_cols = pivot_cols[:]
    group_cols.extend(["Pos", "Rank"])

    # selecting max on cfp, per Pos\Rank (Major\minor?)- than will be graded per Pos, and summed
    # TODO: use all minor variants (not only max)
    filtered_data = filtered_data.groupby(group_cols)['counts_for_position'].aggregate(max).unstack().reset_index()
    if 'Minor' not in filtered_data.columns:
        return 0

    filtered_data["Total"] = filtered_data["Major"] + filtered_data["Minor"]

    filtered_data["pdp"] = filtered_data.apply(lambda row: pairwise_differences_proportion(row), axis=1)

    if any(pivot_cols):
        bysample_diversity = filtered_data.groupby(pivot_cols)['pdp'].agg(['count', 'sum']).reset_index()
        bysample_diversity["Pi"] = bysample_diversity["sum"] * 1.0 / bysample_diversity["count"]
        output_cols = pivot_cols[:]
        output_cols.append("Pi")
        pis = bysample_diversity[output_cols]
    else:
        pis = filtered_data['pdp'].mean()

    return pis


def pi_rates_summary():
    # freq_files = glob.glob('/Users/omer/PycharmProjects/SternLab/RG_data_analysis/freq_files_ZA04_2/*')
    freq_files = glob.glob('/Users/omer/PycharmProjects/SternLab/RG_data_analysis/ET86_2s/*.freqs')
    pi_diversity_rates = pd.DataFrame(
        columns=['sample_id', 'global', 'gag', 'pol', 'env'])

    # TODO- use append instead of loc[i]
    i = 0
    for file in freq_files:
        sample_id = os.path.splitext(os.path.basename(file))[0]
        if sample_id in excluded_samples:
            print('Excluded sample: {} - Skipping'.format(sample_id))
            continue
        print('Handling sample: ' + sample_id)

        freq_df = pd.read_csv(file, sep='\t')
        global_pi_rate = pis_calc(data=freq_df, min_read_count= 1000, freq_threshold= 0.01)

        # hxb2_file = glob.glob('/Users/omer/PycharmProjects/SternLab/RG_data_analysis/freq_files_HXB2_2/{}/*.freqs'.format(sample_id))[0]
        # freq_df_hxb2 = pd.read_csv(hxb2_file, sep='\t')
        # et86_file = glob.glob('/Users/omer/PycharmProjects/SternLab/RG_data_analysis/ET86_2s/{}.freqs'.format(sample_id))[0]
        # freq_df_et86 = pd.read_csv(et86_file, sep='\t')
        # global_pi_rate2 = pis_calc(data=freq_df_et86, min_read_count= 1000, freq_threshold= 0)

        gag_pi_rate = pis_calc(data=freq_df, min_read_count= 1000, freq_threshold= 0.01, interval= gag_ET86_interval)
        pol_pi_rate = pis_calc(data=freq_df, min_read_count= 1000, freq_threshold= 0.01, interval= pol_ET86_interval)
        env_pi_rate = pis_calc(data=freq_df, min_read_count= 1000, freq_threshold= 0.01, interval= env_ET86_interval)

        row = [sample_id] + [global_pi_rate] + [gag_pi_rate] + [pol_pi_rate] + [env_pi_rate]
        # print(row)
        pi_diversity_rates.loc[i] = row
        i = i + 1

    pi_diversity_rates = pi_diversity_rates.sort_values(by='pi_diversity_global', ascending=False)
    print(pi_diversity_rates)
    pi_diversity_rates.to_csv('/Users/omer/PycharmProjects/SternLab/RG_data_analysis/pi_rates_ET86_2.csv', index=False)
    return pi_diversity_rates

def main1():
    # pi_rates_summary()
    pi_diversity_plots()


def mergre_summary_tables():
    dates_vl = pd.read_csv('/Users/omer/PycharmProjects/SternLab/RG_data_analysis/dates_vl_stats.csv', sep=',')
    samples_format_conversion = pd.read_csv('/Users/omer/PycharmProjects/SternLab/RG_data_analysis/samples_format_conversion.csv', sep=',')
    pi_rates = pd.read_csv('/Users/omer/PycharmProjects/SternLab/RG_data_analysis/pi_rates_ZA04_2.csv', sep=',')
    coverage_stats = pd.read_csv('/Users/omer/PycharmProjects/SternLab/RG_data_analysis/coverage_stats_ET86_2.csv', sep=',')

    dates = dates_vl.set_index('FASTQ_name')
    conv = samples_format_conversion.set_index('table_fastq')
    join1 = dates.join(conv).set_index('sample_id')
    pi_rates = pi_rates.set_index('sample_id')
    cov = coverage_stats.set_index('sample_id')
    join2 = pi_rates.join(cov)

    final = join1.join(join2)
    final['sample_date'] = pd.to_datetime(final['sample_date'], format='%d/%m/%Y')
    final = final.sort_values(by=['ind_id', 'sample_date'])

    print(final)
    final.to_csv(path_or_buf='/Users/omer/PycharmProjects/SternLab/RG_data_analysis/final_ET86_pol_cov.csv')

    print(final.loc['130945_S2'])


def pi_diversity_plots():
    pd.set_option('display.width', 600)  # TODO- remove
    pd.set_option('display.max_columns', 16)  # TODO- remove

    # joining diversity values with patient & date info
    samples_to_patient_and_dates = pd.read_csv('/Users/omer/PycharmProjects/SternLab/RG_data_analysis/final_ZA04.csv',
                                sep=',')[['sample_id', 'ind_id', 'sample_date']]
    samples_to_patient_and_dates = samples_to_patient_and_dates.set_index('sample_id')
    pi_rates = pd.read_csv('/Users/omer/PycharmProjects/SternLab/RG_data_analysis/pi_rates_ET86_2.csv', sep=',').set_index('sample_id')

    pis_by_ind = samples_to_patient_and_dates.join(pi_rates)

    # sorting by patient + sampling date
    pis_by_ind['sample_date'] = pd.to_datetime(pis_by_ind['sample_date'], format='%d/%m/%Y')
    ordered_pis_by_ind = pis_by_ind.sort_values(by=['ind_id', 'sample_date'])

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
                        value_vars= ('global', 'gag', 'pol', 'env'),
                        var_name='regions',  value_name='pi_diversity'
                        )
    ordered_pis_by_ind = ordered_pis_by_ind.sort_values(by=['ind_id', 'years_since_infection'])
    # print(ordered_pis_by_ind[ordered_pis_by_ind['ind_id'] == 16207])

    g = sns.relplot(
        x='years_since_infection',
        y='pi_diversity',
        col='ind_id',
        hue='regions',
        data=ordered_pis_by_ind,
        col_wrap=5,
        kind='line',
        facet_kws={'sharex': True, 'legend_out':True},
        )
    g.set(yscale="log")
    g.set_ylabels("Pi diversity")
    g.set_xlabels("ET first sample (Years)")
    # g.set_xticklabels(rotation=45, fontsize=14)

    # extracting plot
    # plt.show()
    plt.savefig(fname= '/Users/omer/PycharmProjects/SternLab/RG_data_analysis/pi_trends_ET86_2.pdf')
    # g.savefig('/Users/omer/PycharmProjects/SternLab/RG_data_analysis/pi_rates_ET86_2.png')

def aggregation_tries():
    summary_table = pd.read_csv('/Users/omer/PycharmProjects/SternLab/RG_data_analysis/final_ZA04.csv',
                                               sep=',').set_index('sample_id')
    # a= summary_table.groupby('ind_id')['sample_date'].aggregate(min).unstack().reset_index()
    a = summary_table[['ind_id', 'sample_date']].groupby('ind_id').agg(lambda x: x.iloc[0]).set_index('ind_id')
    print(a)
    sfil= summary_table[['ind_id', 'sample_date', 'pi_diversity']].set_index('ind_id')
    join= sfil.join(a, lsuffix='sample_date', rsuffix='fisrt_date')
    print(join)

    # # print(summary_table)
    # # a = summary_table[summary_table["ind_id"] == 12796][['ind_id', 'sample_date', 'pi_diversity']]
    # a = summary_table[['ind_id', 'sample_date', 'pi_diversity']]
    # fd = a[['ind_id', 'sample_date']].groupby('ind_id').agg(lambda x: x.iloc[0])
    # print(a)
    # print(fd)


    # diversity_trends2 = summary_table.groupby('ind_id').apply(list)
    # diversity_trends_x = summary_table.groupby('ind_id')['sample_date'].apply(list)
    # # diversity_trends = summary_table.groupby('ind_id').agg({'sample_date':'list','pi_diversity':'sum'})


if __name__ == "__main__":
    main1()
