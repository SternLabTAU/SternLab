import os

import pandas as pd
import glob
import numpy as np
import matplotlib.pyplot as plt
import seaborn as sns; sns.set_context("poster")
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
    # remove low coverage
    filtered_data = filtered_data[filtered_data["Read_count"] > min_read_count]
    # remove low frequency
    filtered_data = filtered_data[filtered_data["Freq"] >= freq_threshold]
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
    freq_files = glob.glob('/Users/omer/PycharmProjects/SternLab/RG_data_analysis/freq_files_ZA04_2/*')
    # freq_files = glob.glob('/Users/omer/PycharmProjects/SternLab/RG_data_analysis/ET86_2s/*.freqs')
    pi_diversity_rates = pd.DataFrame(
        columns=['sample_id', 'pi_diversity_global', 'pi_diversity_gag', 'pi_diversity_pol', 'pi_diversity_env'])
    gag_ET86_interval = (170, 1684)
    pol_ET86_interval = (1456, 4488)
    env_ET86_interval = (5637, 8196)
    excluded_samples = ('X84335_S20','504214_S58','504184_S28','504186_S30','84864_S47','504206_S50','504190_S34','504191_S35','504192_S36','504198_S42','X84434_S3','X145364-R_S95')

    # TODO- use append instead of loc[i]
    i = 0
    for file in freq_files:
        sample_id = os.path.splitext(os.path.basename(file))[0]
        if sample_id in excluded_samples:
            print('Excluded sample: {} - Skipping'.format(sample_id))
            continue
        print('Handling sample: ' + sample_id)

        freq_df = pd.read_csv(file, sep='\t')
        global_pi_rate = pis_calc(data=freq_df, min_read_count= 1000, freq_threshold= 0)

        hxb2_file = glob.glob('/Users/omer/PycharmProjects/SternLab/RG_data_analysis/freq_files_HXB2_2/{}/*.freqs'.format(sample_id))[0]
        freq_df_hxb2 = pd.read_csv(hxb2_file, sep='\t')
        global_pi_rate2 = pis_calc(data=freq_df_hxb2, min_read_count= 1000, freq_threshold= 0)

        # global_pi_rate2 = pis_calc(data=freq_df, min_read_count= 100, freq_threshold= 0.01)
        # global_pi_rate3 = pis_calc(data=freq_df, min_read_count= 100, freq_threshold= 0.01)
        # global_pi_rate4 = pis_calc(data=freq_df, min_read_count= 1000, freq_threshold= 0.01)
        # gag_pi_rate = pis_calc(data=freq_df, min_read_count= 100, freq_threshold= 0.01, interval= gag_ET86_interval)
        # pol_pi_rate = pis_calc(data=freq_df, min_read_count= 100, freq_threshold= 0.01, interval= pol_ET86_interval)
        # env_pi_rate = pis_calc(data=freq_df, min_read_count= 100, freq_threshold= 0.01, interval= env_ET86_interval)

        # row = [sample_id] + [global_pi_rate] + [gag_pi_rate] + [pol_pi_rate] + [env_pi_rate]
        # row = [sample_id] + [global_pi_rate] + [global_pi_rate2] + [global_pi_rate3] + [global_pi_rate4]
        row = [sample_id] + [global_pi_rate] + [global_pi_rate2] + [0] + [0]
        # print(row)
        pi_diversity_rates.loc[i] = row
        i = i + 1

    pi_diversity_rates = pi_diversity_rates.sort_values(by='pi_diversity_global', ascending=False)
    print(pi_diversity_rates)
    pi_diversity_rates.to_csv('/Users/omer/PycharmProjects/SternLab/RG_data_analysis/pi_rates_ZA04_check.csv', index=False)
    return pi_diversity_rates

def main1():
    pi_rates_summary()
    # pi_diversity_plots()


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
    samples_to_patient_and_dates = pd.read_csv('/Users/omer/PycharmProjects/SternLab/RG_data_analysis/final_ZA04.csv',
                                sep=',').set_index('sample_id')[['sample_id', 'ind_id', 'sample_date']]
    pi_rates = pd.read_csv('/Users/omer/PycharmProjects/SternLab/RG_data_analysis/pi_rates_ET86.csv', sep=',').set_index('sample_id')

    pis_with_dates = samples_to_patient_and_dates.join(pi_rates)

    pis_with_dates['sample_date'] = pd.to_datetime(pis_with_dates['sample_date'], format='%d/%m/%Y')
    pis_with_dates = pis_with_dates.sort_values(by=['ind_id', 'sample_date'])

    # TODO- switch "sample_date" to "time_since_infection"
    # pis_with_dates["time_since_infection"] = pis_with_dates.groupby('ind_id')['sample_date'].aggregate(min).unstack()

    pis_with_dates.melt(id_vars= ('ind_id', 'sample_date'),
                        value_vars= ('pi_diversity_global', 'pi_diversity_gag', 'pi_diversity_pol', 'pi_diversity_env'),
                        var_name='regions',  value_name='Pi diversity'
                        )
    pis_with_dates = pis_with_dates.sort_values(by=['ind_id', 'sample_date'])

    g = sns.relplot(
        x='sample_date',
        y='pi_diversity',
        col='ind_id',
        hue='regions',
        data=pis_with_dates,
        col_wrap=5,
        facet_kws={'sharex': False},
        kind='line')
    g.set(yscale="log")
    g.set_xticklabels(rotation=45, fontsize=14)

    # plt.savefig(fname= '/Users/omer/PycharmProjects/SternLab/RG_data_analysis/pi_trends_ET86.pdf')
    plt.show()

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
