import os

import pandas as pd
import numpy as np
import glob
import matplotlib.pyplot as plt
import seaborn as sns;
from scipy import stats

from RG_HIVC_analysis.constants import excluded_samples

sns.set_context("poster")

def create_unified_samples_to_patient_and_dsi():
    extension = 'tsv'
    all_filenames = [i for i in glob.glob('/Users/omer/PycharmProjects/SternLab/RG_HIVC_analysis/ZN_input/tables/samples_*.{}'.format(extension))]

    # combine all files in the list
    combined_csv = pd.concat([pd.read_csv(f) for f in all_filenames])
    # export to csv
    combined_csv.to_csv("samples_to_patient_and_dsi.csv", index=False, encoding='utf-8-sig')


def generate_unified_filtered_verbose_freqs_df(min_read_count = 100, freq_threshold=0.01):
    folder_path = '/Users/omer/PycharmProjects/SternLab/RG_HIVC_analysis/ET86_4s/with_muts_and_con_as_ref'
    freq_files_with_muts = glob.glob(f'{folder_path}/*.freqs')

    freq_dfs_with_id = []
    samples_to_patient_and_dates = pd.read_csv('/Users/omer/PycharmProjects/SternLab/RG_HIVC_analysis/output_tables/samples_to_patient_and_dsi.csv', sep='\t')
    samples_to_patient_and_dates.set_index('id', inplace= True)

    for file in freq_files_with_muts:
        sample_id = os.path.splitext(os.path.basename(file))[0]

        # filtering samples
        if sample_id in excluded_samples:
            print('Excluded sample: {} - Skipping'.format(sample_id))
            continue
        print('Handling sample: ' + sample_id)
        freq_df = pd.read_csv(file, sep='\t')

        # filtering freq file
        # TODO- additional filters?
        freq_df = freq_df[freq_df["Read_count"] > min_read_count]
        freq_df['Freq'] = np.where(freq_df['Freq'] >= freq_threshold, freq_df['Freq'], 0)
        freq_df = freq_df[freq_df['Mutation_type'] != 'consensus']

        # adding id & dsi columns
        patient_id = samples_to_patient_and_dates.loc[f'{sample_id}', 'patient']
        days_since_infection = samples_to_patient_and_dates.loc[f'{sample_id}', 'days since infection']
        freq_df['ind_id'] = patient_id
        freq_df['sample_id'] = sample_id
        freq_df['years_since_infection'] = str(np.round((days_since_infection) / float(365),2))

        # if patient_id == 29447:
        # print(freq_df)
        freq_dfs_with_id.append(freq_df)


    unified_freq_df_with_ids_ysi = pd.concat(freq_dfs_with_id)
    print(unified_freq_df_with_ids_ysi.shape)
    unified_freq_df_with_ids_ysi.to_csv( f'{folder_path}/unified_freqs_filtered_verbose.csv', index=False)
    return unified_freq_df_with_ids_ysi


def plot_mutation_rate_distribution(unified_freq_df):
    pd.set_option('display.width', 600)  # TODO- remove
    pd.set_option('display.max_columns', 16)  # TODO- remove

    # g = sns.catplot(
    #     x='years_since_infection',
    #     y='Freq',
    #     col='ind_id',
    #     hue='Mutation_type',
    #     data=unified_freq_df,
    #     col_wrap=5,
    #     kind='box',
    #     facet_kws={'sharex': False, 'legend_out':True},
    #     )

    g = sns.boxplot(x='ind_id', y='Freq',hue='Mutation_type', data=unified_freq_df)

    # plot adjustments
    # g.set(yscale="log")
    plt.yscale('log')
    # plt.ylim(10**-6, 1)
    # g.set_ylabels("Muts rate")
    # g.set_xlabels("ET first sample (Years)")
    # g.set_xticklabels(rotation=45, fontsize=14)

    # extract plot
    plt.show()
    # plt.savefig(fname= '/Users/omer/PycharmProjects/SternLab/RG_HIVC_analysis/figures/muts_rates_distribution.png')
    # g.savefig('')


def plot_mutation_rate_conf_interval(unified_freq_df):
    add_weighted_freq(unified_freq_df)

    def weighted_varaint(x, **kws):
        var, count = map(np.asarray, zip(*x))
        return var.sum() / count.sum()

    # basic presentation
    # g = sns.factorplot(x="sample_id",
    #                    y="count_and_weight",
    #                    data=unified_freq_df,
    #                    hue="Mutation_type",
    #                    hue_order=["missense", "synonymous", "stop"],
    #                    # palette="tab20",
    #                    join=False,
    #                    orient="v",
    #                    estimator=weighted_varaint,
    #                    dodge=1)

    # advanced, chronological presentation
    g = sns.catplot(x= "years_since_infection",
                       y= "count_and_weight",
                       hue="Mutation_type",
                       hue_order=["missense", "synonymous", "stop"],
                       # palette="tab20",
                       col='ind_id',
                       col_wrap=5,
                       n_boot=1000,
                       # join=True,
                       orient="v",
                       kind='point',
                       estimator=weighted_varaint,
                       dodge=1,
                       data=unified_freq_df)

    # TODO- optional:
    # gateway to value extraction
    # ax = sns.pointplot(x= "ind_id",
    #                    y= "count_and_weight",
    #                    data=unified_freq_df,
    #                    hue="Mutation_type",
    #                    # palette="tab20",
    #                    join=False,
    #                    orient="v",
    #                    estimator=weighted_varaint,
    #                    dodge=1)


    # plot adjustments
    g.set(yscale="log")
    # plt.ylim(10**-6, 1)
    g.set_ylabels("mutation_rate")
    # g.set_xlabels("ET first sample (Years)")
    # g.set_xticklabels(rotation=45, fontsize=11)

    # extract plot
    # plt.savefig(fname= '/Users/omer/PycharmProjects/SternLab/RG_HIVC_analysis/figures/mutation_rates_ET86_4s.png')
    plt.show()
    # g.savefig('')


def add_weighted_freq(unified_freq_df):
    unified_freq_df["count_for_position"] = unified_freq_df["Freq"] * unified_freq_df["Read_count"]
    # unified_freq_df["count_for_position"] = np.where(unified_freq_df["Prob"] < 0.95, 0, unified_freq_df["no_variants"]) #TODO- relevant?
    unified_freq_df["count_and_weight"] = list(zip(unified_freq_df.count_for_position, unified_freq_df.Read_count))

    # unified_freq_df.to_csv('/Users/omer/PycharmProjects/SternLab/RG_HIVC_analysis/ET86_2s/with_muts_coordinated/examine_freqs_with_weights.csv', index=False)
    print(unified_freq_df.shape)


def verify_chosen_stats(unified_freq_df):
    pd.set_option('display.max_rows', 250)
    chosen_stat = unified_freq_df[(unified_freq_df['ind_id'] == 22097) & (unified_freq_df['Mutation_type'] == 'synonymous')]
    # The mut measure is defined as rate per position, that's why avreage on pos is ok.
    # aggregation on all combinations of positions+muts able to generate specific kind of mutation
    # what if there's only one pair pos+muts that can generate stop mutation? its value is the mut value? and if there's 2?
    #          -> final value we generate is a repesentitive of what?
    # total coverage collected on those positions- serves as sanity measure
    chosen_stat = chosen_stat[chosen_stat['sample_id'] == '504188_S32']
    # chosen_stat = chosen_stat[chosen_stat['sample_id'] == '504211_S55']
    chosen_stat = chosen_stat[['Freq', 'Read_count']]
    chosen_stat['count_for_position'] = chosen_stat['Freq'] * chosen_stat['Read_count']

    print(chosen_stat)

    print('mean is {}'.format(chosen_stat['Freq'].mean()))
    print('weighted_mean is {}'.format(chosen_stat['count_for_position'].sum()/ chosen_stat['Read_count'].sum()))


def potential_analysis(unified_freq_df):
    pd.set_option('display.width', 600)
    pd.set_option('display.max_columns', 16)
    pd.set_option('display.max_rows', 250)
    potential_analysis = unified_freq_df.groupby(['sample_id', 'Mutation_type'])['Pos'].agg(['count'])
    print(potential_analysis)


def calc_regression_lines(unified_freq_df):
    # calc slopes etc. for all patients
    regression_coeefs_summary = pd.DataFrame(columns=['ind_id','mut_type','slope', 'p_value'])

    # quick an dirty, but precise
    for mut_type in ["synonymous", "missense"]:
        for ind_id in ['12796','13003','15664','16207','17339','19937','22097','22763','22828','23271','26892','28545','28841','29219','29447','31254','34253','47939']:
            regressed_data = unified_freq_df[(unified_freq_df['ind_id'] == int(ind_id)) & (unified_freq_df['Mutation_type'] == mut_type)]

            regressed_data["count_for_position"] = regressed_data["Freq"] * regressed_data["Read_count"]
            regressed_data = regressed_data.groupby(['years_since_infection'])['count_for_position', 'Read_count'].agg(['sum'])
            regressed_data.columns = ['count_for_position_sum', 'Read_count_sum']
            regressed_data = regressed_data.reset_index()
            regressed_data['weighted_freq'] = regressed_data['count_for_position_sum'] / regressed_data['Read_count_sum']

            # get coeffs of linear fit
            slope, intercept, r_value, p_value, std_err = stats.linregress(regressed_data['years_since_infection'],
                                                                           regressed_data['weighted_freq'])

            regression_coeefs_summary = regression_coeefs_summary.append({'ind_id': ind_id, 'mut_type': mut_type, 'slope': "{:.2e}".format(slope), 'p_value': "{0:0.4f}".format(p_value)}, ignore_index=True)
            print(ind_id+','+mut_type+': '+ "{:.2e}".format(slope))

    # export to summary file
    regression_coeefs_summary.to_csv('/Users/omer/PycharmProjects/SternLab/RG_HIVC_analysis/output_tables/mutation_rates_regression_stats.csv')


    # regression attempt
    # g = sns.lmplot(x= "years_since_infection",
    #                y= "count_and_weight",
    #                hue="Mutation_type",
    #                hue_order=["missense", "synonymous", "stop"],
    #                col='ind_id',
    #                col_wrap=5,
    #                n_boot=100,
    #                # join=True,
    #                # kind='point',
    #                x_estimator=weighted_varaint,
    #                truncate=True,
    #                logx= True,
    #                data=unified_freq_df)


if __name__ == "__main__":
    # generate unified
    # unified_freq_df = generate_unified_filtered_verbose_freqs_df()

    # get unified
    unified_freq_df = pd.read_csv('/Users/omer/PycharmProjects/SternLab/RG_HIVC_analysis/ET86_4s/with_muts_and_con_as_ref/unified_freqs_filtered_verbose.csv')
    # unified_freq_df = pd.read_csv('/Users/omer/PycharmProjects/SternLab/RG_HIVC_analysis/ET86_2s/with_muts_coordinated/unified_freqs_with_ids_ysi.csv')

    # freqs sanity check
    # print(unified_freq_df[(unified_freq_df['Freq'] > 0.5) & (unified_freq_df['Mutation_type'] == 'stop')].head().to_string())
    print(len(unified_freq_df[(unified_freq_df['Freq'] > 0.5) & (unified_freq_df['Mutation_type'] == 'stop')]))

    # additional filters
    # G->A only
    # print(len(unified_freq_df))
    # unified_freq_df = unified_freq_df[unified_freq_df['Mutation'] == 'GA']
    # print(len(unified_freq_df))

    # plot_mutation_rate_distribution(unified_freq_df)
    # plot_mutation_rate_conf_interval(unified_freq_df)

    # manual verification analysis
    # verify_chosen_stats(unified_freq_df)
    # potential_analysis(unified_freq_df)

    # regression line extraction
    calc_regression_lines(unified_freq_df)
