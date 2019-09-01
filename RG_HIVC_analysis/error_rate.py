import os

import pandas as pd
import numpy as np
import glob
import matplotlib.pyplot as plt
import seaborn as sns;

from RG_HIVC_analysis.constants import excluded_samples

#sns.set_context("poster")

def create_unified_samples_to_patient_and_dsi():
    extension = 'tsv'
    all_filenames = [i for i in glob.glob('/Users/omer/PycharmProjects/SternLab/RG_HIVC_analysis/ZN_input/tables/samples_*.{}'.format(extension))]

    # combine all files in the list
    combined_csv = pd.concat([pd.read_csv(f) for f in all_filenames])
    # export to csv
    combined_csv.to_csv("samples_to_patient_and_dsi.csv", index=False, encoding='utf-8-sig')


def generate_unified_filtered_verbose_freqs_df(min_read_count = 100, freq_threshold=0.01):
    folder_path = '/Users/omer/PycharmProjects/SternLab/RG_HIVC_analysis/ET86_4s'
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
        # freq_df = freq_df[freq_df['Mutation_type'] != 'consensus']

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


def plot_error_rate_distribution(unified_freq_df):
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


def plot_error_rate_conf_interval(unified_freq_df):
    add_weighted_freq(unified_freq_df)

    def weighted_varaint(x, **kws):
        var, count = map(np.asarray, zip(*x))
        return var.sum() / count.sum()

    g = sns.factorplot(x= "sample_id",
                       y= "count_and_weight",
                       data=unified_freq_df,
                       hue="Mutation_type",
                       hue_order=["missense", "synonymous", "stop"],
                       # palette="tab20",
                       join=False,
                       orient="v",
                       estimator=weighted_varaint,
                       dodge=1)

    # TODO- optional:
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
    # g.set_ylabels("Muts rate")
    # g.set_xlabels("ET first sample (Years)")
    g.set_xticklabels(rotation=45, fontsize=11)

    # extract plot
    plt.savefig(fname= '/Users/omer/PycharmProjects/SternLab/RG_HIVC_analysis/figures/error_rate_pol_4s_by_sample2_GA_only.png')
    plt.show()
    # g.savefig('')


def add_weighted_freq(unified_freq_df):
    unified_freq_df["count_for_position"] = unified_freq_df["Freq"] * unified_freq_df["Read_count"]
    # unified_freq_df["count_for_position"] = np.where(unified_freq_df["Prob"] < 0.95, 0, unified_freq_df["no_variants"]) #TODO- relevant?
    unified_freq_df["count_and_weight"] = list(zip(unified_freq_df.count_for_position, unified_freq_df.Read_count))

    unified_freq_df.to_csv('/Users/omer/PycharmProjects/SternLab/RG_HIVC_analysis/ET86_2s/with_muts_coordinated/examine_freqs_with_weights.csv', index=False)
    print(unified_freq_df.shape)


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

    # plot_error_rate_distribution(unified_freq_df)
    plot_error_rate_conf_interval(unified_freq_df)
