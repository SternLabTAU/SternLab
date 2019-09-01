import pandas as pd
import glob
import os
import matplotlib.pyplot as plt
import seaborn as sns;


from RG_HIVC_analysis.constants import gag_ET86_start, gag_ET86_end, pol_ET86_start, pol_ET86_end, env_ET86_start, \
    env_ET86_end


def create_coverage_stats_table():
    freq_files = glob.glob('/sternadi/home/volume1/shared/analysis/HIV_ravi_gupta/run3/**/*.freqs')
    coverage_stats = pd.DataFrame(columns=['sample_id', 'coverage_median', 'long_coverage', 'gag_cov', 'pol_cov', 'env_cov', 'file_size'])

    i=0
    for file in freq_files:
        sample_id = str(file).split("/")[8]
        cov_median = extract_coverage_median(file)
        long_cov, gag_cov, pol_cov, env_cov = extract_long_cov(file)
        file_size = os.path.getsize("/sternadi/datasets/volume2/HIV_ravi_gupta_processed/"+sample_id+"/"+sample_id+"_R1.fastq.gz")
        file_size += os.path.getsize(
            "/sternadi/datasets/volume2/HIV_ravi_gupta_processed/" + sample_id + "/" + sample_id + "_R2.fastq.gz")
        file_size /= 1024

        row = [sample_id] + [cov_median] + [long_cov] + [gag_cov] + [pol_cov] + [env_cov] + [file_size]
        # print(row)
        coverage_stats.loc[i] = row
        i=i+1

    coverage_stats = coverage_stats.sort_values(by='coverage_median', ascending=False)
    return coverage_stats

def extract_coverage_median(freq_file):
    df = pd.read_csv(freq_file, sep='\t')
    df_filtered = df.loc[df['Read_count'] > 100]
    cov_median = df_filtered['Read_count'].median()
    return cov_median

def extract_long_cov(freq_file):
    df = pd.read_csv(freq_file, sep='\t')
    df = df.drop_duplicates("Pos")
    df = df.loc[df['Read_count'] > 1000]
    long_cov = df['Pos'].count()

    gag_cov = df['Pos'].loc[(df['Pos'] > gag_ET86_start) & (df['Pos'] < gag_ET86_end)].count()
    pol_cov = df['Pos'].loc[(df['Pos'] > pol_ET86_start) & (df['Pos'] < pol_ET86_end)].count()
    env_cov = df['Pos'].loc[(df['Pos'] > env_ET86_start) & (df['Pos'] < env_ET86_end)].count()
    return long_cov, gag_cov, pol_cov, env_cov


def compare_coverage_mad_increase():
    freq_files_run2 = glob.glob('/sternadi/home/volume1/shared/analysis/HIV_ravi_gupta/run2/**/*.freqs')
    medians = pd.DataFrame(columns=['sample_id', 'coverage_run2', 'coverage_run3', 'delta_percent'])

    i=0
    for file_run2 in freq_files_run2:
        sample_id = str(file_run2).split("/")[8]
        cov_mad2 = extract_cov_mad(file_run2)
        file_run3 = glob.glob(file_run2.replace('run2','run3'))[0]
        cov_mad3 = extract_cov_mad(file_run3)
        delta_percent = (cov_mad3/cov_mad2 - 1) * 100

        row = [sample_id] + [cov_mad2] + [cov_mad3] + [delta_percent]
        print(row)
        # medians.loc[i] = row
        i=i+1

    medians.sort_values(by='coverage_run3', ascending=False)

def extract_cov_mad(freq_file):
    df = pd.read_csv(freq_file, sep='\t')
    df_filtered = df.loc[df['Read_count'] > 100]
    cov_mad = df_filtered['Read_count'].mad()
    return cov_mad



def create_general_stats_table():
    summary_files = glob.glob('/sternadi/home/volume1/shared/analysis/HIV_ravi_gupta/run3/**/pipeline_summary.txt')
    general_stats = pd.DataFrame(columns=['sample_id', 'reads_mapped_to_ref', 'reads_in_freq_count', 'reads_mapped_once', 'bases_called'])

    i = 0
    for summary_file in summary_files:
        sample_id = str(summary_file).split("/")[8]
        with open(summary_file, 'r') as file:
            content = file.read().split()
            reads_mapped_to_ref = int(content[content.index('reference:')+1])
            reads_in_freq_count = int(content[content.index('count:')+1])
            reads_mapped_once = int(content[content.index('once:')+1])
            bases_called = int(content[content.index('called:')+1])


        row = [sample_id] + [reads_mapped_to_ref] + [reads_in_freq_count] + [reads_mapped_once] + [bases_called]
        # print(row)
        general_stats.loc[i] = row
        i = i + 1

    general_stats = general_stats.sort_values(by='reads_mapped_to_ref', ascending=False)
    return general_stats

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


def get_coverage_distribution_from_unified_freq_df():
    unified_freq_df_2s = pd.read_csv('/Users/omer/PycharmProjects/SternLab/RG_HIVC_analysis/ET86_2s/unified_freqs_filtered_verbose.csv')
    unified_freq_df_2s_muts = pd.read_csv('/Users/omer/PycharmProjects/SternLab/RG_HIVC_analysis/ET86_2s/with_muts_coordinated/unified_freqs_filtered_verbose.csv')
    unified_freq_df_4s = pd.read_csv('/Users/omer/PycharmProjects/SternLab/RG_HIVC_analysis/ET86_4s/unified_freqs_filtered_verbose.csv')
    unified_freq_df_4s_muts = pd.read_csv('/Users/omer/PycharmProjects/SternLab/RG_HIVC_analysis/ET86_4s/with_muts_and_con_as_ref/unified_freqs_filtered_verbose.csv')

    freq_series = [unified_freq_df_2s, unified_freq_df_2s_muts, unified_freq_df_4s, unified_freq_df_4s_muts]

    for f in freq_series:
        # coverage_median = f.groupby('ind_id').mean()['Read_count']
        # print(coverage_median)

        # g = sns.distplot(f['Read_count'])
        g = sns.boxplot(x='ind_id', y='Read_count', data=f)

        plt.show()

    # b = unified_freq_df[unified_freq_df["Read_count"] > 1000]
    # long_cov = b.groupby('ind_id').count()['Read_count']
    # print(long_cov)



if __name__ == "__main__":
    # coverage_stats = create_coverage_stats_table()
    # coverage_stats.to_csv(path_or_buf='/sternadi/home/volume1/shared/analysis/HIV_ravi_gupta/coverage_stats_ZA04_2.csv', index=False)
    get_coverage_distribution_from_unified_freq_df()
