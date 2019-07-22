import matplotlib.pyplot as plt
import seaborn as sns;
import pandas as pd

from pbs_runners import script_runner

sns.set_context("poster")

def run_pipeline_vs_founder(pipeline):
    for sample in pipeline:
        freqs_filename= sample.split('_')[0] + '.freqs'
        output_filename= sample + '.fasta'
        script_runner('python /sternadi/home/volume3/omer/SternLab/NGS_analysis/make_reference_from_consensus.py '
                            '-f sternadi/home/volume1/shared/data/ref_genomes/HXB2.fasta '
                            '-p /sternadi/home/volume1/shared/analysis/HIV_ravi_gupta/run2/$sample/$freqs_filename '
                            '-o /sternadi/home/volume1/shared/analysis/HIV_ravi_gupta/refs/$output_filename '
                            '-i $PBS_ARRAY_INDEX', jnum=76, alias='make_ref_from_con_{}'.format(sample), load_python=True)


def divergence_plots():

    div_rates = pd.read_csv('/Users/omer/PycharmProjects/SternLab/RG_HIVC_analysis/divergence_rates.csv', sep=',')

    g = sns.relplot(
        x='years_since_infection',
        y='divergence',
        col='ind_id',
        data=div_rates,
        col_wrap=5,
        kind='line',
        color='firebrick',
        facet_kws={'sharex': True, 'legend_out':True},
        )
    g.set(yscale="log")
    g.set_ylabels("divergence")
    g.set_xlabels("ET first sample (Years)")
    # g.set_xticklabels(rotation=45, fontsize=14)

    # extracting plot
    plt.show()
    # plt.savefig(fname= '/Users/omer/PycharmProjects/SternLab/RG_HIVC_analysis/figures/divergence_per_patient.pdf')
    # g.savefig('/Users/omer/PycharmProjects/SternLab/RG_HIVC_analysis/figures/divergence_per_patient_2.png')


if __name__ == "__main__":
    divergence_plots()
