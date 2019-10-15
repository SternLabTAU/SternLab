import pandas as pd
import numpy as np
import math
import matplotlib.pyplot as plt
import itertools
import sys
import argparse

from scipy.stats import chi2_contingency, fisher_exact


import seaborn as sns;

sns.set_context("poster")

# TODO- finish reviewing all different pipeline output types.
# TODO- learn blast files & algorithm (2nd iteration- understand pipeline's full flow)

# TODO- go over inclarities in flow (with maoz)

def mutations_association(args):
    # position to handle (will process interval of 250 positions ahead)
    input_x=args.position

    freqs=pd.read_csv(args.freqs_file, sep="\t")
    freqs = freqs[freqs['Pos'] == np.round(freqs['Pos'])] #remove insertions #TODO- maoz- why not simply ref != '-'
    if (input_x < freqs["Pos"].min()) or (input_x > freqs["Pos"].max()):
        sys.exit()

    # blast files (all .fasta.blast files joined together)
    all_mappings=pd.read_csv(args.blast_output, names=["read_id","start","end"], sep="\t")
    # summary of all observed mutations from ref, including mappings to origin reads
    all_mutations=pd.read_csv(args.mutations_all, names=["pos","read_id","mutant","read_positions"], sep="\t")  #TODO- what mutations are included in mutations_all? is there a threshold?

    cons = freqs[(freqs["Rank"] == 0)
                 & (freqs["Base"] != "-")]
    cons.insert(0, "pos", pd.to_numeric(cons.loc[:,"Pos"]))

    all_mutations = pd.merge(all_mutations, cons[["pos","Ref"]], on="pos") # adding Ref\Cons to all_mutations
    #remove C>A and G>T
    #all_mutations = all_mutations[~(((all_mutations["Ref"]=="C")&(all_mutations["mutant"]=="A")) | ((all_mutations["Ref"]=="G")&(all_mutations["mutant"]=="T")))]

    #variants=all_mutations["pos"].unique()

    variants_combinations=range(input_x+1,input_x+2) # x-> (x+1,x+2) instead of (x+1,x+250)

    for y in variants_combinations:
        #x=pair[0]
        x=input_x
        #y=pair[1]
        maps_for_two_pos = all_mappings[(all_mappings["start"] <= x) & (all_mappings["end"] >= y)] # reads surrounding the [x,y] interval
        merge_read_id = pd.DataFrame({"read_id": maps_for_two_pos["read_id"].unique()})
        merge_x = all_mutations[all_mutations["pos"] == x][["pos", "read_id"]]
        merged = pd.merge(merge_read_id, merge_x, on="read_id", how="left")
        merge_y = all_mutations[all_mutations["pos"] == y][["pos", "read_id"]]
        merged = pd.merge(merged, merge_y, on="read_id", how="left")

        x_label = "pos_" + str(x)
        y_label = "pos_" + str(y)
        merged[x_label] = np.where(merged["pos_x"] == x, 1, 0)
        merged[y_label] = np.where(merged["pos_y"] == y, 1, 0)
        ct = pd.crosstab(merged[x_label], merged[y_label])
        if ct.shape == (2,2):
            fisher_test = fisher_exact(ct, alternative='greater') ## TODO- review fisher's test
            print('\t'.join([str(x) for x in [x, y, fisher_test[0], fisher_test[1], ct[1][1]*1.0/(ct[0][0]+ct[0][1]+ct[1][0]+ct[1][1])]]))
        else:
            print('\t'.join([str(x) for x in [x, y, 0.0, 1.0, 0.0]])) # statistic ('odds ratio'), p-value, *shared_freq*
            # list of all connections of pair mutations-
            # 1. should be grouped into cliques with one joint freq?- and mapped to freq plot with color per clique?
            # 2. adjacent positions should be registered with joint freq, and given as input to aa_mut script?


def freq_plot(unified_freq_df):
    # advanced, chronological presentation
    g = sns.relplot(x= "Pos",
                    y= "Freq",
                    col='sample_id',
                    # col_order='years_since_infection',
                    # # hue='sample_id',
                    # col_wrap=5,
                    # join=True,
                    data=unified_freq_df)


    # plot adjustments
    g.set(yscale="log")
    plt.ylim(10**-4, 1)
    # g.set_ylabels("mutation_rate")
    # g.set_xlabels("ET first sample (Years)")
    # g.set_xticklabels(rotation=45, fontsize=11)

    # extract plot
    plt.show()
    # plt.savefig(fname= '/Users/omer/PycharmProjects/SternLab/RG_HIVC_analysis/figures/mutation_rates_ET86_4s.png')
    # g.savefig('')

def dist_plot(freq_df):
    # g = sns.jointplot(x= "Pos", y= "Freq", data=freq_df, kind='kde')

    # plot 1
    # f, ax = plt.subplots(figsize=(10, 6))
    # sns.kdeplot(freq_df.Pos, freq_df.Freq, ax=ax)
    # sns.rugplot(freq_df.Freq, vertical=True, ax=ax)

    # plot 2
    samples = ['X83354_S92','504201_S45','X100748_S73','504223_S67']
    # f, axes = plt.subplots(2, math.ceil(len(samples)/2), figsize=(6, 6))
    # for i in range(len(samples)):
    #     sns.jointplot(x='Pos', y='Freq',
    #                  data=freq_df[freq_df['sample_id'].isin([samples[i]])],
    #                  kind= 'kde',
    #                  legend=False,
    #                  ax=axes[round(i/len(samples))][i % math.ceil(len(samples)/2)])


    # plot 3
    plt.xscale('log')
    f, axes = plt.subplots(1, len(samples), figsize=(len(samples)*6, 6))
    for i in range(len(samples)):
        sns.distplot(freq_df[freq_df['sample_id'].isin([samples[i]])]['Freq'],
                     ax=axes[i])
        # sns.distplot(freq_df.Freq);

    # plt.ylim(10**-4, 1)
    # ax.set_xticks(range(0, 9000, 1000))
    plt.show()


def main_plots():
    # get freq file\s
    # freq_df = pd.read_csv('/Users/omer/PycharmProjects/SternLab/RG_HIVC_analysis/ET86_4s/TASPX119494_S74.freqs', sep='\t')
    freq_df = pd.read_csv('/Users/omer/PycharmProjects/SternLab/RG_HIVC_analysis/ET86_4s/unified_freqs_filtered_verbose.csv')

    # filters
    freq_df = freq_df[freq_df['Rank'] != 0]
    freq_df = freq_df[freq_df['Freq'] != 0]
    # freq_df = freq_df[freq_df['ind_id'].isin(['29447'])]
    # freq_df = freq_df[freq_df['sample_id'].isin(['87415_S75', 'TASPX119494_S74'])]
    # freq_df = freq_df[freq_df['sample_id'].isin(['87415_S75'])]

    # plot
    # freq_plot(freq_df)
    dist_plot(freq_df)


def main_mutations_association():
    global args
    parser = argparse.ArgumentParser()
    parser.add_argument("blast_output", type=str, help="all BLAST output for this sample")
    parser.add_argument("mutations_all", type=str, help="mutations_all.txt file (filtered from text)")
    parser.add_argument("position", type=int, help="The position to consider pairs")
    parser.add_argument("freqs_file", type=str, help="freqs file")
    args = parser.parse_args(sys.argv[1:])

    # inferring co-occurence for a single position, with all the rest
    mutations_association(args)


if __name__ == "__main__":
    # main_plots()
    main_mutations_association()