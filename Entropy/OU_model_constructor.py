from utils import *
import glob
import seaborn as sns
import matplotlib.pyplot as plt
import statsmodels.stats.multitest as multi


def pre_process_files(super_folder, mapping, maxNorm=False):
    """
    pre process the mapping file to fit each phylogenetic tree exists under super folder
    :param super_folder: a folder contains phylogenies. can be nested in different directories under super folder.
    :param mapping: a mapping of each refseq id to a quantitative trait
    :param max_norm: normalize each separate file by max value division
    :return: saves a csv file in the same directory of the phylogenetic tree.
    """


    # get all tree files
    all_trees = []

    for root, dirs, files in tqdm(os.walk(super_folder)):
        tree = [f for f in files if 'phyml_tree' in f]

        if tree != []:
            all_trees.append(os.path.join(root, tree[0]))

    # for each tree calculate its own data frame
    for t in tqdm(all_trees):
        if tree_2_string(t) == '':
            continue
        alias = os.path.basename(t).split('.')[0].strip()
        tree = Phylo.read(t, 'newick')
        term_names = [term.name for term in tree.get_terminals()]
        splitted_names = [term.name.split('.')[0] for term in tree.get_terminals()]
        cur_mapping = mapping[mapping['refseq_id'].isin(splitted_names)]
        data = pd.DataFrame({'node_name': term_names})
        data['refseq_id'] = data['node_name'].apply(lambda x: x.split('.')[0])
        merged = pd.merge(data, cur_mapping, on='refseq_id')

        merged.drop_duplicates(['node_name'], inplace=True)
        merged = merged.dropna()
        if maxNorm:
            features = [c for c in merged.columns if c not in ['virus_name', 'family', 'refseq_id', 'node_name']]
            for f in features:
                merged[f] = merged[f] / merged[f].max()

            merged.to_csv(os.path.join(os.path.dirname(t), 'entropies_{}.csv'.format(alias)), index=False)

    print('Preprocess is done successfully')

def run_OU_model(super_folder, features):
    """
    run the OU model R script sequentially for all phylogenies on all numerical traits
    :param super_folder:
    :return:
    """

    # get all tree files
    all_trees = []

    for root, dirs, files in tqdm(os.walk(super_folder)):
        tree = [f for f in files if 'phyml_tree' in f]

        if tree != []:
            all_trees.append(os.path.join(root, tree[0]))

    for t in tqdm(all_trees):
        if tree_2_string(t) == '':
            continue
        alias = os.path.basename(t).split('.')[0].strip()
        print(alias)
        data = os.path.join(os.path.dirname(t), 'entropies_{}.csv'.format(alias))


        # run separately for each feature
        for f in tqdm(features):
            out = os.path.join(os.path.dirname(t), 'OU_summary_{}_{}.csv'.format(f,alias))
            try:
                os.system('Rscript OU_model.R -f {} -t {} -v {} -o {}'.format(data, t, f, out))
            except:
                print(alias)


def get_family(x):
    return os.path.basename(x).split('_')[-1].split('.')[0]

def get_feature(x):
    return '_'.join(os.path.basename(x).split('_')[2:-1])


def merge_OU_results(ou_output_dir, out, mapping):
    """
    merges all the ou results from the
    :param ou_output_dir: a path to a folder containing csv's from the r script run
    :param out: a path for the resulted csv to be saved in
    :return: a data frame with information for all files, united.
    """

    all_files = glob.glob(os.path.join(ou_output_dir, '*.csv'))
    dfs = []
    for f in tqdm(all_files):
        df = pd.read_csv(f)
        family = get_family(f)
        feature = get_feature(f)
        df['family'] = family
        df['feature'] = feature
        dfs.append(df)

    result = pd.concat(dfs)

    mapping = mapping[[c for c in mapping.columns if c not in ['values', 'statistics', 'significant', 'feature', 'Model']]]
    result = pd.merge(result, mapping, on='family')

    # add model columns
    g = result[result['statistics'] == 'Pvalue'][['family', 'feature', 'values']].drop_duplicates()
    g['Model'] = g['values'].apply(lambda x: 'BM' if x > 0.05 else 'OU')
    g.drop(['values'], inplace=True, axis=1)
    result = pd.merge(result, g, on=['family', 'feature'])
    cols = ['baltimore_1', 'baltimore_2', 'kingdom', 'domain']
    for c in cols:
        result[c] = result[c].apply(lambda x: x.strip())
    result = result.drop_duplicates()
    result.to_csv(out, index=False)


def fdr_correct_by_features(df, features=None):

    # start by taking all needed features
    if features != None:
        df = df[df['feature'].isin(features)]

    raw_pvalues = df[df['statistics'] == 'Pvalue'][['family', 'feature', 'values']]
    raw_pvalues['corrected_Pvalue'] = multi.fdrcorrection(raw_pvalues['values'].values)[1]
    raw_pvalues.drop(['values'], axis=1, inplace=True)
    df = pd.merge(df, raw_pvalues, on=['family', 'feature'])
    df['corrected_Model'] = df['corrected_Pvalue'].apply(lambda x: 'BM' if x > 0.05 else 'OU')

    return df


def test_simulate_bm():

    #simulations_data = r'/Volumes/STERNADILABHOME$/volume1/daniellem1/Entropy/data/OU_model/simulations'
    simulations_data = r'/Volumes/STERNADILABHOME$/volume1/daniellem1/Entropy/data/OU_model/monte_carlo_simulations'
    real_data = r'/Volumes/STERNADILABHOME$/volume1/daniellem1/Entropy/data/OU_model/BM-OU_sampled_trees'

    simulated_files = glob.glob(os.path.join(simulations_data,'*.csv'))
    real_data_files = glob.glob(os.path.join(real_data, '*k5*.csv'))

    dfs = []
    for f in tqdm(simulated_files):

        family = f.split('_')[-1].split('.')[0].strip()
        family_real_data = glob.glob(os.path.join(real_data, '*k5*{}*.csv'.format(family)))
        if family_real_data ==[]:
            print(family)
            continue

        simulation = pd.read_csv(f)
        real = pd.read_csv(family_real_data[0])
        lower_cutoff = np.percentile(simulation['lr_bm'].values, 5)
        upper_cutoff = np.percentile(simulation['lr_bm'].values, 95)

        real_value = real[real['statistics'] == 'chiSquare']['values'].values[0]

        significance_bm = lower_cutoff <= real_value <= upper_cutoff

        pvalue = real[real['statistics'] == 'Pvalue']['values'].values[0]

        lower_cutoff = np.percentile(simulation['lr_ou'].values, 5)
        upper_cutoff = np.percentile(simulation['lr_ou'].values, 95)


        significance_ou = lower_cutoff <= real_value <= upper_cutoff

        df = pd.DataFrame({'isBM': significance_bm, 'isOU': significance_ou, 'family': family, 'pvalue': pvalue}, index=[0])

        dfs.append(df)

    result = pd.concat(dfs)

    result['corrected_pvalue'] = multi.fdrcorrection(result['pvalue'])[1]
    result['model'] = result['corrected_pvalue'].apply(lambda x: 'BM' if x > 0.05 else 'OU')

    result['significance'] = result.apply(lambda row: True if (row['isBM'] == True and row['model'] == 'BM') or
                                                              (row['isBM'] == False and row['model'] == 'OU') else False, axis=1)
    return result





# call for pre-process data for each viral family
# mapping = pd.read_csv(r'/Volumes/STERNADILABHOME$/volume1/daniellem1/Entropy/data/entropies.csv')
# pre_process_files(r'/Volumes/STERNADILABHOME$/volume1/daniellem1/Entropy/data/Phylogeny/family', mapping, maxNorm=True)

# run the R script using python iteratively
# features = [c for c in mapping.columns if c not in ['virus_name', 'family', 'refseq_id', 'node_name']]
# run_OU_model(r'/Volumes/STERNADILABHOME$/volume1/daniellem1/Entropy/data/Phylogeny/family', features)

####### Auxiliary functions for plotting ###########



def replace(group):
    mean, std = group.mean(), group.std()
    outliers = (group - mean).abs() > 3*std
    group[outliers] = mean        # or "group[~outliers].mean()"
    return group

def plot_alphas(df, features, hue=None, out=None, correct=False):


    if correct: # correct df for multiple testing
        df = fdr_correct_by_features(df, features)
        df['Model'] = df['corrected_Model']

    only_alphas = df[(df['statistics'] == 'alpha')]
    only_alphas = only_alphas[only_alphas['Model'] == 'OU']
    only_alphas = only_alphas[~(np.abs(only_alphas['values'] - only_alphas['values'].mean()) >
                                (3 * only_alphas['values'].std()))]
    only_alphas['values'] = only_alphas.groupby('feature')['values'].transform(replace)

    #only_alphas = only_alphas[only_alphas['values'] <= 20]


    only_alphas.reset_index(inplace=True)
    if hue != None:
        only_alphas = only_alphas[only_alphas[hue].isin(only_alphas[hue].value_counts()
                                                        [only_alphas[hue].value_counts() > 1].index)]
        g = sns.FacetGrid(only_alphas[only_alphas['feature'].isin(features)],
                          row="feature", hue=hue, height=2, aspect=5, palette='Dark2')
        g.map(sns.kdeplot, "values", shade=True)

    else:
        sns.distplot(only_alphas[only_alphas['feature'].isin(features)]['values'],
                                 hist= True, rug=True, kde=True, color='#C08113')
        plt.title("Alpha values distribution", fontsize=20)

    plt.legend()
    plt.xlabel('Alpha values', fontsize=20)
    sns.despine()


    if out != None:
        plt.savefig(out, dpi=400, bbox_inches='tight')

    plt.gcf().clear()





def plot_sigmas(df, features, hue=None, out=None, correct=False):

    if correct: # correct df for multiple testing
        df = fdr_correct_by_features(df, features)
        df['Model'] = df['corrected_Model']

    only_sigmas = df[((df['statistics'] == 'OUsigma') & (df['Model'] == 'OU')) |
                     ((df['statistics'] == 'BMsigma') & (df['Model'] == 'BM'))]
    only_sigmas = only_sigmas[~(np.abs(only_sigmas['values'] - only_sigmas['values'].mean()) >
                                (3 * only_sigmas['values'].std()))] # filter outliers
    only_sigmas['values'] = only_sigmas.groupby('feature')['values'].transform(replace)
    # only_sigmas = only_sigmas[only_sigmas['values'] <= 1]
    only_sigmas = only_sigmas[only_sigmas['feature'].isin(features)]
    only_sigmas.reset_index(inplace=True)

    if hue != None:
       only_sigmas = only_sigmas[only_sigmas[hue].isin(only_sigmas[hue].value_counts()
                                                       [only_sigmas[hue].value_counts() > 1].index)]

       g = sns.FacetGrid(only_sigmas,row="feature", hue=hue, height=2, aspect=5, palette='Dark2')
       g.map(sns.kdeplot, "values", shade=True)

    else:
        g = sns.distplot(only_sigmas['values'], hist=True, rug=True, kde=True, color='#C08113')

        plt.title('Sigma values distribution', fontsize=20)
    plt.legend()
    plt.xlabel('Sigma values', fontsize=20)
    sns.despine()


    if out != None:
        plt.savefig(out, dpi=400, bbox_inches='tight')

    plt.gcf().clear()


def plot_model_by_feature(df, statistic,v_features=None, hue=None, model='BM', lst_features=None, out=None, correct=False):

    if correct: # correct df for multiple testing
        if v_features != None:
            df = fdr_correct_by_features(df, v_features)
        if lst_features != None:
            df = fdr_correct_by_features(df, lst_features)
        df['Model'] = df['corrected_Model']

    df = df[(df['statistics'] == statistic) & (df['Model'] == model)]
    if lst_features != None and hue != None:
        df = df[df[hue].isin(lst_features)]
    if v_features != None:
        df = df[df['feature'].isin(v_features)]

    #df = df[~(np.abs(df['values'] - df['values'].mean()) > (3 * df['values'].std()))]
    # keep only the ones that are within +3 to -3 standard deviations in the column 'Data'.

    #df['values'] = df.groupby('feature')['values'].transform(replace)
    df = df[df[hue].isin(df[hue].value_counts()[df[hue].value_counts() > 1].index)] # remove outliers
    grid = sns.FacetGrid(df, palette="Dark2", hue=hue, height=8)
    grid.map(sns.distplot, 'values', rug=False, hist=False, kde=True, kde_kws={'shade':True})
    plt.legend(loc='best')
    plt.title('{} model {} values by {}'.format(model, statistic, hue), fontsize=22)
    plt.xlabel('{} values'.format(statistic), fontsize=20)
    plt.ylabel('Density', fontsize=20)
    plt.tight_layout()

    if out != None:
        plt.savefig(out, dpi=400, bbox_inches='tight')

    plt.gcf().clear()



def boxplot_by_model(df, statistic, feature, subset=None, hue=None, out=None, correct=True):

    if correct: # correct df for multiple testing
        if subset == None and feature == 'feature':
            subset = list(set(df['feature']))
        df = fdr_correct_by_features(df, subset)
        df['Model'] = df['corrected_Model']

    # filter data
    if 'sigma' in statistic:
        df = df[((df['statistics'] == 'OUsigma') & (df['Model'] == 'OU')) |
                     ((df['statistics'] == 'BMsigma') & (df['Model'] == 'BM'))]
    if 'alpha' in statistic:
        df = df[(df['statistics'] == 'alpha')]

    # remove outliers by feature??
    #df['values'] = df.groupby(feature)['values'].transform(replace)

    # filter any if needed
    if subset != None:
        df = df[df[feature].isin(subset)]

    sns.boxplot(x='Model', y = 'values', hue=hue, data=df, palette='Accent')

    plt.title("{} by Model for {} distribution".format(feature, statistic), fontsize=20)
    plt.xlabel("Model", fontsize=20)
    plt.ylabel("{} values".format(statistic), fontsize=20)
    plt.yscale('log')
    sns.despine()

    if out != None:
        plt.savefig(out, dpi=400, bbox_inches='tight')

    plt.gcf().clear()

def boxplot_valus_by_x(df, statistic, feature, subset=None, hue='Model', out=None):

    if 'sigma' in statistic:
        df = df[((df['statistics'] == 'OUsigma') & (df['Model'] == 'OU')) |
                ((df['statistics'] == 'BMsigma') & (df['Model'] == 'BM'))]
    if 'alpha' in statistic:
        df = df[(df['statistics'] == 'alpha')]

    # remove outliers by feature??
    # df['values'] = df.groupby(feature)['values'].transform(replace)

    # filter any if needed
    if subset != None:
        df = df[df['feature'].isin(subset)]  # do not consider all k's

    sns.boxplot(x=feature, y='values', hue=hue, data=df, palette='Dark2')

    plt.title("{} by Model for {} distribution".format(feature, statistic), fontsize=20)
    plt.xlabel(feature, fontsize=20)
    plt.ylabel("{} values".format(statistic), fontsize=20)
    plt.yscale('log')
    sns.despine()


    if out != None:
        plt.savefig(out, dpi=400, bbox_inches='tight')

    plt.gcf().clear()

def multiple_boxplots(df, statistic, x, subset, out=None):

    if 'sigma' in statistic:
        df = df[((df['statistics'] == 'OUsigma') & (df['Model'] == 'OU')) |
                ((df['statistics'] == 'BMsigma') & (df['Model'] == 'BM'))]
    if 'alpha' in statistic:
        df = df[(df['statistics'] == 'alpha')]

    df = df[df['feature'].isin(subset)]  # do not consider all k's

    # remove outliers by feature??
    # df['values'] = df.groupby(feature)['values'].transform(replace)

    fig, axes = plt.subplots(len(subset), 1, sharex=True, sharey=True)
    for i, f in enumerate(subset):
        cur_df = df[df['feature']==f]
        sns.boxplot(x=x, y='values', hue='Model', data=cur_df, palette='Dark2', ax=axes[i])
        axes[i].set_title(subset[i])
        if i != len(subset) - 1:
            axes[i].axes.get_xaxis().set_visible(False)

    plt.yscale('log')
    plt.suptitle('{} by model for {} attribute'.format(statistic,x), fontsize=20)
    plt.xlabel(x, fontsize=20)
    sns.despine()

    if out != None:
        plt.savefig(out, dpi=400, bbox_inches='tight')

    plt.gcf().clear()

def plot_pie_chart_by_feature(df, feature, column, out=None):


    colors = ['gold', 'yellowgreen', 'lightcoral', 'lightskyblue', '#68635A', '#844196', 'darkorange']

    df = df[df['feature'] == feature].drop_duplicates('family')
    grouped = df.groupby(['Model', column])['family'].count().reset_index()
    labels_2_color = {}
    print(len(set(grouped[column])))
    for i, label in enumerate(set(grouped[column])):
        labels_2_color[label] = colors[i]
    for mdl in set(grouped['Model']):
        mdl_df = grouped[grouped['Model'] == mdl]
        labels = mdl_df[column]
        adapted_colors = [labels_2_color[l] for l in labels]
        counts = mdl_df['family']
        explode = [0.1] * (len(counts) -1) + [0]
        patches, texts, autotexts = plt.pie(counts, labels=labels, colors=adapted_colors, startangle=140, autopct='%1.1f%%')
        plt.axis('equal')
        plt.title("{} Model by {} for {}".format(mdl, column, feature), fontsize=20)
        for text in texts:
            text.set_fontsize(12)
        for patch in patches:
            patch.set_edgecolor('k')
        if out != None:
            plt.savefig(os.path.join(out, '{}_piechart_{}_{}'.format(mdl,column,feature)), dpi=400, bbox_inches='tight')
        plt.gcf().clear()


def plot_pie_chart_by_model(df, feature, column, out=None):


    colors = ['#D35C37', '#97B8C2']

    for c in set(df[column]):
        cur_df = df[(df['feature'] == feature) & (df[column] == c)].drop_duplicates('family')
        counts_df = cur_df['Model'].value_counts().reset_index(name='counts').sort_values('index')
        labels = counts_df['index']
        counts = counts_df['counts']
        patches, texts, autotexts = plt.pie(counts, labels=labels, colors=colors, startangle=140,
                                            autopct='%1.1f%%')
        plt.axis('equal')
        plt.title("BM\OU Model by {}, {} for {}".format(c, column, feature), fontsize=20)
        for text in texts:
            text.set_fontsize(12)
        for patch in patches:
            patch.set_edgecolor('k')
        if out != None:
            plt.savefig(os.path.join(out, 'Model_by_{}_{}_{}'.format(c,column,feature)), dpi=400,
                        bbox_inches='tight')
        plt.gcf().clear()

def plot_minus_log_pval(df, hue, feature, out=None, correct=False):

    df = df[(df['statistics'] == 'Pvalue') & (df['feature'] == feature)].reset_index()
    if correct: # correct df for multiple testing
        df = fdr_correct_by_features(df, [feature])
        df['Model'] = df['corrected_Model']

    df['minus_log_pval'] = df['values'].apply(lambda x: -np.log(x+0.000001))
    df['idx'] = df.index

    # plot
    sns.lmplot(x='idx', y='minus_log_pval', data=df, fit_reg=False, palette='cubehelix', hue=hue, legend=False)
    plt.axhline(y=-np.log(0.05), color='#2F4E89', linestyle='-')
    plt.xticks(list(df.index), list(df['family']), rotation=90)
    plt.title('OU\BM Pvalues by  {}, {}'.format(hue, feature), fontsize=20)
    plt.ylabel('-log(P-value)', fontsize=20)
    plt.xlabel('Family', fontsize=20)
    plt.legend(bbox_to_anchor=(1, 0.65))
    plt.tight_layout()

    if out != None:
        plt.savefig(os.path.join(out, 'minus_log_pval_{}_{}'.format(feature, hue)), dpi=400,
                    bbox_inches='tight')
    plt.gcf().clear()



def generate_all_plots(df, out):
    """
    generates all OU\BM plots and saves them under the folder out
    :param df: a data frame with information about models outputs
    :param out: the folder in which the results will be saved
    :return: saves all plots
    """

    correct = True
    ks = ['k{}'.format(i) for i in range(1,6)]
    ctl = ['codon_position_1', 'codon_position_2', 'codon_position_3']
    ks_n_ctl = ['k5', 'reading_frame', 'codon_position_3']

    all_features = ['baltimore_1', 'baltimore_2', 'domain', 'kingdom', 'feature',
                    'statistics', 'Model']
    # alphas
    plot_alphas(df, ks, hue='Model', out=os.path.join(out, 'ks_alpha_model.png'), correct=correct)
    plot_alphas(df, ks, hue='baltimore_1', out=os.path.join(out, 'ks_alpha_baltimore1.png'))
    plot_alphas(df, ks, hue='baltimore_2', out=os.path.join(out, 'ks_alpha_baltimore2.png'))
    plot_alphas(df, ks, hue='kingdom', out=os.path.join(out, 'ks_alpha_kingdom.png'))
    plot_alphas(df, ks, hue='domain', out=os.path.join(out, 'ks_alpha_domain.png'))

    plot_alphas(df, ctl, hue='Model', out=os.path.join(out, 'ctl_alpha_model.png'), correct=correct)
    plot_alphas(df, ctl, hue='baltimore_1', out=os.path.join(out, 'ctl_alpha_baltimore1.png'))
    plot_alphas(df, ctl, hue='baltimore_2', out=os.path.join(out, 'ctl_alpha_baltimore2.png'))
    plot_alphas(df, ctl, hue='kingdom', out=os.path.join(out, 'ctl_alpha_kingdom.png'))
    plot_alphas(df, ctl, hue='domain', out=os.path.join(out, 'ctl_alpha_domain.png'))

    plot_alphas(df, ks_n_ctl, hue='Model', out=os.path.join(out, 'ks_n_ctl_alpha_model.png'), correct=correct)
    plot_alphas(df, ks_n_ctl, hue='baltimore_1', out=os.path.join(out, 'ks_n_ctl_alpha_baltimore1.png'))
    plot_alphas(df, ks_n_ctl, hue='baltimore_2', out=os.path.join(out, 'ks_n_ctl_alpha_baltimore2.png'))
    plot_alphas(df, ks_n_ctl, hue='kingdom', out=os.path.join(out, 'ks_n_ctl_alpha_kingdom.png'))
    plot_alphas(df, ks_n_ctl, hue='domain', out=os.path.join(out, 'ks_n_ctl__alpha_domain.png'))

    print('Done with kde alphas!\n')
    # sigmas
    plot_sigmas(df, ks, hue='Model', out=os.path.join(out, 'ks_sigma_model.png'), correct=correct)
    plot_sigmas(df, ks, hue='baltimore_1', out=os.path.join(out, 'ks_sigma_baltimore1.png'))
    plot_sigmas(df, ks, hue='baltimore_2', out=os.path.join(out, 'ks_sigma_baltimore2.png'))

    plot_sigmas(df, ctl, hue='Model', out=os.path.join(out, 'ctl_sigma_model.png'))
    plot_sigmas(df, ctl, hue='baltimore_1', out=os.path.join(out, 'ctl_sigma_baltimore1.png'))
    plot_sigmas(df, ctl, hue='baltimore_2', out=os.path.join(out, 'ctl_sigma_baltimore2.png'))

    plot_sigmas(df, ks_n_ctl, hue='Model', out=os.path.join(out, 'ks_n_ctl_sigma_model.png'), correct=correct)
    plot_sigmas(df, ks_n_ctl, hue='baltimore_1', out=os.path.join(out, 'ks_n_ctl_sigma_baltimore1.png'))
    plot_sigmas(df, ks_n_ctl, hue='baltimore_2', out=os.path.join(out, 'ks_n_ctl_sigma_baltimore2.png'))

    print('Done with kde sigma!\n')
    #model by feature
    plot_model_by_feature(df, 'alpha', hue='feature', model='OU', lst_features=ks,
                          out=os.path.join(out, 'feature_by_OU_alpha_ks.png'), correct=correct)
    plot_model_by_feature(df, 'alpha', hue='feature', model='OU', lst_features=ctl,
                          out=os.path.join(out, 'feature_by_OU_alpha_ctl.png'),correct=correct)
    plot_model_by_feature(df, 'alpha', hue='feature', model='OU', lst_features=ks_n_ctl,
                          out=os.path.join(out, 'feature_by_OU_alpha_ks_n_ctl.png'),correct=correct)

    plot_model_by_feature(df, 'OUsigma', hue='feature', model='OU', lst_features=ks,
                          out=os.path.join(out, 'feature_by_OU_sigma_ks.png'),correct=correct)
    plot_model_by_feature(df, 'OUsigma', hue='feature', model='OU', lst_features=ctl,
                          out=os.path.join(out, 'feature_by_OU_sigma_ctl.png'),correct=correct)
    plot_model_by_feature(df, 'OUsigma', hue='feature', model='OU', lst_features=ks_n_ctl,
                          out=os.path.join(out, 'feature_by_OU_sigma_ks_n_ctl.png'),correct=correct)

    plot_model_by_feature(df, 'BMsigma', hue='feature', model='BM', lst_features=ks,
                          out=os.path.join(out, 'feature_by_BM_sigma_ks.png'),correct=correct)
    plot_model_by_feature(df, 'BMsigma', hue='feature', model='BM', lst_features=ctl,
                          out=os.path.join(out, 'feature_by_BM_sigma_ctl.png'),correct=correct)
    plot_model_by_feature(df, 'BMsigma', hue='feature', model='BM', lst_features=ks_n_ctl,
                          out=os.path.join(out, 'feature_by_BM_sigma_ks_n_ctl.png'),correct=correct)
    #model by baltimore_1
    plot_model_by_feature(df, 'alpha', hue='baltimore_1', model='OU',v_features=['k5'],
                          out=os.path.join(out, 'baltimore_1_by_OU_alpha_k5.png'),correct=correct)
    plot_model_by_feature(df, 'alpha', hue='baltimore_1', model='OU', v_features=['k5'],
                          out=os.path.join(out, 'baltimore_1_by_OU_alpha_k5.png'),correct=correct)
    plot_model_by_feature(df, 'alpha', hue='baltimore_1', model='OU', v_features=['k5'],
                          out=os.path.join(out, 'baltimore_1_by_OU_alpha_k5.png'),correct=correct)

    plot_model_by_feature(df, 'OUsigma', hue='baltimore_1', model='OU', v_features=['k5'],
                          out=os.path.join(out, 'baltimore_1_by_OU_sigma_k5.png'),correct=correct)
    plot_model_by_feature(df, 'OUsigma', hue='baltimore_1', model='OU', v_features=['k5'],
                          out=os.path.join(out, 'baltimore_1_by_OU_sigma_k5.png'),correct=correct)
    plot_model_by_feature(df, 'OUsigma', hue='baltimore_1', model='OU', v_features=['k5'],
                          out=os.path.join(out, 'baltimore_1_by_OU_sigma_k5.png'),correct=correct)

    plot_model_by_feature(df, 'BMsigma', hue='baltimore_1', model='BM', v_features=['k5'],
                          out=os.path.join(out, 'baltimore_1_by_BM_sigma_k5.png'),correct=correct)
    plot_model_by_feature(df, 'BMsigma', hue='baltimore_1', model='BM', v_features=['k5'],
                          out=os.path.join(out, 'baltimore_1_by_BM_sigma_k5.png'),correct=correct)
    plot_model_by_feature(df, 'BMsigma', hue='baltimore_1', model='BM', v_features=['k5'],
                          out=os.path.join(out, 'baltimore_1_by_BM_sigma_k5.png'),correct=correct)

    # model by baltimore_2
    plot_model_by_feature(df, 'alpha', hue='baltimore_2', model='OU',v_features=['k5'],
                          out=os.path.join(out, 'baltimore_2_by_OU_alpha_k5.png'),correct=correct)
    plot_model_by_feature(df, 'alpha', hue='baltimore_2', model='OU', lst_features=ctl,
                          out=os.path.join(out, 'baltimore_2_by_OU_alpha_k5.png'),correct=correct)
    plot_model_by_feature(df, 'alpha', hue='baltimore_2', model='OU', v_features=['k5'],
                          out=os.path.join(out, 'baltimore_2_by_OU_alpha_k5.png'),correct=correct)

    plot_model_by_feature(df, 'OUsigma', hue='baltimore_2', model='OU', v_features=['k5'],
                          out=os.path.join(out, 'baltimore_2_by_OU_sigma_k5.png'),correct=correct)
    plot_model_by_feature(df, 'OUsigma', hue='baltimore_2', model='OU', v_features=['k5'],
                          out=os.path.join(out, 'baltimore_2_by_OU_sigma_k5.png'),correct=correct)
    plot_model_by_feature(df, 'OUsigma', hue='baltimore_2', model='OU', v_features=['k5'],
                          out=os.path.join(out, 'baltimore_2_by_OU_sigma_k5.png'),correct=correct)

    plot_model_by_feature(df, 'BMsigma', hue='baltimore_2', model='BM', v_features=['k5'],
                          out=os.path.join(out, 'baltimore_2_by_BM_sigma_k5.png'),correct=correct)
    plot_model_by_feature(df, 'BMsigma', hue='baltimore_2', model='BM', v_features=['k5'],
                          out=os.path.join(out, 'baltimore_2_by_BM_sigma_k5.png'),correct=correct)
    plot_model_by_feature(df, 'BMsigma', hue='baltimore_2', model='BM', v_features=['k5'],
                          out=os.path.join(out, 'baltimore_2_by_BM_sigma_k5.png'),correct=correct)

    # model by kingdom
    plot_model_by_feature(df, 'alpha', hue='kingdom', model='OU', v_features=['k5'],
                          out=os.path.join(out, 'kingdom_by_OU_alpha_k5.png'),correct=correct)
    plot_model_by_feature(df, 'alpha', hue='kingdom', model='OU', v_features=['k5'],
                          out=os.path.join(out, 'kingdom_by_OU_alpha_k5.png'),correct=correct)
    plot_model_by_feature(df, 'alpha', hue='kingdom', model='OU', v_features=['k5'],
                          out=os.path.join(out, 'kingdom_by_OU_alpha_k5.png'),correct=correct)

    plot_model_by_feature(df, 'OUsigma', hue='kingdom', model='OU', v_features=['k5'],
                          out=os.path.join(out, 'kingdom_by_OU_sigma_k5.png'),correct=correct)
    plot_model_by_feature(df, 'OUsigma', hue='kingdom', model='OU', v_features=['k5'],
                          out=os.path.join(out, 'kingdom_by_OU_sigma_k5.png'),correct=correct)
    plot_model_by_feature(df, 'OUsigma', hue='kingdom', model='OU', v_features=['k5'],
                          out=os.path.join(out, 'kingdom_by_OU_sigma_k5.png'),correct=correct)

    plot_model_by_feature(df, 'BMsigma', hue='kingdom', model='BM', v_features=['k5'],
                          out=os.path.join(out, 'kingdom_by_BM_sigma_k5.png'),correct=correct)
    plot_model_by_feature(df, 'BMsigma', hue='kingdom', model='BM', v_features=['k5'],
                          out=os.path.join(out, 'kingdom_by_BM_sigma_k5.png'),correct=correct)
    plot_model_by_feature(df, 'BMsigma', hue='kingdom', model='BM', v_features=['k5'],
                          out=os.path.join(out, 'kingdom_by_BM_sigma_k5.png'),correct=correct)

    # model by domain
    plot_model_by_feature(df, 'alpha', hue='domain', model='OU', v_features=['k5'],
                          out=os.path.join(out, 'domain_by_OU_alpha_k5.png'),correct=correct)
    plot_model_by_feature(df, 'alpha', hue='domain', model='OU', v_features=['k5'],
                          out=os.path.join(out, 'domain_by_OU_alpha_k5.png'),correct=correct)
    plot_model_by_feature(df, 'alpha', hue='domain', model='OU', v_features=['k5'],
                          out=os.path.join(out, 'domain_by_OU_alpha_k5.png'),correct=correct)

    plot_model_by_feature(df, 'OUsigma', hue='domain', model='OU', v_features=['k5'],
                          out=os.path.join(out, 'domain_by_OU_sigma_k5.png'),correct=correct)
    plot_model_by_feature(df, 'OUsigma', hue='domain', model='OU', v_features=['k5'],
                          out=os.path.join(out, 'domain_by_OU_sigma_k5.png'),correct=correct)
    plot_model_by_feature(df, 'OUsigma', hue='domain', model='OU', v_features=['k5'],
                          out=os.path.join(out, 'domain_by_OU_sigma_k5.png'),correct=correct)

    plot_model_by_feature(df, 'BMsigma', hue='domain', model='BM', v_features=['k5'],
                          out=os.path.join(out, 'domain_by_BM_sigma_k5.png'),correct=correct)
    plot_model_by_feature(df, 'BMsigma', hue='domain', model='BM', v_features=['k5'],
                          out=os.path.join(out, 'domain_by_BM_sigma_k5.png'),correct=correct)
    plot_model_by_feature(df, 'BMsigma', hue='domain', model='BM', v_features=['k5'],
                          out=os.path.join(out, 'domain_by_BM_sigma_k5.png'),correct=correct)

    print('Done with model by all features!\n')
    # # plot boxplots
    boxplot_by_model(df, 'alpha', 'feature', hue='feature', out=os.path.join(out, 'alpha_by_model_n_feature.png'),correct=correct)
    boxplot_by_model(df, 'sigma', 'feature', hue='feature', out=os.path.join(out, 'sigma_by_model_n_feature.png'),correct=correct)
    #
    # boxplot_valus_by_x(df, 'alpha', feature='baltimore_1', out=os.path.join(out, 'alpha_x_baltimore_1.png'))
    # boxplot_valus_by_x(df, 'alpha', feature='baltimore_2', out=os.path.join(out, 'alpha_x_baltimore_2.png'))
    # boxplot_valus_by_x(df, 'alpha', feature='kingdom', out=os.path.join(out, 'alpha_x_kingdom.png'))
    # boxplot_valus_by_x(df, 'alpha', feature='domain', out=os.path.join(out, 'alpha_x_domain.png'))
    #
    # boxplot_valus_by_x(df, 'sigma', feature='baltimore_1', out=os.path.join(out, 'sigma_x_baltimore_1.png'))
    # boxplot_valus_by_x(df, 'sigma', feature='baltimore_2', out=os.path.join(out, 'sigma_x_baltimore_2.png'))
    # boxplot_valus_by_x(df, 'sigma', feature='kingdom', out=os.path.join(out, 'sigma_x_kingdom.png'))
    # boxplot_valus_by_x(df, 'sigma', feature='domain', out=os.path.join(out, 'sigma_x_domain.png'))
    #
    # print('Done with boxplots!\n')
    # # plot multi boxplots alpha
    # multiple_boxplots(df, 'alpha', x='baltimore_1', subset=ks, out=os.path.join(out, 'alpha_multi_baltimore_1_ks.png'))
    # multiple_boxplots(df, 'alpha', x='baltimore_2', subset=ks, out=os.path.join(out, 'alpha_multi_baltimore_2_ks.png'))
    # multiple_boxplots(df, 'alpha', x='kingdom', subset=ks, out=os.path.join(out, 'alpha_multi_kingdom_ks.png'))
    # multiple_boxplots(df, 'alpha', x='domain', subset=ks, out=os.path.join(out, 'alpha_multi_domain_1_ks.png'))
    #
    # multiple_boxplots(df, 'alpha', x='baltimore_1', subset=ctl,
    #                   out=os.path.join(out, 'alpha_multi_baltimore_1_ctl.png'))
    # multiple_boxplots(df, 'alpha', x='baltimore_2', subset=ctl,
    #                   out=os.path.join(out, 'alpha_multi_baltimore_2_ctl.png'))
    # multiple_boxplots(df, 'alpha', x='kingdom', subset=ctl, out=os.path.join(out, 'alpha_multi_kingdom_ctl.png'))
    # multiple_boxplots(df, 'alpha', x='domain', subset=ctl, out=os.path.join(out, 'alpha_multi_domain_1_ctl.png'))
    #
    # multiple_boxplots(df, 'alpha', x='baltimore_1', subset=ks_n_ctl,
    #                   out=os.path.join(out, 'alpha_multi_baltimore_1_ks_n_ctl.png'))
    # multiple_boxplots(df, 'alpha', x='baltimore_2', subset=ks_n_ctl,
    #                   out=os.path.join(out, 'alpha_multi_baltimore_2_ks_n_ctl.png'))
    # multiple_boxplots(df, 'alpha', x='kingdom', subset=ks_n_ctl,
    #                   out=os.path.join(out, 'alpha_multi_kingdom_ks_n_ctl.png'))
    # multiple_boxplots(df, 'alpha', x='domain', subset=ks_n_ctl,
    #                   out=os.path.join(out, 'alpha_multi_domain_1_ks_n_ctl.png'))
    #
    # # multi boxplots sigma
    # multiple_boxplots(df, 'sigma', x='baltimore_1', subset=ks, out=os.path.join(out, 'sigma_multi_baltimore_1_ks.png'))
    # multiple_boxplots(df, 'sigma', x='baltimore_2', subset=ks, out=os.path.join(out, 'sigma_multi_baltimore_2_ks.png'))
    # multiple_boxplots(df, 'sigma', x='kingdom', subset=ks, out=os.path.join(out, 'sigma_multi_kingdom_ks.png'))
    # multiple_boxplots(df, 'sigma', x='domain', subset=ks, out=os.path.join(out, 'sigma_multi_domain_1_ks.png'))
    #
    # multiple_boxplots(df, 'sigma', x='baltimore_1', subset=ctl,
    #                   out=os.path.join(out, 'sigma_multi_baltimore_1_ctl.png'))
    # multiple_boxplots(df, 'sigma', x='baltimore_2', subset=ctl,
    #                   out=os.path.join(out, 'sigma_multi_baltimore_2_ctl.png'))
    # multiple_boxplots(df, 'sigma', x='kingdom', subset=ctl, out=os.path.join(out, 'sigma_multi_kingdom_ctl.png'))
    # multiple_boxplots(df, 'sigma', x='domain', subset=ctl, out=os.path.join(out, 'sigma_multi_domain_1_ctl.png'))
    #
    # multiple_boxplots(df, 'sigma', x='baltimore_1', subset=ks_n_ctl,
    #                   out=os.path.join(out, 'sigma_multi_baltimore_1_ks_n_ctl.png'))
    # multiple_boxplots(df, 'sigma', x='baltimore_2', subset=ks_n_ctl,
    #                   out=os.path.join(out, 'sigma_multi_baltimore_2_ks_n_ctl.png'))
    # multiple_boxplots(df, 'sigma', x='kingdom', subset=ks_n_ctl,
    #                   out=os.path.join(out, 'sigma_multi_kingdom_ks_n_ctl.png'))
    # multiple_boxplots(df, 'sigma', x='domain', subset=ks_n_ctl,
    #                   out=os.path.join(out, 'sigma_multi_domain_1_ks_n_ctl.png'))


    print('Done with multi boxplots!\n')
    print('Done!!!!!!!')





ou_results = r'/Volumes/STERNADILABHOME$/volume1/daniellem1/Entropy/data/OU_model/ou_results_unite_sampled.csv'

sns.set_style('white')
df = pd.read_csv(ou_results)
df = df[df['num_sequences_in_tree'] >= 10]
df = df[df['statistics'] != 'chiSquare']
df['values'] = df['values'].astype(float)

out = r'/Users/daniellemiller/Google Drive/Msc Bioinformatics/Projects/entropy/most_updated/OU_BM/plots_sampling_corrected'
# generate_all_plots(df, out)


wanted_cols = ['baltimore_1', 'baltimore_2', 'kingdom', 'domain']
wanted_rows = [c for c in set(df['feature']) if 'shift' not in c]
# #
# for r in tqdm(wanted_rows):
#     for c in tqdm(wanted_cols):
# #         # plot_pie_chart_by_feature(df, feature=r, column=c, out=out)
# #         plot_pie_chart_by_model(df, feature=r, column=c, out=out)
#         plot_minus_log_pval(df, hue=c, feature=r, out=out)


x = test_simulate_bm()
x['feature'] = 'k5'
x.to_csv(r'/Volumes/STERNADILABHOME$/volume1/daniellem1/Entropy/data/OU_model/simulations_significance_bm_k5.csv', index=False)



# sns.distplot(sim['chi_sqr'], color='#2FA3BD', hist=False, kde_kws={'shade':True)
# plt.axvline(x=lower_cutoff, color = '#DF2115')
# plt.axvline(x=upper_cutoff, color = '#DF2115')
# plt.axvline(x=real_value, color = 'olive')
# sns.despine()
# plt.title('Picornaviridae BM simulations distribution', fontsize=22)
# plt.xlabel('Chi Square statistic', fontsize=20)
# plt.savefig(r'/Volumes/STERNADILABHOME$/volume1/daniellem1/Entropy/data/OU_model/picorna_k5_bm_sim.png', dpi=400, bbox_inches='tight')
