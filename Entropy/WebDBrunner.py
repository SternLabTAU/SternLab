
import argparse
from selection_model_analysis import *
from stats_utils import *

'''
This is a runner used for generating all the data for the web DB of the entropy project

This runner will re-generate all the entropy profiles for k=5, both joint entropy and Shannon entropy

'''


def generate_entropy_profiles_per_family(index, run_type):
    """
    this function runs the entropy analysis for all families in the virus host db dataset
    :param index: an index of the run - used by the PBS queue system
    :param run_type: the type of run
    :return:
    """

    df = pd.read_csv(r'/sternadi/home/volume1/daniellem1/Entropy/data/OU_model/simulations_significance_bm_k5.csv')

    out = r'/sternadi/home/volume1/daniellem1/Entropy/WebDB'

    families = df['family'].values
    family = families[args.index - 1]

    fasta = r'/sternadi/home/volume1/daniellem1/Entropy/data/Phylogeny/family/{}/{}.fasta'.format(family, family)

    run_type = args.type

    w = 200

    if run_type == 1:   # shannon entropy
        cur_out = os.path.join(out, 'Entropy')
        get_entropy_profile(fasta, w, cur_out)

    if run_type == 2:
        cur_out = os.path.join(out, 'JointEntropy')
        get_joint_entropy_profile(fasta, w, cur_out)


def run_statistic_analysis_for_drop_sig(family, index, m=10**4, w=100):
    """
    run stretch finder on each sequence in a family data frame
    :param family: the family name
    :param index: an index to be used by the PBS system
    :param w: the size of the linear window stretch (default =100)
    :return: saves the data frame into a file
    """

    family_df = pd.read_csv(r'/sternadi/home/volume1/daniellem1/Entropy/WebDB/JointEntropy/{}_Joint_profile.csv'.format(family))
    col = 'seq_{}'.format(args.index - 1)
    seq = np.array(family_df[col].dropna())
    res = stretchFinder(seq, w, m)
    out_dir = r'/sternadi/home/volume1/daniellem1/Entropy/WebDB/DropsStatisticsJoint/{}'.format(family)

    # create a directory per family if do not exist, if so, don't touch
    if not os.path.exists(out_dir):
        os.makedirs(out_dir)

    out = r'/sternadi/home/volume1/daniellem1/Entropy/WebDB/DropsStatisticsJoint/{}/{}_seq_{}_stats.csv'.format(family, family, col)


    res.to_csv(out, index=False)


def merge_stats_by_family(family, joint=False):
    """
    the method merges all stats csv files into a single file
    :param family: the family name
    :param joint: indicator for joint entropy. if true run the analysis on joint entropy
    :return: saves a csv with all p values
    """

    # define the folder from which all stats are going to be generated
    family_folder = r'/sternadi/home/volume1/daniellem1/Entropy/WebDB/DropsStatistics/{}'.format(family)
    if joint:
        family_folder = r'/sternadi/home/volume1/daniellem1/Entropy/WebDB/DropsStatisticsJoint/{}'.format(family)

    stats = glob.glob(os.path.join(family_folder, '*stats.csv'))
    mapping = {}
    for f in stats:
        print(f)
        df = pd.read_csv(f)
        seq_id = re.findall(r'seq_\d+', f)[0]
        mapping[seq_id] = df['corrected_pvalue']

    # create a complete dataframe with all p values
    stats_by_seq = pd.DataFrame(dict([(k, pd.Series(v)) for k, v in mapping.items()]))

    output = os.path.join(family_folder, 'all_stats.csv')

    stats_by_seq.to_csv(output, index=False)
    print('Done!')
    return stats_by_seq






def main(args):

    # parse parameters
    run_type = args.type
    index = args.index
    family = args.family

    if run_type in [1,2]:   # entropy profile calculation
        generate_entropy_profiles_per_family(run_type, index)

    elif run_type == 3:   #statistical analysis
        run_statistic_analysis_for_drop_sig(family, index)

    else:
        df = pd.read_csv(r'/sternadi/home/volume1/daniellem1/Entropy/data/OU_model/simulations_significance_bm_k5.csv')
        families = df['family'].values
        family = families[args.index - 1]
        print(family)
        merge_stats_by_family(family)
        merge_stats_by_family(family, joint=True)


if __name__ == "__main__":
    parser = argparse.ArgumentParser()
    parser.add_argument("-i", "--index", type=int, help="array index", required=True)
    parser.add_argument("-t", "--type", type=int, help="indicator for the function to execute", required=True)
    parser.add_argument("-f", "--family", type=str, help="family name", required=False)
    args = parser.parse_args()

    main(args)




