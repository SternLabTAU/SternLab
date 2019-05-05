from selection_model_analysis import *
import argparse
import glob
from scipy.stats import spearmanr


def main(args):


    #run profile
    # out = r'/sternadi/home/volume1/daniellem1/Entropy/data/profile'
    # zikv_dengv = r'/sternadi/home/volume1/daniellem1/Entropy/data/Phylogeny/structure/degue_zika_3utr_only_id.fasta'
    # df = pd.read_csv(r'/sternadi/home/volume1/daniellem1/Entropy/data/OU_model/simulations_significance_bm_k5.csv')
    # families = df['family'].values

    #
    # if args.index == 61:
    #     get_entropy_profile(zikv_dengv, 60, out)
    #     deltaG_profile(zikv_dengv, 60, out)
    #     get_joint_entropy_profile(zikv_dengv, 60, out)
    # else:
    #     family = families[args.index - 1]
    #     fasta = r'/sternadi/home/volume1/daniellem1/Entropy/data/Phylogeny/family/{}/{}.fasta'.format(family, family)
    #     get_joint_entropy_profile(fasta, 200, out)

    # calculate correlations of the secondary structures analysis
    ent = []
    g = []
    joint = []
    idx = args.index
    db = glob.glob(r'/sternadi/home/volume1/daniellem1/Entropy/data/Phylogeny/structure/known_structures/DBS/{}*'.format(idx))[0]


    for rec in tqdm(SeqIO.parse(db, "fasta")):
        seq = str(rec.seq).lower()
        try:
            joint_ent = joint_entropy(seq, str(get_reverse_complement(seq)), 5)
            entropy = entropy_by_kmer(seq,5)
        except:
            continue
        delta = deltaG_calculator(seq)
        ent.append(entropy)
        joint.append(joint_ent)
        g.append(delta)

    db_type = os.path.basename(db).split('_')[1]

    if '.fasta' in db_type:
        db_type = db_type.split('.fas')[0]

    df = pd.DataFrame({'entropy':ent, 'mfe':g, 'joint':joint, 'type':db_type})

    df.to_csv(r'/sternadi/home/volume1/daniellem1/Entropy/data/Phylogeny/structure/known_structures/DBS/ent_2_g_{}.csv'.format(db_type))


if __name__ == "__main__":
    parser = argparse.ArgumentParser()
    parser.add_argument("-i", "--index", type=int, help="array index")

    args = parser.parse_args()

    main(args)