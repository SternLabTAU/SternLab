from selection_model_analysis import *
import argparse


def main(args):
    out = r'/sternadi/home/volume1/daniellem1/Entropy/data/profile'
    zikv_dengv = r'/sternadi/home/volume1/daniellem1/Entropy/data/Phylogeny/structure/degue_zika_3utr_only_id.fasta'
    df = pd.read_csv(r'/sternadi/home/volume1/daniellem1/Entropy/data/OU_model/simulations_significance_bm_k5.csv')
    families = df['family'].values

    family = families[args.index - 1]
    fasta = r'/sternadi/home/volume1/daniellem1/Entropy/data/Phylogeny/family/{}/{}.fasta'.format(family, family)
    get_joint_entropy_profile(fasta, 200, out)


if __name__ == "__main__":
    parser = argparse.ArgumentParser()
    parser.add_argument("-i", "--index", type=int, help="array index")

    args = parser.parse_args()

    main(args)