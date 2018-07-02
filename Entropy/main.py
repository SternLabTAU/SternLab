from utils import *
from entropy_selection import *
import argparse


def main():

    cds = r'/sternadi/home/volume1/daniellem1/Entropy/data/virushostdb.cds.fna'
    genomics = r'/sternadi/home/volume1/daniellem1/Entropy/data/virushostdb.genomic.fna'


    main_dir = r'/sternadi/home/volume1/daniellem1/Entropy/data/Phylogeny/family'
    # entropies = r'/sternadi/home/volume1/daniellem1/Entropy/data/entropies.csv'
    entropies = r'/sternadi/home/volume1/daniellem1/Entropy/data/codon_and_joint_integrated.csv'
    selection = r'/sternadi/home/volume1/daniellem1/Entropy/data/control/entropy_selection_test.csv'
    # selection = r'/sternadi/home/volume1/daniellem1/Entropy/data/normalization/entropy_selection_test.csv'


    mapping = pd.read_csv(entropies)
    selection_stats = pd.read_csv(selection)

    out = r'/sternadi/home/volume1/daniellem1/Entropy/data/control'

    #normalize mapping by each family maximum value
    # features = [c for c in mapping.columns if c not in ['virus_name', 'family', 'refseq_id']]
    # for f in features:
    #     mapping[f] = mapping[f] / mapping.groupby('family')[f].transform(np.max)


    # run_entropy_selection_test(main_dir, mapping, out)
    # f1 = 'codon_position_3'
    f1 = 'third codon position'
    f2 = 'k5'
    # f2 = 'codon_position_2'
    test_selection_validity(selection_stats, f1, f2, out)
    

#     k = args.kmer
#     coding = args.coding
#     if coding == 1:
#         mapping = refseq_2_cds(cds)
#         cds_mapping = entropy_by_cds(mapping, out)
#     else:
#         basic_genome_entropy = genome_2_entropy(genomics,k, out=out, rc_joint=True)


# if __name__ == "__main__":
#     parser = argparse.ArgumentParser()
#     parser.add_argument("-k", "--kmer", type=int,
#                         help="the kmer size", required=True)
#     parser.add_argument("-c", "--coding", type=int,
#                         help="lthe type of file", default=0)

#     args = parser.parse_args()

#     main(args)

if __name__ == "__main__":
    main()