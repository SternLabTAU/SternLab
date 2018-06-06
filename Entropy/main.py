from utils import *
from entropy_selection import *
import argparse


def main(args):

    cds = r'/sternadi/home/volume1/daniellem1/Entropy/data/virushostdb.cds.fna'
    genomics = r'/sternadi/home/volume1/daniellem1/Entropy/data/virushostdb.genomic.fna'


    out = r'/sternadi/home/volume1/daniellem1/Entropy/data/'

    k = args.kmer
    coding = args.coding
    if coding == 1:
        mapping = refseq_2_cds(cds)
        cds_mapping = entropy_by_cds(mapping, out)
    else:
        basic_genome_entropy = genome_2_entropy(genomics,k, out=out, rc_joint=True)



if __name__ == "__main__":
    parser = argparse.ArgumentParser()
    parser.add_argument("-k", "--kmer", type=int,
                        help="the kmer size", required=True)
    parser.add_argument("-c", "--coding", type=int,
                        help="lthe type of file", default=0)

    args = parser.parse_args()

    main(args)