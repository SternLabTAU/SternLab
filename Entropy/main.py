from utils import *

def main():

    cds = r'/sternadi/home/volume1/daniellem1/Entropy/data/virushostdb.cds.fna'
    genomics = r'/sternadi/home/volume1/daniellem1/Entropy/data/virushostdb.genomic.fna'
    virushost_path = r'/sternadi/home/volume1/daniellem1/Entropy/data/virus_host_db.xlsx'

    out = r'/sternadi/home/volume1/daniellem1/Entropy/data/cds'
    fasta_out = r'/sternadi/home/volume1/daniellem1/Entropy/data/Phylogeny/genus'

    entropy_2_refseq(genomic=None, cds, out, 3)

    #split_fasta_by_feature(genomics, fasta_out, virushost_path, 'genus')




if __name__ == '__main__':
    main()