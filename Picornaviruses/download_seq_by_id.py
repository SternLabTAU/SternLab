from optparse import OptionParser
from Bio import Entrez
from Bio import SeqIO
import pandas as pd


def main():
    parser = OptionParser("usage: %prog[options]")
    parser.add_option("-g", "--gene_hits", dest="gene_hits", help="gene hits file")
    parser.add_option("-o", "--output", dest="output", help="common utput path. Will be completed by the id name")

    (options, args) = parser.parse_args()
    gene_hits = options.gene_hits
    output = options.output

    download_seqs(gene_hits, output)


def download_seqs(gene_hits, output_prefix):
    '''
    Download all records from the NCBI site by id. Take only id's where the collection date is specified.
    :param gene_hits: the the file of all hits of the gene
    :param output_prefix: the common path for all the output files. Will be completed by the id.
    :return: 
    '''
    table = pd.read_csv(gene_hits)
    # ignore id's that do not have a collection date
    table = table[table['date'] == "Yes"]
    for id_line in list(table.id):
        print("hi")
        id = id_line.split("|")[3]
        output = output_prefix + id + ".fasta"
        handle = Entrez.efetch(db="nucleotide", id=id, rettype="gb", retmode="fasta")
        record = SeqIO.read(handle, "gb")
        SeqIO.write(record, output, "fasta")
        handle.close()

if __name__ == "__main__":
    main()