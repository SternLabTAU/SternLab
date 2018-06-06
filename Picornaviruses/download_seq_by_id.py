from optparse import OptionParser
from Bio import Entrez
from Bio import SeqIO
import pandas as pd
import os
from os import path


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
        id = id_line.split("|")[3]
        output = output_prefix + id + ".fasta"
        handle = Entrez.efetch(db="nucleotide", id=id, rettype="gb", retmode="fasta")
        record = SeqIO.read(handle, "gb")
        SeqIO.write(record, output, "fasta")
        handle.close()


def unite_seqs_to_file(folder, united_path):
    '''
    After downloading all sequences to seperate files, this function unites them into a single fasta file that will be
    accepted by Prank
    :param folder: the folder that contains all downloaded seqs
    :param united_path: output file path
    '''
    with open(united_path, "w") as united_handle:
        for file in os.listdir(folder):
            try:
                seq = SeqIO.read(folder + "/" + file, "fasta")
                SeqIO.write(seq, united_handle, "fasta")
            except ValueError:
                print("more than one record in "+ file)
    united_handle.close()



if __name__ == "__main__":
    main()