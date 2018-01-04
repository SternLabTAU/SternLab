import pandas as pd
import xml.etree.ElementTree as ET
from optparse import OptionParser
from Bio import Entrez


def main():
    parser = OptionParser("usage: %prog[options]")
    parser.add_option("-g", "--gene_hits", dest="gene_hits", help="gene hits file")
    parser.add_option("-o", "--output", dest="output", default="same as input", help="Output path. When left empty rewrited input")

    (options, args) = parser.parse_args()

    gene_hits = options.gene_hits
    if options.output == "same as input":
        output = gene_hits
    else:
        output = options.output

    add_collection_date(gene_hits, output)


def add_collection_date(gene_hits, output):
    '''
    Reads through all gene hits and checks through entraz if they have a collection date listed
    :param gene_hits: the gene hits file
    :param output: output file
    '''
    table = pd.read_csv(gene_hits)
    table["date"] = ""
    for id_line in list(table.id):
        id = id_line.split("|")[3]
        handle = Entrez.efetch(db="nucleotide", id=id, rettype="gb", retmode="text")
        record = handle.read()
        handle.close()
        if "collection_date" in record:
            table.loc[table['id'] == id_line, "date"] = "Yes"
        else:
            table.loc[table['id'] == id_line, "date"] = "No"
    table.to_csv(output)


if __name__ == "__main__":
    main()