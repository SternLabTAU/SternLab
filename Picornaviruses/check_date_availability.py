import pandas as pd
import xml.etree.ElementTree as ET
from optparse import OptionParser
from Bio import Entrez



def add_collection_date(gene_hits, output):
    '''
    Reads through all gene hits and checks through entraz if they have a collection date listed
    :param gene_hits: the gene hits file
    :param output: output file
    '''
    cnt = 0
    table = pd.read_csv(gene_hits)
    table["date"] = ""
    for id_line in list(table.id):
        print(cnt)
        cnt += 1
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