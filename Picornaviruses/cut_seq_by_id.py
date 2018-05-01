import xml.etree.ElementTree as ET
from Bio import Entrez
import pandas as pd
from Bio import Entrez
from Bio import SeqIO


def find_gene_range(gene, id_file):
    '''
    Reads through seq id file and adds the start and end locations for the specified gene according to the annotation
    information available in the record.
    :param gene: gene of interest
    :param id_file: sequence id file, where start and end positions will be added
    :return:
    '''
    table = pd.read_excel(id_file)
    table["start" + str(gene)] = ""
    table["end" + str(gene)] = ""
    res_lst = []
    for seq_id in list(table.id):
        print(seq_id)
        handle = Entrez.efetch(db="nucleotide", id=seq_id, rettype="gb", retmode="fasta")
        record = SeqIO.read(handle, "genbank")
        handle.close()
        for feature in record.features:
            if feature.type == "mat_peptide":
                for description in feature.qualifiers['product']:
                    if gene in description or gene.lower() in description:
                        start = int(feature.location.start)
                        end = int(feature.location.end)
                        table.loc[table['id'] == seq_id, "start" + str(gene)] = start
                        table.loc[table['id'] == seq_id, "end" + str(gene)] = end
                try:
                    for description in feature.qualifiers['note']:
                        if gene in description or gene.lower() in description:
                            start = int(feature.location.start)
                            end = int(feature.location.end)
                            table.loc[table['id'] == seq_id, "start" + str(gene)] = start
                            table.loc[table['id'] == seq_id, "end" + str(gene)] = end
                except KeyError:
                    print(" - no note for this gene \n")
    table.to_excel(id_file)


def add_cut_seq_to_file(id_file, seqs_file):
    '''
    Writes trimmed sequence to given file, according to the start and end positions.
    :param seqs_file: output file
    '''
    table = pd.read_excel(id_file)
    f = open(seqs_file, "w")
    for line in table.iterrows():
        print(line[1].id)
        seq_id = line[1].id # + "|" + str(line[1].year)
        start = int(line[1].startVP1)
        end = int(line[1].endVP1)
        handle = Entrez.efetch(db="nucleotide", id=line[1].id, rettype="gb", retmode="fasta")
        record = SeqIO.read(handle, "genbank")
        handle.close()
        seq = record.seq[start:end]
        f.write(">")
        f.write(seq_id + "\n")
        f.write(str(seq) + "\n")
    f.close()





