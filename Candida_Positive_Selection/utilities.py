import pandas as pd

FASTA_FILE_PATH = "C:\\Users\\Ella\\Google Drive\\Lab\\orf_genomic_assembly_21.fasta"
GENE_TABLE_PATH = "C:\\Users\\Ella\\Google Drive\\Lab\\gene_table.csv"
PATIENT1 = "C:\\Users\\Ella\\Google Drive\\Lab\\patient1.csv"


def create_chromosomal_map():
    """
    :return: A map that maps each chromosome to the indexes where it appears in the data set that matches
    Patient1's data.
    """
    chromosome_map = {}
    data = pd.read_csv(PATIENT1)
    last_chromosome = None
    for record in data.iterrows():
        if record[1].loc["Chromosome"] != last_chromosome:
            if last_chromosome:
                chromosome_map[last_chromosome] = (start, end)
            start = record[0]
            end = start
            last_chromosome = record[1].loc["Chromosome"]
        else:
            end += 1
    chromosome_map[last_chromosome] = (start, end)
    return chromosome_map

def switch_strand_direction(string):
    """
    :return: The complementary strand of string
    """
    string = string[::-1]
    string = string.replace("A", "M")
    string = string.replace("T", "A")
    string = string.replace("M", "T")
    string = string.replace("G", "M")
    string = string.replace("C", "G")
    string = string.replace("M", "C")
    return string


def confirm_gene_names():
    """
    This function checks that all the gene names that appear in the Data File (patient1 samples)
    also appear in the reference genome Fasta file
    :return: list of genes that appear in the Data but not in the fasta file
    """
    genes = pd.read_csv("C:\\Users\\Ella\\Google Drive\\Lab\\gene_table.csv", usecols=[0], names=["names"])
    genes = genes.names.tolist()
    unmatched = []
    genes_from_data = []
    data = pd.read_csv(PATIENT1, usecols=[20, 21], names=[0, 1])
    for row in range(1, len(data)):
        if data[1][row] == "upstream" or data[1][row] == "downstream":
            continue
        elif data[0][row] not in genes_from_data:
            genes_from_data.append(data[0][row])
    for gene in genes_from_data:
        if gene not in genes:
            unmatched.append(gene)
    return unmatched
