import pandas as pd


def create_gene_table(FASTA_FILE_PATH, GENE_TABLE_PATH):
    """
    This function builds a new csv file of all the genes from the reference genome with and listing of
    their chromosomal locations.
    :return:
    """

    table = extract_data_from_fasta(FASTA_FILE_PATH)
    table.to_csv(GENE_TABLE_PATH)


def extract_data_from_fasta(FASTA_FILE_PATH):
    """
    This function extracts the gene specific data (location + chromosome + direction) from the fasta file.
    :return: A DataFrame of genes with their chromosomes and chromosomal locations.
    """
    #   read info into gene_data map
    gene_data = {}
    with open(FASTA_FILE_PATH, "r") as file:
        for line in file:
            if ">" in line:
                gene_name, chromosome, start, end, direction = extract_args_from_description(line)
                gene_data[gene_name] = (chromosome, start, end, direction)

    file.close()

    #   create DataFrame and fill it with the information from the map
    table = pd.DataFrame(index=gene_data.keys(), columns=["Chromosome", "Start", "End", "Direction"])
    table.columns.name = "Gene"
    for gene in gene_data.keys():
        table["Chromosome"][gene] = gene_data[gene][0]
        table["Start"][gene] = gene_data[gene][1]
        table["End"][gene] = gene_data[gene][2]
        table["Direction"][gene] = gene_data[gene][3]
    return table


def extract_args_from_description(line):
    """
    Extracts the information about a gene from its descriptive line in the reference genome file.
    :param line: The description line of the fasta file of the reference genome
    :return: gene_name, chromosome, start, end, direction (W or C strand)
    """
    gene_name = line.split(" ")[1]
    chromosome = line.split(" ")[5].split(":")[0]
    start = int((line.split(" ")[5].split(":")[1].split("-"))[0])
    end = line.split(" ")[5].split(":")[1].split("-")[1]
    direction = end[-1]
    end = int(end[:len(end) - 1])
    return gene_name, chromosome, start, end, direction
