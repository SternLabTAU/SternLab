from Candida_Positive_Selection import create_gene_table
from Candida_Positive_Selection import utilities
import pandas as pd
from Bio import SeqIO


def create_gene_files(FASTA_FILE_PATH, GENE_TABLE_PATH, PATIENT1,PATH_PREFIX, ERRORS):
    """
    Creates a new file for each gene, that contains its versions in each of the 17 samples.
    It check if the gene length is a multiple of three and if there is only one codon at the end- if not, it exclued it
    and adds it to the file fo excluded gene with the reason for being excluded.
    For each gene, the list loc_for_change is created - it contains all positions in the gene that need to be altered in
    at least one of the samples, or is empty in case all the samples have the exact same gene sequence.
    Then, if loc_for_change is not empty, the function will create a new file for the gene and
    start adding all of its 17 versions according to the information in Patient1's data.
    """
    #   add new column for the changed positions of each gene
    gene_info = pd.read_csv(GENE_TABLE_PATH)
    gene_info.rename(columns={"Unnamed: 0": "Gene"}, inplace=True)
    gene_info["Changed Positions"] = "None"
    gene_info["Changed Positions Count"] = 0
    gene_info["Gene Error"] = ""
    gene_info["Error info"] = ""

    with open(FASTA_FILE_PATH, "r") as file:
        for record in SeqIO.parse(file, "fasta"):
            description_line = record.description
            print(description_line)
            gene_seq = str(record.seq)
            gene_name, chromosome, start, end, direction = create_gene_table.extract_args_from_description(description_line)

            loc_for_change = find_changed_locations(chromosome, start, end, PATIENT1)
            if loc_for_change == ["identical"]:
                gene_info.loc[gene_info['Gene'] == gene_name, "Changed Positions"] = "Identical gene in all samples"
            elif len(loc_for_change) != 0:
                # update Changed Positions and Changed Positions count in Gene Table
                gene_info.loc[gene_info['Gene'] == gene_name, "Changed Positions"] =\
                    "-".join(str(x) for x in loc_for_change)
                gene_info.loc[gene_info['Gene'] == gene_name, "Changed Positions Count"] = len(loc_for_change)

                write_new_fasta(gene_name, chromosome, start, end, direction, loc_for_change,
                                gene_seq, PATH_PREFIX, PATIENT1, ERRORS)

    gene_info.to_csv(GENE_TABLE_PATH, index=False)
    file.close()


def find_changed_locations(chromosome, gene_start, gene_end, PATIENT1):
    """
    Goes over Patient1's data and returns a list of all the genome locations that needs modification within the
    gene who's chromosome and location were passed as arguments.
    If all the changes compared to the reference genome are the same in all of the samples, it will return an empty
    list, because that gene is not relevant to the FoG phenotype.
    :param chromosome: The gene's chromosome extracted from the reference fasta file
    :param gene_start: The gene's start position extracted from the reference fasta file
    :param gene_end: The gene's end position extracted from the reference fasta file
    :return: A list of all positions that needs modification within the gene.
    """
    data = pd.read_csv(PATIENT1)
    data.columns = ["chr", "pos"] + list(data.columns)[2:]
    data = data.loc[(data["chr"] == chromosome) & \
                    (data["pos"] <= gene_end) & (data["pos"] >= gene_start)]
    # check if there is at least one sample that is different than the others
    diff_between_samples = False
    sample1 = [x.upper() for x in list(data["1"])]
    for i in [x for x in range(2,18) if x != 10]:
        if sample1 != [x.upper() for x in list(data[str(i)])]:
            diff_between_samples = True

    if diff_between_samples:
        return list(data["pos"])
    else:
        return ["identical"]


def write_new_fasta(gene, chromosome, start, end, direction, loc_for_change, gene_seq, PATH_PREFIX, PATIENT1, ERRORS):
    """
    Creates a new file and writed into it the original version of the gene and its modification for each of the 17
    samples.
    :param gene: The name of the gene
    :param chromosome: The chromosome in which the gene appears
    :param start: The gene's start location on the chromosome
    :param end: The gene's end location on the chromosome
    :param direction: W or C strand
    :param loc_for_change: All the positions of the gene that need modification in at least one of the sratnds
    :param gene_seq: The gene full sequence, as a string
    :return:
    """
    new_path = PATH_PREFIX + chromosome + "-" + gene + ".fasta"
    new_file = open(new_path, "w")
    new_file.write(">0\n")
    new_file.write(gene_seq)
    new_file.write("\n")
    for sample in (x for x in range(1, 18) if x != 10):
        sample_specific_gene = generate_sample_specific_gene\
            (sample, gene_seq, direction, start, chromosome, loc_for_change, PATIENT1, ERRORS)
        new_file.write(">" + str(sample) + "\n")
        new_file.write(sample_specific_gene)
        new_file.write("\n")
    new_file.close()


def generate_sample_specific_gene(sample, gene_sec, direction, start, chromosome, loc_for_change, PATIENT1, ERRORS):
    """
    Generates the correct sequence for the sample according to the information from Patient1's data set.
    :param sample: The sample number between 1 to 17
    :param gene_sec: The original gene sequence from the reference genome
    :param direction: W or C strand
    :param start: The start position of the gene in the chromosome
    :param chromosome: The chromosome where the gene is located in
    :param loc_for_change: All locations in the gene that needs modification in at least one of the samples
    :param PATIENT1: The path to the data of patient 1
    :return: The correct sequence of the sample
    """
    patient1_data = pd.read_csv(PATIENT1)
    new_cols = ["chr", "pos"] + list(patient1_data.columns)[2:]
    patient1_data.columns = new_cols
    #   if the seq is the C strand, turn it to the W strand, modify it an then turn back to C
    modified_gene = utilities.switch_strand_direction(gene_sec) if direction == "C" else gene_sec
    for location in loc_for_change:
        position_to_change = location - start
        modification_line = patient1_data.loc[(patient1_data["chr"] == chromosome) & (patient1_data["pos"] == location)].reset_index()
        # extract the base of the reference genome and confirm it matches the expected reference base from Patient1's data.
        # if not, add this information to the ERRORS file
        reference_in_loc = modification_line.get_value(0, "Reference Base (SC5314)")
        if modified_gene[position_to_change].upper() != reference_in_loc and sample == 1:
            print(modified_gene[position_to_change], reference_in_loc + "\n")
            errors = open(ERRORS, "w")
            log = "Mismatch between reference genome and Patient1's data in location: " + location + "\n" +\
                  modified_gene[position_to_change] + "and" + reference_in_loc + "\n"
            errors.write(log)
            errors.write("gene seq:\n")
            errors.write(gene_sec)
            errors.close()

        modification = modification_line.get_value(0, str(sample))
        #   exclude cases of n,n N,N and - modifications
        if len(modification) == 1 and modification != "-":
            modified_gene = modified_gene[:position_to_change] + modification.upper() + modified_gene[position_to_change+1:]

    modified_gene = utilities.switch_strand_direction(modified_gene) if direction == "C" else modified_gene

    return modified_gene



