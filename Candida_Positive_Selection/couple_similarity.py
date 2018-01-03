import pandas as pd

def compare_couple(sample1, sample2,
                   position_file = '/sternadi/home/volume2/ella/Candida/positions_info.csv',
                   gene_seq_prefix =  '/sternadi/home/volume2/ella/Candida/Genes/NoRef/'):
    """" Returns the number of positions where sample1 and sample2 are the same between all the positively
    selected positions.
    :position_file: the path to the file that contains all positively selected genes.
    :gene_seq_prefix: The path to the folder that contains the gene specific subfolders
    """
    position_data = pd.read_csv(position_file)
    similar = 0
    for line in position_data.iterrows():
        chromosome = line[1]["Chromosome"]
        gene = line[1]["Gene"]
        pos = line[1]["Position"]
        gene_file = list(open(gene_seq_prefix + chromosome + "-" + gene + "/" + chromosome + "-" + gene + ".fasta"))
        sample1_codon = extract_codon_from_sample(gene_file, pos, sample1)
        sample2_codon = extract_codon_from_sample(gene_file, pos, sample2)
        if sample1_codon == sample2_codon:
            similar += 1
    return similar


def extract_codon_from_sample(gene_file, position, sample):
    """
    Return the codon in "position" of the gene file given as a list.
    :param gene_file: The list of the read fasta file
    :param position: Codon position to be returned
    :param sample: Sample of interest
    :return: Tree letter codon
    """
    if sample < 10:
        seq = gene_file[2*sample-1]
    else:
        seq = gene_file[2*sample - 2 -1]
    start = (position-1) * 3
    return seq[start:start+3]


def compare_all_couples(output_path):
    all_couples_results = pd.DataFrame(columns=["Couple", "Similarity_score"])
    for i in range(1,18):
        for j in range(1,18):
            if i != 10 and j != 10 and i < j:
                val = compare_couple(i, j)
                new_line = {
                    "Couple": [(i,j)],
                    "Similarity_score": [val]}
                new_line = pd.DataFrame(new_line, columns=["Couple", "Similarity_score"])
                all_couples_results= all_couples_results.append(new_line)
    all_couples_results.to_csv(output_path)
