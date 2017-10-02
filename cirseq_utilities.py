


import re
import glob
import pandas as pd
import numpy as np
import time
from Bio.Seq import Seq
from Bio import Entrez
from Bio import SeqIO



def parse_reads(freqs):

    ''' this method returns a vector of reads corresponding to genome positions.
input:
        freqs file
output:
        an integer vector containing for each position in the genome it's num of reads.
'''

    path = freqs
    df = pd.read_csv(path, sep='\t')

    # remove all duplicates from Pos except the first occurrence
    # remove all x.number duplicates
    df[["Pos"]] = df[["Pos"]].astype(int)
    df = df.drop_duplicates("Pos")

    pos = df["Pos"]  # a vector of all positions
    reads = df["Read_count"]
    med_reads = reads.median()
    print(med_reads)
    return pos, reads


def get_repeats_num(tmp_cirseq_dir, out_dir):
    """
    :param tmp_cirseq_dir: tmp directory path of the cirseq pipeline analysis
    :param out_dir: the output directory path
    :return: a dictionary of the repeat stats from the cirseq pipeline analysis
    """
    files = glob.glob(tmp_cirseq_dir + "/*.fasta.blast.freqs.stats")
    repeat_summery = {}
    for file in files:
        pattern = re.compile("(\d+\t{1}\d+\n{1})", re.MULTILINE)
        text = open(file, "r").read()
        reads = pattern.findall(text)
        for r in reads:
            key = int(r.split("\t")[0])
            value = int(r.split("\t")[1].split("\n")[0])
            if key in repeat_summery:
                repeat_summery[key] += value
            else:
                repeat_summery[key] = value
    np.save(out_dir + '/repeat_summery.npy', repeat_summery)
    print("repeat_summery.npy is in your folder")
    return repeat_summery


def get_read_and_repeat_length(in_dir, out_dir):
    """
    :param in_dir: directory path of the cirseq pipeline analysis
    :param out_dir: the output directory path
    :return: DataFrame of the freqs file with the reads and repeats length, and saves it into csv file
    """
    files = glob.glob(in_dir + "/*.fasta.blast")
    arr_names = ["sseqid", "qstart", "qend", "sstart", "send", "sstrand", "length", "btop"]
    wdf = pd.DataFrame()
    for file in files:
        data = pd.read_csv(file, sep="\t",  header=None, names=arr_names)
        grouped_and_filtered = data.groupby("sseqid").filter(lambda x: min(x["qend"]) - max(x["qstart"]) > 0).groupby("sseqid")["length"].agg([np.count_nonzero, np.sum, np.max])
    wdf = pd.DataFrame.append(wdf, grouped_and_filtered)
    wdf.to_csv(out_dir + '/length_df.csv', sep=',', encoding='utf-8')
    print("length_df.csv is in your folder")
    return wdf


def find_mutation_type(freqs_file, ncbi_id):
    """
    This function adds Mutation type to the freqs file
    :param freqs_file:  The path of the relevant freqs file
    :return:DataFrame of the freqs file with mutation type column, save it into txt file
    """
    start_time = time.time()
    file_name = freqs_file
    data = pd.read_table(freqs_file)
    data.reset_index(drop=True, inplace=True)
    start_pos, end_pos = find_coding_region(ncbi_id)
    data = data.loc[data['Pos'] >= start_pos]
    check_12mer(data)
    data['Codon'] = ""
    data['Ref_Protein'] = ""
    data['Potential_Protein'] = ""
    data['Mutation Type'] = ""
    data['Pos'] = data['Pos'].astype(int)
    data.reset_index(drop=True, inplace=True)
    for kmer in range(data.first_valid_index(), (len(data)), 12):
        print ("going over kmer: " + str(kmer) + "/" + str((len(data))))
        kmer_df = data[:].loc[kmer:kmer+11]
        print("finding codons....")
        find_codon(kmer_df)
        print("translating codons....")
        translate_codon(kmer_df)
        print("Sets the Mutation type....")
        kmer_df['Mutation Type'] = kmer_df[['Ref_Protein', 'Potential_Protein']].apply(
            lambda protein: check_mutation_type(protein[0], protein[1]), axis=1)
    print("After a long for loop")
    file_name += ".with.mutation.type.func2.freqs"
    data.to_csv(file_name, sep='\t', encoding='utf-8')
    print("The File is ready in the folder")
    print("--- %s sec ---" % (time.time() - start_time))
    return data


def freqs_to_dataframe(freqs_file):
    """
    This function returns arranged DataFrame without deletions
    :param freqs_file: The path of the relevant freqs file
    :return: DataFrame without deletions
    """
    data = pd.read_table(freqs_file)
    data = data[data.Ref != '-']
    data = data[data.Base != '-']
    data.reset_index(drop=True, inplace=True)
    return data


def find_coding_region(ncbi_id):
    """
    :param ncbi_id: NCBI_ID number
    :return:the start and the end positions of the Coding region
    """
    try:
        Entrez.email = "A.N.Other@example.com"
        handle = Entrez.efetch(db="nucleotide", rettype="gb", retmode="text", id=ncbi_id)
        ncbi_gb = SeqIO.read(handle, "gb")
        handle.close()
        location = list(ncbi_gb.features[1].location)
        start_pos = location[0]
        end_pos = location[-4]
        return start_pos, end_pos
    except:
        print('Failed to fetch record.')


def check_12mer(data):
    """
      :param data: pandas DataFrame of the freqs file
      :return: None - checks that DataFrame is in 12 kmer
      """
    if (len(data)) % 12 != 0:
        print('The data is not in 12 mer. arrange the data')
    else:
        print('The data is in 12 mer. All OK.')
        print('The length of data is:', (len(data)))


def find_codon(data):
    """
    :param data: pandas DataFrame of the freqs file
    :return: pandas DataFrame with codon column
    """
    data['Codon'] = ""
    x = 4
    first_index_kmer = data.first_valid_index()
    for i in range(first_index_kmer, first_index_kmer+x, 1):
        data['Codon'][i] = data["Base"].loc[i] + data["Ref"].loc[first_index_kmer+x] + data["Ref"].loc[(first_index_kmer+(x*2))]
    for i in range((first_index_kmer+x), (first_index_kmer+(x*2)), 1):
            data['Codon'][i] = data["Ref"].loc[first_index_kmer] + data["Base"].loc[i] + data["Ref"].loc[first_index_kmer+(x*2)]
    for i in range((first_index_kmer+(x*2)), (first_index_kmer+(x * 3)), 1):
            data['Codon'][i] = data["Ref"].loc[first_index_kmer] + data["Ref"].loc[first_index_kmer+x] + data["Base"].loc[i]
    return data


def translate(seq):
    """
    :param seq: string of ATGC
    :return: protein sequence of seq
    """
    seq = Seq(str(seq))
    protein = seq.translate()
    return protein


def translate_codon(data):
    """
    :param data: pandas DataFrame of the freqs file
    :return: pandas DataFrame with Reference protein column and Potential protein column
    """
    data['Ref_Protein'] = ""
    data['Potential_Protein'] = ""
    protein_ref = Seq(data["Codon"][data.first_valid_index()]).translate()
    data['Ref_Protein'].loc[data.first_valid_index():data.first_valid_index()+11] = protein_ref
    data['Ref_Protein'].replace('(', '', inplace=True)
    data['Ref_Protein'].replace(')', '', inplace=True)

    data["Potential_Protein"] = data["Codon"].apply(translate)
    data["Potential_Protein"].replace('(', '', inplace=True)
    data["Potential_Protein"].replace(')', '', inplace=True)
    return data


def check_mutation_type(protein1, protein2):
    """
    :param protein1: amino acid 1
    :param protein2: amino acid 2
    :return: The mutation type
    """
    Mutation_Type = ""
    if protein1 == protein2:
        Mutation_Type = "Synonymous"
    elif protein1 != protein2:
        Mutation_Type = "Non-Synonymous"
    if protein2 == '*':
        Mutation_Type = "Premature Stop Codon"
    return Mutation_Type
