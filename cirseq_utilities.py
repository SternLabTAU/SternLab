#! /usr/local/python/anaconda_python-3.6.1


import re
import glob
import pandas as pd
import numpy as np
import time
from Bio.Seq import Seq
from Bio import Entrez
from Bio import SeqIO
import collections
import pathlib
import os




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

    #Debug
    # freqs_file = "/Volumes/STERNADILABHOME$/volume3/okushnir/AccuNGS/180503_OST_FINAL_03052018/merged/RV-p11/q30_3UTR_new/RV-p11.freqs"
    # ncbi_id = "NC_001490"

    start_time = time.time()
    file_name = freqs_file
    data = freqs_to_dataframe(freqs_file)
    data.reset_index(drop=True, inplace=True)
    orign_data = data
    start_pos, end_pos = find_coding_region(ncbi_id)
    # for half gene -  PAY ATTENTION!!!
    #start_pos = 3608


    data = data.loc[data['Pos'] >= start_pos]
    data = data.loc[data['Pos'] <= end_pos]
    if check_12mer(data) != 1:
        raise Exception("The data is not in 12 mer. arrange the data, start_pos=" + str(start_pos) + " end_pos=" + str(end_pos))
    data['Codon'] = ""
    data['Ref_AA'] = ""
    data['Potential_AA'] = ""
    data['Type'] = ""
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
        kmer_df['Type'] = kmer_df[['Ref_AA', 'Potential_AA']].apply(
            lambda protein: check_mutation_type(protein[0], protein[1]), axis=1)
    print("After a long for loop")
    top_data = orign_data.loc[orign_data['Pos'] < start_pos]
    top_data['Codon'] = ""
    top_data['Ref_AA'] = ""
    top_data['Potential_AA'] = ""
    top_data['Type'] = ""
    top_data['Pos'] = top_data['Pos'].astype(int)
    bottom_data = orign_data.loc[orign_data['Pos'] > end_pos]
    bottom_data['Codon'] = ""
    bottom_data['Ref_AA'] = ""
    bottom_data['Potential_AA'] = ""
    bottom_data['Type'] = ""
    bottom_data['Pos'] = bottom_data['Pos'].astype(int)
    frames = [top_data, data, bottom_data]
    data_final = pd.concat(frames)
    file_name = file_name[0:-5]
    file_name += "with.mutation.type.freqs"
    data_final.to_csv(file_name, sep='\t', encoding='utf-8')
    print("The File is ready in the folder")
    print("--- %s sec ---" % (time.time() - start_time))
    print("start_pos:%i" % start_pos)
    print("end_pos:%i" % end_pos)
    return data_final


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
        start_pos = 0
        end_pos = 0
        # if ncbi_id == "V01149":
        for feature in ncbi_gb.features:
            if feature.type == 'CDS':
                start_pos = feature.location.start + 1
                end_pos = feature.location.end
            #     location = list(ncbi_gb.features[5].location)
            # else:
            #     location = list(ncbi_gb.features[1].location)
            # start_pos = location[0] + 1
            # end_pos = location[-1] + 1
        return start_pos, end_pos
    except:
        print('Failed to fetch record.')


def check_12mer(data):
    """
      :param data: pandas DataFrame of the freqs file
      :return: 1/0; 1- the data in 12mer, 0 - ths data is NOT in 12mer checks that DataFrame is in 12 kmer
      """
    if (len(data)) % 12 != 0:
        # print('The data is not in 12 mer. arrange the data')
        flag = 0
    else:
        # print('The data is in 12 mer. All OK.')
        # print('The length of data is:', (len(data)))
        flag = 1
    return flag


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
    data['Ref_AA'] = ""
    data['Potential_AA'] = ""
    protein_ref = Seq(data["Codon"][data.first_valid_index()]).translate()
    data['Ref_AA'].loc[data.first_valid_index():data.first_valid_index()+11] = protein_ref
    data['Ref_AA'].replace('(', '', inplace=True)
    data['Ref_AA'].replace(')', '', inplace=True)

    data["Potential_AA"] = data["Codon"].apply(translate)
    data["Potential_AA"].replace('(', '', inplace=True)
    data["Potential_AA"].replace(')', '', inplace=True)
    return data


def check_mutation_type(aa1, aa2):
    """
    :param aa1: amino acid 1
    :param aa2: amino acid 2
    :return: The mutation type
    """
    Mutation_Type = ""
    if aa1 == aa2:
        Mutation_Type = "Synonymous"
    elif aa1 != aa2:
        Mutation_Type = "Non-Synonymous"
    if aa2 == '*':
        Mutation_Type = "Premature Stop Codon"
    return Mutation_Type


def filter_by_coverage_mutation(type_file, output_dir):
    """
    Creates DataFrame of frequencies with Mutations X->Y column
    :param type_file: freqs file after it has ben processed with find_mutation_type func
    :param output_dir: where to save the returned DF
    :return: DF with Mutation and mutation_type column with minimal coverage
    """
    #data = pd.read_table(freqs_file)
    data = type_file
    data.reset_index(drop=True, inplace=True)
    flag = '-' in data.Base.values
    if flag is True:
        data = data[data.Ref != '-']
        data = data[data.Base != '-']
        data.reset_index(drop=True, inplace=True)
        # data['abs_counts'] = data['Freq'] * data["Read_count"]  # .apply(lambda x: abs(math.log(x,10))/3.45)
        # data['Frequency'] = data['abs_counts'].apply(lambda x: 1 if x == 0 else x) / data["Read_count"]
        # raise Exception("This script does not support freqs file with deletions, for support please contact Maoz ;)"
    data['Base'].replace('T', 'U', inplace=True)
    data['Ref'].replace('T', 'U', inplace=True)
    min_read_count = 100000
    data = data[data['Pos'] == np.round(data['Pos'])]  # remove insertions
    data['Pos'] = data[['Pos']].apply(pd.to_numeric)
    data = data[data['Read_count'] > min_read_count]
    data['mutation_type'] = data['Ref'] + data['Base']
    data = data[data['Ref'] != data['Base']]
    data = data[data["Base"] != "-"]
    data['abs_counts'] = data['Freq'] * data["Read_count"]  # .apply(lambda x: abs(math.log(x,10))/3.45)
    data['Frequency'] = data['abs_counts'].apply(lambda x: 1 if x == 0 else x) / data["Read_count"]
    data["Mutation"] = data["Ref"] + "->" + data["Base"]
    data = data[["Pos", "Base", "Freq", "Ref", "Read_count", "Rank", "Prob", "Codon", "Ref_AA", "Potential_AA", "Type",
                 "abs_counts", "Frequency", "mutation_type", "Mutation", "label", "passage", "replica"]]
    data.to_csv(output_dir + "data_mutation.csv", sep=',', encoding='utf-8')
    return data



# def main():


# if __name__ == "__main__":
#     main()