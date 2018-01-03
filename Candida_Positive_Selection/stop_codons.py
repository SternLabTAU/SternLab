import pandas as pd
import os
from Bio import SeqIO
from os import path
import shutil
import seqFileTools
import PAML_utilities
from Candida_Positive_Selection import run_PAML
from Candida_Positive_Selection import create_gene_table
from optparse import OptionParser


def main():
    parser = OptionParser("usage: %prog[options]")
    parser.add_option("-f", "--fasta", dest="FASTA_FILE_PATH", help="Reference genome fasta file")
    parser.add_option("-g", "--gene_table", dest="GENE_TABLE_PATH", help="Gene table desired location")
    parser.add_option("-i", "--gene_folder", dest="GENES_FOLDER", help="The folder that contains all candida genes")

    (options, args) = parser.parse_args()

    FASTA_FILE_PATH = options.FASTA_FILE_PATH
    GENE_TABLE_PATH = options.GENE_TABLE_PATH
    GENES_FOLDER = options.GENES_FOLDER

    find_problematics(FASTA_FILE_PATH, GENE_TABLE_PATH, GENES_FOLDER)


def find_problematics(FASTA_FILE_PATH, GENE_TABLE_PATH, GENES_PATH):
    """
    Reads all genes from the referance genome fasta file, and checks for each gene if it is valid.
    It not - it is moved to the "excluded" folder
    :param FASTA_FILE_PATH:
    :param GENE_TABLE_PATH:
    :return:
    """
    gene_info = pd.read_csv(GENE_TABLE_PATH)
    gene_info.rename(columns={"Unnamed: 0": "Gene"}, inplace=True)
    gene_info["Gene Error"] = ""
    gene_info["Error info"] = ""
    cnt = 0

    with open(FASTA_FILE_PATH, "r") as file:
        for record in SeqIO.parse(file, "fasta"):
            description_line = record.description
            gene_seq = str(record.seq)
            gene_name, chromosome, start, end, direction = create_gene_table.extract_args_from_description(description_line)
            filename = chromosome + "-" + gene_name + ".fasta"

            # if it is not a gene that was modified in patient 1 (total of 826) - skip it
            if not path.isfile(GENES_PATH + filename):
                continue

            # Check if the gene seq is valid of not
            result, info = valid_seq(gene_seq.upper())
            print(result)
            if result != "valid":
                print(gene_name)
                gene_info.loc[gene_info['Gene'] == gene_name, "Gene Error"] = result
                gene_info.loc[gene_info['Gene'] == gene_name, "Error info"] = info

                #os.rename(GENES_PATH + filename,  GENES_PATH + "excluded/" + filename)

            cnt += 1

    gene_info.to_csv(GENE_TABLE_PATH, index=False)
    file.close()
    print("verified " + str(cnt) + "gene files")


def valid_seq(gene_seq):
    """
    This function recieves a gene sequence and checks if its length is not a multiple of three, or if it contains more
    than one stop codon (or if the stop codon is not in the end of the file)
    :param gene_seq: The sequence of the gene to be checked.
    :return:    "valid", None -                                         valid gene sequence
                "len problem", None -                                   invalid length (not a multiple of 3)
                "single stop- not last", position of stop codon -       in case of one stop codon that is not in the end of
                                                                        the sequence
                "no stop codons", None -                                in case there are no stop codons at all in the sequence.
                "multiple stops", [(position, stop sequence ref] -      in case there are more than one stop codon
    """
    l = len(gene_seq)
    if l % 3 != 0:
        return "len problem", l
    codon_index = 0
    stops = []
    while codon_index < l:
        if gene_seq[codon_index:codon_index + 3] in ["TAA", "TAG", "TGA"]:
            stops += [(codon_index, gene_seq[codon_index:codon_index + 3])]
        codon_index += 3
    if len(stops) == 1 and stops[0][0] != l-3:
        return "single stop- not last", stops[0][0]
    elif len(stops) == 0:
        return "no stop codons", None
    elif len(stops) > 1:
        return "multiple stops", "-".join(str(x) for x in stops)
    else:
        return "valid", None


def check_end(folder):
    """Return a list of all the sequences that have no stop codon in at least one of the samples"""
    bads = []
    for filename in os.listdir(folder):
        if filename.endswith(".fasta"):
            file = open("/sternadi/home/volume2/ella/Candida/Genes/" + filename, "r")
            # skip the reference gene
            file.readline()
            file.readline()
            for line in file:
                if ">" not in line:
                    if line[-4:-1] not in ["TAA", "TAG", "TGA"]:
                        if filename not in bads:
                            print(filename)
                            bads.append(filename)
        file.close()
    return bads


def rerun_middle_stops(folder):
    """
    This function runs over the original gene folder and for every gene that is in the list "run_again" it creates its gene
    folder empty, modifies the sequence so it has gaps instead of stop codons in the middle and then runs PAML again.
    After creating the new Fasta file it will confirm that there is a difference between the sequences.
    :param folder: The general Genes folder
    :return number of files changed
    """
    cnt = 0
    run_again = ["Ca21chr2-orf19.813","Ca21chr6-RBF1","Ca21chr4-PGA31","Ca21chr5-orf19.1935","Ca21chr4-GST1","Ca21chr2-orf19.894","Ca21chr7-orf19.5139","Ca21chr1-CHS3","Ca21chr4-orf19.2680","Ca21chr3-orf19.6008","Ca21chr1-orf19.4984","Ca21chr4-RAM1","Ca21chr2-orf19.1768","Ca21chr5-orf19.4337","Ca21chrR-orf19.1737","Ca21chrR-orf19.6382","Ca21chr1-orf19.7278","Ca21chr1-orf19.6209","Ca21chr3-FGR23","Ca21chr4-ZCF27","Ca21chr4-ERG26","Ca21chr5-orf19.937"    ]
    for filename in os.listdir(folder):
        if filename[:-6] in run_again:
            seqs = []
            cnt += 1
            gene_folder = "/sternadi/home/volume2/ella/Candida/Genes/NoRef/" + filename[:-6]
            shutil.rmtree(gene_folder)
            file = open("/sternadi/home/volume2/ella/Candida/Genes/" + filename, "r")
            #skip the reference gene
            file.readline()
            file.readline()
            os.makedirs(gene_folder)
            output = open(gene_folder + "/" + filename, "w")

            # write all sequences with A in the begining of the name, and without the stop codon at the end
            # also ommits the stop codons from the middle and completes it with gaps.
            for line in file:
                if ">" in line:
                    output.write(line.replace(">",">A"))
                else:
                    new_line = return_gapped_seq(line[:-4].upper())
                    output.write(new_line + "\n")
                    seqs += [new_line]
            output.close()
            file.close()

            #confirm all the sequences are not the same:
            if all_the_same(seqs):
                print("all the same in " + filename)
                continue

            #create phylip
            seqFileTools.convert_fasta_to_phylip(
                gene_folder + "/" + filename, outfile=None)

            #create CLT
            M8a_CTL = gene_folder + "/" + filename.replace(".fasta", "-M8a.CLT")
            M8_CTL = gene_folder + "/" + filename.replace(".fasta", "-M8.CLT")
            PAML_utilities.candida_write_ctl_codeml_file(
                M8a_CTL,
                gene_folder + "/" + filename.replace(".fasta", ".phy"),
                "/sternadi/home/volume2/ella/Candida/Trees/SNP_seq_without_ref.phy_phyml_tree.txt",
                gene_folder + "/" + filename.replace(".fasta", "-M8a-Results.txt"),
                1)
            PAML_utilities.candida_write_ctl_codeml_file(
                M8_CTL,
                gene_folder + "/" + filename.replace(".fasta", ".phy"),
                "/sternadi/home/volume2/ella/Candida/Trees/SNP_seq_without_ref.phy_phyml_tree.txt",
                gene_folder + "/" + filename.replace(".fasta", "-M8-Results.txt"),
                0)

            #run codeml
            run_PAML.candida_codeml_runner(M8a_CTL, M8_CTL, filename[:-6])

    print(str(cnt) + " files were ran again")


def mark_old_outputs(folder):
    """This marks all the outputs of files that had stop codons in the middle,
    before running PAML on them again with gaps"""
    stops_old_cmls = ["cml.o4917128","cml.o4916712","cml.o4917318","cml.o4917112","cml.o4917041","cml.o4917143","cml.o4916972","cml.o4916724","cml.o4917240","cml.o4916583","cml.o4916750","cml.o4916988","cml.o4916882","cml.o4916888","cml.o4916873","cml.o4916958","cml.o4917305","cml.o4916831","cml.o4916731","cml.o4917248","cml.o4916699","cml.o4917266"]
    for f in os.listdir(folder):
        if f in stops_old_cmls:
            p = path.abspath(f)
            os.rename(folder +f , folder +f +"-old")


def return_gapped_seq(gene_seq):
    """Recieves a sequence, finds its first stop codon, and completes it with gaps"""
    l = len(gene_seq)
    codon_index = 0
    while codon_index < l:
        if gene_seq[codon_index:codon_index + 3] in ["TAA", "TAG", "TGA"]:
            break
        codon_index += 3
    new_seq = gene_seq[:codon_index] + "-"*(l-codon_index)
    return new_seq


def all_the_same(lst):
    """Check if all the sequences in the list are identical"""
    for seq in lst:
        if seq != lst[0]:
            return False
    return True


if __name__ == "__main__":
    main()