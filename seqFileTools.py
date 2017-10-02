#! /usr/local/python_anaconda/bin/python3.4

import os
from os import path
from Bio import SeqIO
from Bio import AlignIO
from Bio import Alphabet
from Bio.Alphabet import IUPAC
from Bio.Seq import Seq
from Bio.Align import MultipleSeqAlignment
import re
from file_utilities import check_filename


def unalign(filename, in_format="fasta", gap = "-", outfile = None):
    """
    unaligns file
    :param filename: input alignment filename
    :param in_format: input format (default: fasta)
    :param gap: gap type (default: - )
    :return: out file path without gaps
    """
    filename = check_filename(filename)
    alignment = AlignIO.read(filename, in_format, alphabet=Alphabet.Gapped(IUPAC.unambiguous_dna))
    for seq in alignment:
        seq.seq = seq.seq.ungap(gap)
    if outfile == None:
        outfile = path.splitext(filename)[0] + "-unaligned.fasta"
    else:
        outfile = check_filename(outfile, Truefile=None)
    SeqIO.write(alignment, outfile, "fasta")
    print("saved unaligned %s" % outfile)
    return outfile

def convert_fasta_to_phylip(filename, outfile = None):
    if outfile == None:
        outfile = os.path.splitext(filename)[0] + ".phy"
    filename = open(filename, "r")
    fasta = filename.readlines()
    phylip = open(outfile, "w")
    count = 0
    subheader = []
    subcontent = []

    for line in fasta:
        line = line.rstrip()
        if ">" in line:
            line = line.split(">", 1)[1]
            subheader.append(line)
            subcontent.append("")
            count += 1
        else:
            subcontent[count-1] = subcontent[count-1] + line
    filename.close()

    content_length = len(subcontent[0])
    split_content = list(subcontent[0])
    for i in range(content_length):
        if split_content[i] == " ":
            content_length = content_length - 1

    phylip.write("    %i    %i\n" % (count, content_length))
    for i in range(count):
        length = 10 - len(subheader[i])
        phylip.write(subheader[i])
        phylip.write(" " * length)
        phylip.write("  %s\n" %subcontent[i])
    phylip.close()
    return outfile


def format_changer(filename, out_format, outfile= None, in_format="fasta"):
    """
    sequence file format changer
    :param filename: input sequence filename
    :param out_format: output format
    :param outfile: output file (default: None)
    :param in_format: input format (default: fasta)
    :return: out file path in out format
    """
    filename = check_filename(filename)
    if outfile != None:
        outfile = check_filename(outfile, Truefile=False)
    else:
        outfile = path.splitext(filename)[0] + "." + out_format
    alignment = AlignIO.read(filename, in_format, alphabet=Alphabet.Gapped(IUPAC.unambiguous_dna))
    AlignIO.write(alignment, outfile, out_format)
    print("saved %s in format %s" % (outfile, out_format))
    return outfile


def translate_file(filename, outfile = None, in_format="fasta"):
    """
    translates nucleotide seq file to amino acid seq file
    :param filename: input nucleotide sequence filename
    :param outfile: output file (default: None)
    :param in_format: input format (default: fasta)
    :return: out file path of amino acid seq
    """
    filename = check_filename(filename)
    dataset = list(SeqIO.parse(filename, in_format))
    for seq in dataset:
        print (seq.name, len(seq.seq))
        seq.seq = seq.seq.translate(stop_symbol="*")
    if outfile != None:
        outfile = check_filename(outfile, Truefile=False)
    else:
        outfile = path.splitext(filename)[0] + "-translated.fasta"
    SeqIO.write(dataset, outfile, "fasta")
    print("saved translated file in %s" % outfile)
    return outfile


def get_first_seq_file(filename, outfile=None, in_format="fasta"):
    """
    saves the first sequence of the input file in the output file
    :param filename: input sequence filename
    :param outfile: output file (default: None)
    :param in_format: input format (default: None)
    :return: output filename path, first sequence id. return None if file in empty
    """
    filename = check_filename(filename)
    if outfile != None:
        outfile = check_filename(outfile, Truefile=False)
    else:
        outfile = path.splitext(filename)[0] + "_first_seq.fasta"
    dataset = list(SeqIO.parse(filename, in_format))
    if not dataset == []:
        new_dataset = dataset[0]
        SeqIO.write(new_dataset, outfile, "fasta")
        print("saved file with first seq %s" % outfile)
        return outfile, dataset[0].id
    else:
        print("empty file, didn't save first seq. %s" % outfile)
        return None, None


def two_fasta_file_merger_without_duplication(first_file, second_file, output_file):
    """
    merges two fasta files - without duplicates ids
    :param first_file: first sequence fasta filename
    :param second_file: second sequence fasta filename
    :param output_file: output filename
    :return: output filname path of merged fasta file.
    """
    first_file = check_filename(first_file)
    second_file = check_filename(second_file)
    first_dataset = list(SeqIO.parse(first_file, "fasta"))
    first_ids = [seq.id for seq in first_dataset]
    second_dataset = list(SeqIO.parse(second_file, "fasta"))
    merged_dataset_with_duplicates = first_dataset + second_dataset
    #remove duplicate sequences
    merged_dataset_without_duplicates = []
    duplicate_ids = []
    for i in range(0, len(merged_dataset_with_duplicates)):
        if merged_dataset_with_duplicates[i].id not in duplicate_ids:
                merged_dataset_without_duplicates.append(unalign_seq(merged_dataset_with_duplicates[i]))
        if merged_dataset_with_duplicates[i].id in first_ids:
            duplicate_ids.append(merged_dataset_with_duplicates[i].id)
    SeqIO.write(merged_dataset_without_duplicates, output_file, "fasta")
    print("saved merged fasta files in %s" % output_file)
    return output_file

def two_fasta_file_merger(first_file, second_file, output_file):
    """
    merges two fasta files
    :param first_file: first sequence filename
    :param second_file: second sequence filename
    :param output_file: output filename
    :return: ouput filename path of merged fasta file.
    """
    first_file = check_filename(first_file)
    second_file = check_filename(second_file)
    first_file = open(first_file, "r").read()
    second_file = open(second_file, "r").read()
    joined_file = first_file + second_file

    output = open(output_file, "w")
    output.write(joined_file)
    output.close()
    print("saved merged fasta files in %s" % output_file)
    return output_file


def split_fasta_file_per_seq(filename, remove = True):
    """
    splits a fasta file to individual fasta files for each seq
    :param filename: input sequence filename
    :param remove: to remove original file (default True)
    :return: list of output files
    """
    filename = check_filename(filename)
    records = re.split(">", open(filename, "r").read())[1:]
    count = 0
    suffix = os.path.splitext(filename)[0]
    files = []
    for record in records:
        new_filename = suffix + "_" + str(count) + ".fasta"
        files.append(new_filename)
        new_filename = open(new_filename, "w")
        new_filename.write(">" + record)
        new_filename.close()
        count += 1
    if remove:
        os.remove(filename)
    print("saved %i files" % count)
    return files

def remove_description(filename, outfile = None):
    """
    remove seq description from fasta file.
    :param filename: input fasta file
    :param outfile: output fasta file
    :return: output filename without seq description
    """
    filename = check_filename(filename)
    records = list(SeqIO.parse(filename, "fasta"))
    for record in records:
        record.description = ""
        record.id = record.id.split("_")[0]
    if outfile == None:
        outfile = path.dirname(filename) + "/" +  path.basename(filename) + "-NO_DESCRIPTION.fatsa"
    else:
        outfile = check_filename(outfile, Truefile=False)
    SeqIO.write(records, outfile, "fasta")
    return outfile


def cut_alignemnt_by_coordinates(aln_file, coor=[], perfix="cut", in_format="fasta"):
    """
    cuts alignment file by sequnce coordinate
    attention - the coordinates must be normelized to the specific alignment
    :param aln_file: input alignment file
    :param coor: input coordinates (default: [])
    :param perfix: perfix for output file (default: cut)
    :param in_format: input alignment formar (default: fasta)
    :return: output filename of cut alignment
    """
    if coor == []:
        raise Exception("no coordinates")
    aln_file = check_filename(aln_file)
    aln = AlignIO.read(aln_file, in_format, alphabet=Alphabet.Gapped(IUPAC.unambiguous_dna))
    new_aln = aln[:, coor[0]:coor[1]]
    output = aln_file.split(".aln")[0] + "_%s.aln" % perfix
    AlignIO.write(new_aln, output, "fasta")
    print("wrote cut alignemnt in %s" % output)
    return output


def upper_case_seq_file(filename, in_format="fasta", outfile = None):
    """
    changes the sequence to uppercase
    if outfile == None rewrites on the same file
    :param filename: input sequence file
    :param in_format: input sequence format (default: fasta)
    :param output: output file (default: None)
    :return: output filename
    """
    filename = check_filename(filename)
    dataset = list(SeqIO.parse(filename, in_format))
    for seq in dataset:
        seq.seq = Seq(str(seq.seq).upper())
    if outfile == None:
        outfile = filename
    else:
        output = check_filename(outfile, Truefile=False)
    SeqIO.write(dataset, outfile, "fasta")
    return outfile


def numerate_fasta_file(filename, in_format = "fasta", outfile = None):
    """
    change seq ids to numbers in increasing order
    if output == None, overwrites  the original file
    :param filename: input sequence filename
    :param in_format: input format (default: fasta)
    :param output: output filename (default: None)
    :return:
    """
    filename = check_filename(filename)
    records = list(SeqIO.parse(filename, "fasta"))
    num = 1
    for record in records:
        record.description = ""
        record.id = str(num)
        record.name = str(num)
        num += 1
    if outfile == None:
        outfile = filename
    else:
        outfile = check_filename(outfile, Truefile=False)
    SeqIO.write(records, outfile, "fasta")


def check_duplicate_seqs(filename, in_format="fasta"):
    """
    checks if there are duplicated sequences in a sequence file
    :param filename: input sequence filename
    :param in_format: input format (default: fasta)
    :return: list of tupples of duplicate sequences
    """
    filename = check_filename(filename)
    duplicate_ids = []
    dataset = list(SeqIO.parse(filename, in_format))
    for i in range(len(dataset)):
        first = str(dataset[i].seq.ungap("-"))
        for j in range(len(dataset)):
            if i == j:
                continue
            second = str(dataset[j].seq.ungap("-"))
            if first == second:
                print("DUPLICATE: %s - %s" % (dataset[i].id, dataset[j].id))
                duplicate_ids.append((dataset[i].id, dataset[j].id))
    return duplicate_ids


def remove_gapped_positions(aln_file, output = None, in_format = "fasta"):
    """
    removes positions in an alignment which are all gapped
    if output == None - rewrites on the input file
    :param aln_file: input alignment file path
    :param output: output file path (default: None)
    :param in_format: input format (default: fatsa)
    :return: ouptut file path
    """
    aln_file = check_filename(aln_file)
    if output == None:
        output = aln_file
    else:
        output = check_filename(output, Truefile=False)
    aln = AlignIO.read(aln_file, in_format, alphabet=Alphabet.Gapped(IUPAC.unambiguous_dna))
    new_aln = None
    for i in range(len(aln[0])):
        position = aln[:, i]
        if "".join(set(position)) != "-":
            if new_aln == None:
                new_aln = aln[:, i:i+1]
            else:
                new_aln = new_aln + aln[:, i:i+1]

    AlignIO.write(new_aln, output, "fasta")


def sam_to_fasta(input, output=None):
    """
    covert sam file format to fasta file format.
    if output not provided - output will be the input filename.fasta (without the .sam extension)
    :param input: sam file path
    :param output: outpu file name (default:None)
    :return: nothing
    """
    #-v OFS=\\\n
    input = check_filename(input)
    if output == None:
        output = input.split(".sam")[0] + ".fasta"
    output = check_filename(output, Truefile=False)

    os.system("cat %s | awk -F \\\t  -v OFS=delimiter '{print $1,$10}' > %s" % (input, output))
    out_fasta = open(output, "r").read()
    out_fasta = out_fasta.replace("delimiter", "\n")
    out_fasta = re.sub("\nM", "\n>M", out_fasta)
    out_fasta = ">" + out_fasta

    out = open(output, "w")
    out.write(out_fasta)
    out.close()
