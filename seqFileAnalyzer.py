#! /usr/local/python_anaconda/bin/python3.4


from Bio import SeqIO
from Bio import AlignIO
from Bio import Alphabet
from Bio.Alphabet import IUPAC
import collections
from file_utilities import check_filename
import pandas as pd

def count_gaps_and_characters(aln_file, file_format = "fasta"):
    """
    count how many gaps and how many characters there are in an alignemnt
    :param aln_file: input alignment file
    :param file_format: input file format (default: fasta)
    :return: alignment length, number of gap chars, number of non-gap chars
    """
    aln_file = check_filename(aln_file)
    aln = AlignIO.read(aln_file, file_format, alphabet=Alphabet.Gapped(IUPAC.unambiguous_dna))
    total_gaps = 0
    total_not_gaps = 0
    for record in aln:
        local_gaps = record.seq.count("-")
        local_not_gaps = len(record.seq) - local_gaps
        total_gaps += local_gaps
        total_not_gaps += local_not_gaps
    return len(aln), total_gaps, total_not_gaps


def base_frequencies(filename, in_format="fasta"):
    """
    calculates base frequencies in sequence file
    :param filename: input nucleotide sequence filename
    :param in_format: input format (default: fasta)
    :return: freqs dictionary
    """
    filename = check_filename(filename)
    dataset = list(SeqIO.parse(filename, in_format))
    freqs = {"A": 0, "G": 0, "C": 0, "T": 0}
    count = 0
    for seq in dataset:
        a = seq.seq.count("A")
        c = seq.seq.count("C")
        t = seq.seq.count("T")
        g = seq.seq.count("G")
        count += len(seq.seq)
        freqs["A"] += a
        freqs["C"] += c
        freqs["T"] += t
        freqs["G"] += g
    for k in freqs.keys():
        freqs[k] = freqs[k] / float(count)
    print(freqs)
    return freqs


def get_consensus_from_alignment(aln_file, in_format="fasta"):
    """
    constructs a consensus sequence from alignment file
    :param aln_file: alignment file
    :param in_format: file format (default: fasta)
    :return: consensus sequence
    """
    aln_file = check_filename(aln_file)
    aln = AlignIO.read(aln_file, in_format, alphabet=Alphabet.Gapped(IUPAC.unambiguous_dna))
    len_aln = len(aln[0])
    consensus = ""
    for i in range(len_aln):
        count = 0
        max_char = ""
        counter = collections.Counter(aln[:, i])
        for j in counter:
            if counter[j] > count:
                count = counter[j]
                max_char = j
        if max_char == "-":
            continue
        consensus += max_char
    return consensus


def stop_mutation_potential_in_coding_sequence(filename, in_format="fasta"):
    """
    checks the potential to create stop codons in a given sequence
    assumes that coding region starts from the first nucleotide
    :param filename: input file name
    :param in_format: input format
    :return: dictionary of sequence name and stop mutation count
    """
    stop_codons = ["TGA", "TAA", "TAG"]
    filename = check_filename(filename)
    dataset = list(SeqIO.parse(filename, in_format))
    df = pd.DataFrame(columns=["seq_id", "seq_len", "stop_mutation_count"])
    for seq in dataset:
        stop_mutation_count = 0
        seq_len = len(seq.seq)
        for i in range(0, seq_len, 3):
            codon = seq.seq[i:i + 3]
            if codon in stop_codons:
                continue
            for nuc in ["T", "C", "A", "G"]:
                new_codon_1 = codon[:2] + nuc
                new_codon_2 = nuc + codon[1:]
                new_codon_3 = codon[0] + nuc + codon[2]
                if new_codon_1 in stop_codons:
                    stop_mutation_count += 1
                if new_codon_2 in stop_codons:
                    stop_mutation_count += 1
                if new_codon_3 in stop_codons:
                    stop_mutation_count += 1
        df = df.append({"seq_id":seq.id, "seq_len":seq_len, "stop_mutation_count":stop_mutation_count},
                       ignore_index=True)
    return df


def get_consensus_percentage(aln_file, in_format="fasta"):
    """
    gets alignment file and returns the consensus and
    the percentage of each position in the alignment
    the percentage calculation ignores gaps
    :param aln_file: input alignment file path
    :param in_format: input file format (defualt: fasta)
    :return: consensus sequance and consensus percentage
    """
    aln_file = check_filename(aln_file)
    aln = AlignIO.read(aln_file, in_format, alphabet=Alphabet.Gapped(IUPAC.unambiguous_dna))
    len_aln = len(aln[0])
    num_of_seq = len(aln)
    consensus_percentage= {1:0, 0.9:0, 0.8:0, 0.7:0, 0.6:0, 0.5:0, 0.4:0, 0.3:0, 0.2:0}
    consensus = ""
    for i in range(len_aln):
        counter = collections.Counter(aln[:, i])
        count = 0
        max_char = ""
        for j in counter:
            if j == "-":
                continue
            elif counter[j] > count:
                count = counter[j]
                max_char = j
        if "-" not in counter:
            gap_count = 0
        else:
            gap_count = counter["-"]
        percentage = round(count/(num_of_seq-gap_count), 1)
        consensus_percentage[percentage] += 1
        consensus += max_char

    for n in consensus_percentage:
        consensus_percentage[n] = round(consensus_percentage[n] / len_aln, 3)
    return consensus, consensus_percentage


def get_codon_freqs_from_consensus(filename, in_format="fasta"):
    """
    gets codon frequencies of concensus sequence from sequence file
    :param filename: input sequence filename
    :param in_format: input format (default: fasta)
    :return: codons dictionary and text for simulation input
    """
    filename = check_filename(filename)
    dataset = list(SeqIO.parse(filename, in_format))
    consensus = get_consensus_from_alignment(filename)
    codons = {"TTT": 0, "TTC": 0, "TTA": 0, "TTG": 0,
              "TCT": 0, "TCC": 0, "TCA": 0, "TCG": 0,
              "TAT": 0, "TAC": 0, "TAA": 0, "TAG": 0,
              "TGT": 0, "TGC": 0, "TGA": 0, "TGG": 0,
              "CTT": 0, "CTC": 0, "CTA": 0, "CTG": 0,
              "CCT": 0, "CCC": 0, "CCA": 0, "CCG": 0,
              "CAT": 0, "CAC": 0, "CAA": 0, "CAG": 0,
              "CGT": 0, "CGC": 0, "CGA": 0, "CGG": 0,
              "ATT": 0, "ATC": 0, "ATA": 0, "ATG": 0,
              "ACT": 0, "ACC": 0, "ACA": 0, "ACG": 0,
              "AAT": 0, "AAC": 0, "AAA": 0, "AAG": 0,
              "AGT": 0, "AGC": 0, "AGA": 0, "AGG": 0,
              "GTT": 0, "GTC": 0, "GTA": 0, "GTG": 0,
              "GCT": 0, "GCC": 0, "GCA": 0, "GCG": 0,
              "GAT": 0, "GAC": 0, "GAA": 0, "GAG": 0,
              "GGT": 0, "GGC": 0, "GGA": 0, "GGG": 0}

    all_codons = 0
    for i in range(0, len(consensus), 3):
        c = str(consensus[i:i + 3])
        if len(c) < 3:
            continue
        if "N" in c:
            continue
        codons[c] += 1
        all_codons += 1

    for c in codons.keys():
        codons[c] = float(codons[c]) / all_codons
        if codons[c] == 0:
            codons[c] = 0.000001

    to_print = ""
    to_print += "%f %f %f %f // TTT TTC TTA TTG\n" % (codons["TTT"], codons["TTC"], codons["TTA"], codons["TTG"])
    to_print += "%f %f %f %f // TCT TCC TCA TCG\n" % (codons["TCT"], codons["TCC"], codons["TCA"], codons["TCG"])
    to_print += "%f %f %f %f // TAT TAC TAA TAG\n" % (codons["TAT"], codons["TAC"], codons["TAA"], codons["TAG"])
    to_print += "%f %f %f %f // TGT TGC TGA TGG\n" % (codons["TGT"], codons["TGC"], codons["TGA"], codons["TGG"])
    to_print += "\n"

    to_print += "%f %f %f %f // CTT CTC CTA CTG\n" % (codons["CTT"], codons["CTC"], codons["CTA"], codons["CTG"])
    to_print += "%f %f %f %f // CCT CCC CCA CCG\n" % (codons["CCT"], codons["CCC"], codons["CCA"], codons["CCG"])
    to_print += "%f %f %f %f // CAT CAC CAA CAG\n" % (codons["CAT"], codons["CAC"], codons["CAA"], codons["CAG"])
    to_print += "%f %f %f %f // CGT CGC CGA CGG\n" % (codons["CGT"], codons["CGC"], codons["CGA"], codons["CGG"])
    to_print += "\n"

    to_print += "%f %f %f %f // ATT ATC ATA ATG\n" % (codons["ATT"], codons["ATC"], codons["ATA"], codons["ATG"])
    to_print += "%f %f %f %f // ACT ACC ACA ACG\n" % (codons["ACT"], codons["ACC"], codons["ACA"], codons["ACG"])
    to_print += "%f %f %f %f // AAT AAC AAA AAG\n" % (codons["AAT"], codons["AAC"], codons["AAA"], codons["AAG"])
    to_print += "%f %f %f %f // AGT AGC AGA AGG\n" % (codons["AGT"], codons["AGC"], codons["AGA"], codons["AGG"])
    to_print += "\n"

    to_print += "%f %f %f %f // GTT GTC GTA GTG\n" % (codons["GTT"], codons["GTC"], codons["GTA"], codons["GTG"])
    to_print += "%f %f %f %f // GCT GCC GCA GCG\n" % (codons["GCT"], codons["GCC"], codons["GCA"], codons["GCG"])
    to_print += "%f %f %f %f // GAT GAC GAA GAG\n" % (codons["GAT"], codons["GAC"], codons["GAA"], codons["GAG"])
    to_print += "%f %f %f %f // GGT GGC GGA GGG\n" % (codons["GGT"], codons["GGC"], codons["GGA"], codons["GGG"])

    print(to_print)

    return codons, to_print
