#! /usr/local/python_anaconda/bin/python3.4


from Bio import SeqIO
from Bio import AlignIO
from Bio import Alphabet
from Bio.Alphabet import IUPAC
import collections
from file_utilities import check_filename
import pandas as pd
import textwrap
from collections import Counter
from itertools import product
from phyVirus.get_baltimore import get_baltimore_classifiaction



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


def get_major_and_minor_consensus(aln_file, in_format="fasta"):
    """
    calculates major and minor consensus and each position's probability
    - major consensus - the most prominent base (including "-")
    - minor consensus - the most prominent base (not including "-")
    :param aln_file: alignment file path
    :param in_format: input alignment format (default: fasta)
    :return: major_consensus, major_freqs, minor_consensus, minor_freqs
    """
    aln_file = check_filename(aln_file)
    aln = AlignIO.read(aln_file, in_format, alphabet=Alphabet.Gapped(IUPAC.unambiguous_dna))
    len_aln = len(aln[0])
    num_of_seq = len(aln)
    major_consensus = ""
    major_freqs = []
    minor_consensus = ""
    minor_freqs = []
    for i in range(len_aln):
        counter = collections.Counter(aln[:, i])
        major_count = 0
        minor_count = 0
        major_char = ""
        minor_char = ""
        for j in counter:
            if counter[j] > major_count:
                major_count = counter[j]
                major_char = j
                if j != "-":
                    minor_count = counter[j]
                    minor_char = j
            if counter[j] > minor_count and j != "-":
                if j not in ["A", "C", "G", "T"]:
                    minor_count = counter[j]
                    minor_char = "N"
                else:
                    minor_count = counter[j]
                    minor_char = j
        gap_count = counter["-"]
        major_consensus += major_char
        major_freqs.append(round(major_count / (num_of_seq - gap_count), 2))

        minor_consensus += minor_char
        minor_freqs.append(round(minor_count / (num_of_seq - gap_count), 2))

    return major_consensus, major_freqs, minor_consensus, minor_freqs


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

def get_amino_acid_freqs(filename, no_strange_aas = True):
    """
    gets a fasta file of protein seqs and returns a dataframe of amino acid frequencies for each fasta entry
    id - ncbi id
    :param filename: file path of amino acid fasta file
    :param no_strange_aas: if to ignore strange aa letters
    :return: dataframe with aa frequecnies
    """
    filename = check_filename(filename)
    aa = open(filename, "r").read()
    aa = aa.split(">")[1:]
    aa_data = pd.DataFrame()
    for a in (aa):
        ncbi_id = a.split("|")[3]
        seq = "".join(a.split("\n")[1:])
        length = len(seq)
        if no_strange_aas:
            # remove problematic AA's
            if "X" in seq:
                seq = "".join(seq.split("X"))
            if "J" in seq:
                seq = "".join(seq.split("J"))
            if "Z" in seq:
                seq = "".join(seq.split("Z"))
            if "B" in seq:
                seq = "".join(seq.split("B"))
        counter = dict(Counter(seq))
        for i in counter:
            counter[i] = counter[i] / float(length)
        counter["ncbi_id"] = ncbi_id
        aa_data = aa_data.append(counter, ignore_index=True)

    return aa_data


def get_codon_freqs(filename, no_strange_nucs = True):
    """
    gets a fasta file of nucleotide seqs and returns a dataframe of codon frequencies for each fasta entry
    id - ncbi id
    :param filename: file path of amino acid fasta file
    :param no_strange_nucs: if to ignore strange nucleotide letters
    :return: dataframe with codon frequecnies
    """
    alphabet = ["a", "c", "t", "g"]
    codons_possibilities = ["".join(i) for i in product(alphabet, repeat=3)]
    
    filename = check_filename(filename)
    gb = open(filename, "r").read()
    items = gb.split(">")[1:]
    codon_info = pd.DataFrame()
    for k in (items):
        result = {}
        ncbi_id = k.split("|")[3]
        seq = "".join(k.split("\n")[1:])
        if len(seq) % 3 != 0:
            continue
        codons = textwrap.wrap(seq, 3)
        counter = {}
        for c in codons:
            if c in codons_possibilities:
                if not c in counter.keys():
                    counter[c] = 0
                counter[c] += 1
        sum_co = sum(counter.values())
        for codon in counter:
            counter[codon] = counter[codon] / float(sum_co)
        counter['ncbi_id'] = ncbi_id
        codon_info = codon_info.append(counter, ignore_index=True)
    return codon_info



def get_dinucleotide_odds_ratio(fasta_file, in_format="fasta", output_dir=None):
    fasta_file = check_filename(fasta_file)

    dinucs_or = pd.DataFrame()
    dinucs = {}
    for p in product("ACTG", repeat=2):
        dinucs[p[0]+p[1]] = 0
    nucs = {"A":0, "C":0, "G":0, "T":0}
    comp = {"A":"T", "C":"G", "G":"C", "T":"A"}

    basename = fasta_file.split("/")[-1].split(".fasta")[0].split(".aln")[0].split(".aln.best.fas")[0].split(".codon_aln.best.fas")[0]
    family = basename.split("_")[0]
    baltimore = get_baltimore_classifiaction(family)
    if output_dir == None:
        output_base = fasta_file.split(".fasta")[0].split(".aln")[0].split(".aln.best.fas")[0].split(".codon_aln.best.fas")[0]
    else:
        output_base = "%s/%s" % (output_dir, basename)
    output_dinuc = output_base + ".dinuc_odds_ratio"
    output_dinuc_averaged = output_base + ".dinuc_averaged_odds_ratio"

    seqs = list(SeqIO.parse(fasta_file, format=in_format))
    for s in seqs:
        for i in dinucs:
            dinucs[i] = 0
        for i in nucs:
            nucs[i] = 0

        s.seq = str(s.seq).replace("-", "").upper()

        #count and calculate nucleotide freqs
        for i in nucs:
            nucs[i] = s.seq.count(i)
        count = len(s.seq)
        for i in nucs:
            nucs[i] = nucs[i] / count
        #count and calculate dinucleotide freqs
        for i in dinucs:
            dinucs[i] = s.seq.count(i)
        count_dinucs = sum(dinucs.values())
        for i in dinucs:
            dinucs[i] = dinucs[i] / count_dinucs
        #calculate odds ratio
        for i in dinucs:
            if "ds" in baltimore:
                comp_dinuc = comp[i[1]] + comp[i[0]]
                dinucs[i] = (2*(dinucs[i] + dinucs[comp_dinuc]) /
                             ((nucs[i[0]] + nucs[i[1]]) * (nucs[comp_dinuc[0]] + nucs[comp_dinuc[1]])))
            else:
                dinucs[i] = dinucs[i] / (nucs[i[0]] * nucs[i[1]])

        dinucs_or = dinucs_or.append(
            {"baltimore":baltimore, "family": family, "basename": basename, "seq_name": s.id,
             **dinucs},
            ignore_index=True)


    dinucs_average = dinucs_or.mean(axis=0).to_frame().transpose()
    dinucs_average = pd.concat([dinucs_average,
                                pd.DataFrame([{"baltimore":baltimore, "family": family, "basename": basename,}])], axis=1)


    dinucs_or.to_csv(output_dinuc, index=False)
    dinucs_average.to_csv(output_dinuc_averaged, index=False)



def analyze_nuc_frequencies_and_wobble_freqs(fasta_file, in_format="fasta", output=None):
    fasta_file = check_filename(fasta_file)
    if output == None:
        output = fasta_file.split(".fasta")[0] + ".base_freqs.csv"
    else:
        output = check_filename(output, Truefile=False)


    df = pd.DataFrame(columns=["filename", "dir", "base_file", "A", "C", "T", "G", "wob_A", "wob_C", "wob_T", "wob_G",
                               "non_wob_A", "non_wob_G", "non_wob_C", "non_wob_T"])
    seqs = list(SeqIO.parse(fasta_file, format=in_format))
    base_freqs = {"A": 0, "G": 0, "C": 0, "T": 0}
    wobble_freqs = {"wob_A": 0, "wob_G": 0, "wob_C": 0, "wob_T": 0}
    non_wobble_freqs = {"non_wob_A":0, "non_wob_G":0, "non_wob_C":0, "non_wob_T":0}
    count = 0
    wobble_count = 0
    non_wobble_count = 0
    for s in seqs:
        if len(s.seq) % 3 != 0:
            continue
        #general freqs
        s.seq = s.seq.upper()
        a = s.seq.count("A")
        c = s.seq.count("C")
        t = s.seq.count("T")
        g = s.seq.count("G")
        count += len(s.seq)
        base_freqs["A"] += a
        base_freqs["C"] += c
        base_freqs["T"] += t
        base_freqs["G"] += g
        #wobble freqs
        wobble_s = s.seq[2::3]

        a = wobble_s.count("A")
        c = wobble_s.count("G")
        t = wobble_s.count("C")
        g = wobble_s.count("T")
        wobble_count += len(wobble_s)
        wobble_freqs["wob_A"] += a
        wobble_freqs["wob_G"] += c
        wobble_freqs["wob_C"] += t
        wobble_freqs["wob_T"] += g
        #non wobble freqs
        non_wobble_s = s.seq[0::3] + s.seq[1::3]
        a = non_wobble_s.count("A")
        c = non_wobble_s.count("G")
        t = non_wobble_s.count("C")
        g = non_wobble_s.count("T")
        non_wobble_count += len(non_wobble_s)
        non_wobble_freqs["non_wob_A"] += a
        non_wobble_freqs["non_wob_G"] += c
        non_wobble_freqs["non_wob_C"] += t
        non_wobble_freqs["non_wob_T"] += g

    for k in base_freqs.keys():
        base_freqs[k] = base_freqs[k] / float(count)

    for k in wobble_freqs.keys():
        wobble_freqs[k] = wobble_freqs[k] / float(wobble_count)

    for k in non_wobble_freqs.keys():
        non_wobble_freqs[k] = non_wobble_freqs[k] / float(non_wobble_count)

    df = df.append({"filename":fasta_file, "dir":fasta_file.split("/")[-2], "base_file":fasta_file.split("/")[-1], **base_freqs, **wobble_freqs, **non_wobble_freqs}, ignore_index=True)
    df.to_csv(output)


def analyze_nuc_frequencies_and_wobble_freqs_from_aln(aln_file, in_format="fasta", output_dir=None):
    aln_file = check_filename(aln_file)
    base = aln_file.split(".fasta")[0].split(".aln")[0].split(".aln.best.fas")[0].split(".codon_aln.best.fas")[0]
    basename = base.split("/")[-1]
    if output_dir == None:
        output_freqs = base + ".base_freqs_info.csv"
        output_counts = base + ".base_counts_info.csv"
        output_averaged_freqs = base + ".base_freqs_averaged_freqs.csv"
        output_averaged_counts = base + ".base_freqs_averaged_counts.csv"
    else:
        output_freqs = output_dir + basename + ".base_freqs_info.csv"
        output_counts = output_dir + basename + ".base_counts_info.csv"
        output_averaged_freqs = output_dir + basename  + ".base_freqs_averaged_freqs.csv"
        output_averaged_counts = output_dir + basename + ".base_freqs_averaged_counts.csv"

    counts = pd.DataFrame()
    freqs = pd.DataFrame()
    aln = AlignIO.read(aln_file, format=in_format)



    family = basename.split("_")[0]
    baltimore = get_baltimore_classifiaction(family)

    for a in aln:
        a.seq = a.seq.upper()
        wobble = a[2::3]
        non_wobble = a[0:3] + a[1::3]

        all_count = a.seq.count("A") + a.seq.count("C") + a.seq.count("G") + a.seq.count("T")
        wobble_count = wobble.seq.count("A") + wobble.seq.count("C") + wobble.seq.count("G") + wobble.seq.count("T")
        non_wobble_count = non_wobble.seq.count("A") + non_wobble.seq.count("C") + non_wobble.seq.count("G") + non_wobble.seq.count("T")


        counts = counts.append({"baltimore":baltimore, "family":family, "basename":basename, "seq_name":a.id,
                        "A": a.seq.count("A"), "C": a.seq.count("C"), "G": a.seq.count("G"), "T": a.seq.count("T"),
                        "wob_A": wobble.seq.count("A"), "wob_C": wobble.seq.count("C"), "wob_C": wobble.seq.count("C"),
                        "non_wob_A":non_wobble.seq.count("A"), "non_wob_C":non_wobble.seq.count("C"), "non_wob_G":non_wobble.seq.count("G"), "non_wob_T":non_wobble.seq.count("T")},
                               ignore_index=True)

        freqs = freqs.append({"baltimore":baltimore, "family":family, "basename":basename, "seq_name":a.id,
                    "A": a.seq.count("A")/all_count, "C": a.seq.count("C")/all_count,
                           "G": a.seq.count("G")/all_count, "T": a.seq.count("T")/all_count,
                    "wob_A": wobble.seq.count("A")/wobble_count, "wob_C": wobble.seq.count("C")/wobble_count,
                           "wob_G": wobble.seq.count("G")/wobble_count, "wob_T": wobble.seq.count("T")/wobble_count,
                    "non_wob_A":non_wobble.seq.count("A")/non_wobble_count, "non_wob_C":non_wobble.seq.count("C")/non_wobble_count,
                           "non_wob_G":non_wobble.seq.count("G")/non_wobble_count, "non_wob_T":non_wobble.seq.count("T")/non_wobble_count},
                           ignore_index=True)

    averaged_freqs = freqs.mean(axis=0).to_frame().transpose()
    averaged_freqs = pd.concat([averaged_freqs,
                                pd.DataFrame([{"baltimore":baltimore, "family":family, "basename":basename}])],
                                axis=1)
    averaged_counts = counts.mean(axis=0).to_frame().transpose()
    averaged_counts = pd.concat([averaged_counts,
                                pd.DataFrame(
                                    [{"baltimore":baltimore, "family":family, "basename":basename}])],
                               axis=1)

    counts.to_csv(output_counts, index=False)
    freqs.to_csv(output_freqs, index=False)
    averaged_freqs.to_csv(output_averaged_freqs, index=False)
    averaged_counts.to_csv(output_averaged_counts, index=False)
