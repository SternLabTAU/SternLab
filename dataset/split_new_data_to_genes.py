#! /usr/local/python_anaconda/bin/python3.4

from optparse import OptionParser
from file_utilities import check_dirname, check_filename, change_filename
import glob
from seqFileAnalyzer import get_consensus_percentage
import pandas as pd
from PAML_utilities import write_ctl_file

def gene_split():
    for f in files:
        fasta = list(SeqIO.parse(f, "fasta"))
        groups = {}
        for f in fasta:
            if "Gene Symbol" not in f.description:
                gene = "unknown"
            else:
                gene = f.description.split("Gene Symbol:")[1].split("|")[0]
            if gene not in groups.keys():
                groups[gene] = []
            groups[gene].append(f)

        for i in groups.keys():
            filename = "/sternadi/nobackup/volume1/influenza/Orthomyxo_infA_%s.fasta" % i
            SeqIO.write(groups[i], filename, "fasta")



def segment_split():
    files = glob.glob("/sternadi/nobackup/volume1/influenza/*unknown*")
    for f in files:
        inf_type = f.split("_inf")[1].split("_")[0]
        fasta = list(SeqIO.parse(f, "fasta"))
        groups = {}
        for f in fasta:
            if "Segment" not in f.description:
                gene = "unknown"
            else:
                gene = f.description.split("Segment: ")[1].split("|")[0]
            if gene not in groups.keys():
                groups[gene] = []
            groups[gene].append(f)

        for i in groups.keys():
            filename = "/sternadi/nobackup/volume1/influenza/Orthomyxo_inf%i_%s.fasta" % (inf_type, i)
            SeqIO.write(groups[i], filename, "fasta")





if __name__ == "__main__":
    main()


