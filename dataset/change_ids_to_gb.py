#! /usr/local/python_anaconda/bin/python3.4

from optparse import OptionParser
from file_utilities import check_dirname, check_filename, change_filename
import glob
from seqFileAnalyzer import get_consensus_percentage
import pandas as pd
from PAML_utilities import write_ctl_file
from Bio import SeqIO
from collections import Counter

def main():
    cross_files = glob.glob("/sternadi/nobackup/volume1/phyVirus/*txt")
    for cross in cross_files:
        family = cross.split(".txt")[0].split("/")[-1]
        cross = open(cross, "r").readlines()
        cross_data = {l.split(" ")[0]:l.split("|")[0].split(" ")[1] for l in cross}

        files = glob.glob("/sternadi/nobackup/volume1/phyVirus/fasta/%s*" % family)
        for f in files:
            fasta_data = open(f, "r").readlines()
            new_fasta_data = ""
            new_fasta_path = "/new_fasta/".join(f.split("/fasta/"))
            for line in fasta_data:
                if ">" not in line:
                    new_fasta_data += line
                else:
                    id = line.split("\n")[0]
                    new_id = cross_data[id]
                    new_fasta_data += ">%s\n" %new_id

            temp = open(new_fasta_path, "w")
            temp.write(new_fasta_data)
            temp.close()




def change_name_to_gb():
    files = glob.glob("/sternadi/nobackup/volume1/influenza/*fasta")
    for f in files:
        fasta = list(SeqIO.parse(f, "fasta"))
        new_fasta = []
        for seq in fasta:
            seq.id = seq.id.split("|")[0]
            seq.name = ""
            seq.description = ""
            new_fasta.append(seq)
        SeqIO.write(new_fasta, f, "fasta")



def change_dup_names():
    files = glob.glob("/sternadi/nobackup/volume1/phyVirus/new_fasta/*fasta")
    for f in files:
        fasta = list(SeqIO.parse(f, "fasta"))
        ids = [seq.id for seq in fasta]
        counts = Counter(ids)
        dups = {}
        for id in ids:
            if counts[id] > 1:
                dups[id] = 0
        if dups != {}:
            new_fasta = []
            for seq in fasta:
                if seq.id in dups:
                    dups[seq.id] = dups[seq.id] + 1
                    seq.id = seq.id + "_" + str(dups[seq.id])
                    seq.description = seq.id
                new_fasta.append(seq)
            new_file = "fasta_no_dup".join(f.split("new_fasta"))
            SeqIO.write(new_fasta, new_file, "fasta")



if __name__ == "__main__":
    main()


