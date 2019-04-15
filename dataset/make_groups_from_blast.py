#! /usr/local/python_anaconda/bin/python3.4

from optparse import OptionParser
from file_utilities import check_dirname, check_filename, change_filename
import glob
from seqFileAnalyzer import get_consensus_percentage
import pandas as pd
from PAML_utilities import write_ctl_file
from Bio import SeqIO


def main():
    base = "Picorna"
    print(base)
    blast_results = glob.glob("/sternadi/home/volume3/taliakustin/phyVirus2/big_files/splited/%s*txt" % base)
    fasta_file = glob.glob("/sternadi/home/volume3/taliakustin/phyVirus2/%s.fasta" % base)[0]
    outfiles_prefix = "/sternadi/home/volume3/taliakustin/phyVirus2/big_files/groups_flavi/%s_" % base
    groups = {}
    evalue_cutoff = 0.000000000001
    group_counter = 1
    for f in blast_results:
        with open(f, "r") as file:
            data = file.readlines()
            res = []
            for l in data:
                eval = float(l.split("\t")[10])
                if eval <= evalue_cutoff or eval == 0:
                    res.append(l.split("\t")[1])
            exist = False
            selected = []
            count_r = 0


            while count_r < len(res):
                for g in groups:
                    if res[count_r] in groups[g]:
                        exist = True
                        selected.append(g)
                count_r += 1
                
            selected = list(set(selected))
            if exist:
                if len(selected) > 1:
                    print(selected)
                    all_selected = []
                    for i in selected:
                        all_selected = all_selected + groups[i]
                    groups[min(selected)] = all_selected
                    for i in selected:
                        if i == min(selected):
                            continue
                        else:
                            groups[i] = []
                    selected = min(selected)
                elif len(selected) == 1:
                    selected = selected[0]

                for r in res:
                    if r in groups[selected]:
                        continue
                    else:
                        groups[selected].append(r)
            else:
                groups[group_counter] = []
                for r in res:
                    groups[group_counter].append(r)
                print("group+1 = %i" % group_counter)
                group_counter += 1

    for g in groups:
        print(g, len(groups[g]))


    seqs = list(SeqIO.parse(fasta_file, format="fasta"))
    groups_proteins = {}

    for s in seqs:
        for g in groups:
            if s.description in groups[g]:
                if g not in groups_proteins.keys():
                    groups_proteins[g] = []
                groups_proteins[g].append(s)

    for g in groups_proteins:
        print(g, len(groups_proteins[g]))

    for g in groups_proteins:
        if len(groups_proteins[g]) >= 10:
            SeqIO.write(groups_proteins[g], outfiles_prefix + str(g) + ".fasta", format="fasta")



if __name__ == "__main__":
    main()
