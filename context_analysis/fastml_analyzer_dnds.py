#! /usr/local/python_anaconda/bin/python3.4

from optparse import OptionParser
import glob
from file_utilities import check_dirname, check_filename
import pandas as pd
import re
from scipy.stats import chi2_contingency, fisher_exact
from Bio.Seq import Seq
from Bio.Seq import Seq





def main():
    parser = OptionParser("usage: %prog [options]")
    parser.add_option("-d", "--dir", dest="directory", help="directory with fastml output files")

    (options, args) = parser.parse_args()
    dir = options.dir
    dir = check_dirname(dir)
    files = glob.glob(dir + "*prob.marginal.txt")
    basenames = [f.split(".prob.marginal.txt")[0] for f in files]
    df = pd.DataFrame(columns=["Basename", "Mutation", "Branch", "Context", "Mutation_type", "Codon_position", "APOBEC_context"])

    for basename in basenames:
        print(basename)
        prob_marginal = basename + ".prob.marginal.txt"
        seq_marginal = basename + ".seq.marginal.txt"
        tree_ancestor = basename + ".tree.ancestor.txt"

        positions_to_remove = get_position_to_remove(prob_marginal)
        ancestor_info, seqs = get_sequence_and_ancestry_data(tree_ancestor, seq_marginal)
        df = go_over_positions(ancestor_info, seqs, positions_to_remove, basename, df)



    df.to_csv(dir + "fastml_anlysis.csv", index=False)


def go_over_positions(ancestor_info, seqs, positions_to_remove, basename, df, mutations_to_check=["GA"]):
    same_positions = 0
    bad_positions = 0
    positions_with_one_difference = 0
    positions_with_more_than_one_difference = 0

    for son in ancestor_info:
        father = ancestor_info[son]
        if father == "root!":
            continue
        son_seq = seqs[son]
        father_seq = seqs[father]
        pattern = re.compile("^N\d*") #checks if the node is an internal node
        if pattern.findall(son) == []:
            branch = "external"
        else:
            branch = "internal"


        for i in range(0, len(father_seq) - 1):
            if i in positions_to_remove:
                bad_positions += 1
                continue
            context = father_seq[i-1:i+2]

            father_nuc = father_seq[i]
            son_nuc = son_seq[i]

            if son_nuc == father_nuc:
                same_positions += 1
                continue
            if son_nuc not in ["A", "C", "T", "G"]:
                continue
            elif "-" in son_nuc:
                continue
            if i == 0 or i == len(father_seq) - 2:
                continue
            codon_position = i % 3
            if codon_position == 0:
                father_codon = father_seq[i:i+3]
                son_codon = son_seq[i:i+3]
            elif codon_position == 1:
                father_codon = father_seq[i-1:i+2]
                son_codon = son_seq[i-1:i+2]
            else:
                father_codon = father_seq[i-2:i+1]
                son_codon = son_seq[i-2:i+1]
            if "-" in son_codon:
                continue
            changes = [j for j in range(len(father_codon)) if father_codon[j] != son_codon[j]]
            number_of_changes = len(changes)
            mutation = father_nuc + son_nuc
            father_aa = str(Seq(father_codon).translate())
            son_aa = str(Seq(son_codon).translate())
            if father_aa == "*":
                continue
            if number_of_changes == 1:
                positions_with_one_difference += 1

                if father_aa == son_aa:
                    mutation_type = "synonymous"
                    if i%3 == 1:
                        print(i, father_codon, son_codon)
                else:
                    mutation_type = "non-synonymous"

            elif number_of_changes >= 2:
                positions_with_more_than_one_difference +=1
                if i%3 == 0:
                    mutation_type = "non-synonymous"
                elif i%3 == 1:
                    mutation_type = "non-synonymous"
                elif i%3 == 2:
                    mutation_type = "synonymous"

            #print(i, context)
            if  context[2] in ["G", "A"]:
                APOBEC_context = True
            else:

                APOBEC_context =  False



            df = df.append({"Basename": basename, "Mutation": mutation, "Branch": branch, "Context": context,
                                "Mutation_type": mutation_type, "Codon_position": (i % 3) + 1, "APOBEC_context":APOBEC_context}, ignore_index=True)

    return(df)







def get_sequence_and_ancestry_data(tree_ancestor, seq_marginal):
    # create ancestor info - for each son - who is the father
    pattern = re.compile(r'\s+')
    ancestors = open(tree_ancestor, "r").read().split("\n")[2:-2]
    ancestor_info = {}
    for line in ancestors:
        line = re.sub(pattern, '$', line)
        line = line.split("$")
        if len(line[0].split("N")) != 2 and line[0].count("N") == 2:
            line = line[0].split("N")
            son = "N" +  line[1]
            father = "N" + line[2]
        else:
            son = line[0]
            father = line[1]
        ancestor_info[son] = father

    #create seq dictionary - name and sequance
    seqs = open(seq_marginal, "r").read()
    seqs = seqs.split("\n>")[1:]
    seqs = {seq.split("\n")[0]:seq.split("\n")[1] for seq in seqs}

    return(ancestor_info, seqs)


def get_position_to_remove(prob_marginal, cutoff = 0.7):
    positions_to_remove = []
    pm = open(prob_marginal, "r").read()
    poss = pm.split("\n\nmarginal probabilities at ")[1:]
    p = re.compile(": p\([GATC]\)=([01].\d*)")
    for pos in poss:
        r = p.findall(pos)
        r = [float(i) for i in r]
        pos_num = pos.split("\n")[0].split("position: ")[1]
        if min(r) < cutoff:
            positions_to_remove.append(int(pos_num) - 1)
    alignment_len = int(pos_num) -1
    print(alignment_len)
    if float(len(positions_to_remove)) / alignment_len >= 0.3:
        print("warning - more than 30% of positions are removed")
    return positions_to_remove



if __name__ == "__main__":
    main()
