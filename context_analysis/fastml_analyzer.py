#! /usr/local/python_anaconda/bin/python3.4

from optparse import OptionParser
import glob
from file_utilities import check_dirname, check_filename
import pandas as pd
import re
from scipy.stats import chi2_contingency, fisher_exact




def main():
    parser = OptionParser("usage: %prog [options]")
    parser.add_option("-d", "--dirname", dest="dirname", help="dirname that contains fastml result files")
    parser.add_option("-o", "--output", dest="output", help="output file")
    (options, args) = parser.parse_args()

    dirname = options.dirname
    dirname = check_dirname(dirname)

    output = options.output
    if output == None:
        output = dirname + "/fastml_analysis_output.csv"
    output = check_filename(output, Truefile=False)


    files = glob.glob(dirname + "/*/*.fasta")
    files = []
    if files == []:
        files = glob.glob(dirname + "/*seq.joint.txt")
    basenames = [f.split(".")[0] for f in files]

    df = pd.DataFrame(columns = ["family", "group", "mutation", "mutation_count_in_context", "context_count_overall",
                                 "mutation_count_overall"])
    
    for basename in basenames:
        df = run_on_basename(basename, df)

 
    #df.to_csv(output)


def create_mutation_dic(s = None):
    mut_dic = {}
    for i in ["A", "C", "G", "T"]:
        for j in ["A", "C", "G", "T"]:
            if i != j:
                if s != None and i == s:
                    mut_dic[i + j] = 0
                elif s == None:
                    mut_dic[i + j] = 0
    return mut_dic


def create_context_mutation_dic():
    context_dic = {}
    context_mutation_dic = {}
    for i in ["A", "C", "G", "T"]:
        for j in ["A", "C", "G", "T"]:
            for h in ["A", "C", "G", "T"]:
                context_dic[i + j + h] = 0
                context_mutation_dic[i + j + h] = create_mutation_dic(j)
    return context_dic, context_mutation_dic



def run_on_basename(basename, df, cutoff = 0.7):
    prob_marginal = basename + ".prob.marginal.txt"
    seq_marginal = basename + ".seq.marginal.txt"
    tree_ancestor = basename + ".tree.ancestor.txt"

    try:
        prob_marginal = check_filename(prob_marginal)
        seq_marginal = check_filename(seq_marginal)
        tree_ancestor = check_filename(tree_ancestor)
    except:
        return

    # remove position which have low probabilities
    positions_to_remove = []
    pm = open(prob_marginal, "r").read()
    poss = pm.split("\n\nmarginal probabilities at ")[1:]
    p = re.compile(": p\([GATC]\)=([01].\d*)") #gets probabilities under 0.7
    for pos in poss:
        r = p.findall(pos)
        r = [float(i) for i in r]
        pos_num = pos.split("\n")[0].split("position: ")[1]
        if min(r) < cutoff:
            positions_to_remove.append(int(pos_num)-1)


    alignment_len = int(pos_num) -1
    if float(len(positions_to_remove)) / alignment_len >= 0.3:
        print("warning - more than 30% of positions are removed")
    #print("removed %i positions (out of %i)" % (len(positions_to_remove), alignment_len))

    #create ancestor info - for each son - who is the father
    pattern = re.compile(r'\s+')
    ancestors =  open(tree_ancestor, "r").read().split("\n")[2:-2]
    ancestor_info = {}
    for line in ancestors:
        line = re.sub(pattern, '$', line)
        line = line.split("$")
        son = line[0]
        father = line[1]
        ancestor_info[son] = father

    #create seq dictionary - name and sequance
    seqs = open(seq_marginal, "r").read()
    seqs = seqs.split("\n>")[1:]
    seqs = {seq.split("\n")[0]:seq.split("\n")[1] for seq in seqs}




    #CB - close branches
    #FB - far branches
    #AB - all branches
    mutation_count_AB = create_mutation_dic()
    mutation_count_CB = create_mutation_dic()
    mutation_count_FB = create_mutation_dic()
    context_count_AB, context_mutation_count_AB = create_context_mutation_dic()
    context_count_CB, context_mutation_count_CB = create_context_mutation_dic()
    context_count_FB, context_mutation_count_FB = create_context_mutation_dic()


    GA_VS_NOT = pd.DataFrame([[0, 0], [0,0]], index=["GA", "other"], columns=["Apobec_context", "not_apobec_context"])
    GA_that_happand_vs_didnt_happen = pd.DataFrame([[0, 0], [0,0]], index=["Happand", "Didnt_Happen"], columns=["Apobec_context", "not_apobec_context"])


    for son in ancestor_info:
        father = ancestor_info[son]
        if father == "root!":
            continue
        son_seq = seqs[son]
        father_seq = seqs[father]
        pattern = re.compile("^N\d*") #checks if the node is an internal node
        if pattern.findall(son) == []:
            node_type = "external"
        else:
            node_type = "internal"

        # add context information
        for c in context_count_AB.keys():
            p = re.compile(c)
            l = len(p.findall(father_seq))
            context_count_AB[c] += l
            if node_type == "external": #if external node - goes to close branches - CB
                context_count_CB[c] += l
            else: # if internal node - goes to far brances - FB
                context_count_FB[c] += l


        #go over all positions
        for i in range(len(father_seq)-1):
            if i in positions_to_remove:
                continue
            if father_seq[i] != "G":
                continue
            if son_seq[i] == "G":
                if father_seq[i+1] in ["G", "A"]:
                    GA_that_happand_vs_didnt_happen.Apobec_context.Didnt_Happen += 1
                else:
                    GA_that_happand_vs_didnt_happen.not_apobec_context.Didnt_Happen += 1
            else:
                if father_seq[i+1] in ["G", "A"]:
                    GA_that_happand_vs_didnt_happen.Apobec_context.Happand += 1
                else:
                    GA_that_happand_vs_didnt_happen.not_apobec_context.Happand += 1

        oddsratio, pvalueisher_exact = fisher_exact(GA_that_happand_vs_didnt_happen)
        percentage = GA_that_happand_vs_didnt_happen.groupby(level=0).apply(lambda x:
                                                      100 * x / float(x.sum(axis=1)))
    if pvalueisher_exact <= 0.05:
        print(basename, pvalueisher_exact)
        print(percentage)



    else:
        print(basename, pvalueisher_exact)
        print(percentage)
        """
        different_positions = [i for i in range(len(son_seq)) if son_seq[i] != father_seq[i]]
        for diff_pos in different_positions:
            if diff_pos in positions_to_remove:
                continue # position that are supposed to be removed - didn't infer their ancestry well enough
            if alignment_len - diff_pos == 0 or diff_pos == 0:
                continue # position at the beginning or end of the sequence
            if son_seq[diff_pos] == "-":
                continue # positions that are gapped in the son sequence - non relevant mutations
            mutation = father_seq[diff_pos] + son_seq[diff_pos]
            context = father_seq[diff_pos-1] + father_seq[diff_pos] +  father_seq[diff_pos+1] # from ancestral sequence
            if mutation not in mutation_count_AB.keys():
                continue # mutation is not one of the twelve
            if context not in context_count_AB.keys():
                continue # context is not one of the 16
            mutation_count_AB[mutation] += 1
            context_mutation_count_AB[context][mutation] += 1
            if node_type == "external": #if the node is external  - goes to close branches - CB
                mutation_count_CB[mutation] += 1
                context_mutation_count_CB[context][mutation] += 1
            else: # if hte node is internal - goes to far branches - FB
                mutation_count_FB[mutation] += 1
                context_mutation_count_FB[context][mutation] += 1

            if mutation == "GA":
                if father_seq[diff_pos+1] in ["G", "A"]:
                    GA_VS_NOT.Apobec_context.GA += 1
                else:
                    GA_VS_NOT.not_apobec_context.GA += 1
            else:
                if father_seq[diff_pos+1] in ["G", "A"]:
                    GA_VS_NOT.Apobec_context.other += 1
                else:
                    GA_VS_NOT.not_apobec_context.other += 1



    oddsratio, pvalueisher_exact = fisher_exact(GA_VS_NOT)
    percentage = GA_VS_NOT.groupby(level=0).apply(lambda x:
                                                   100 * x / float(x.sum(axis=1)))
    if pvalueisher_exact <= 0.05:
        print(basename, pvalueisher_exact)
        print(percentage)
    else:
        print(basename)
    """

    df = pd.DataFrame(columns=["Close_branches", "Far_branches"])

    for c in context_mutation_count_AB.keys():
        for m in context_mutation_count_AB[c].keys():
            # far mutation
            mutation_count_in_context_FB = context_mutation_count_FB[c][m]
            context_count_overall_FB = context_count_FB[c]

            mutation_count_in_context_CB = context_mutation_count_CB[c][m]
            context_count_overall_CB = context_count_CB[c]

            if m  == "GA":
                if mutation_count_in_context_FB == 0 and mutation_count_in_context_CB == 0:
                    continue

                temp_df = pd.DataFrame({'Close_branches':[mutation_count_in_context_CB],
                    'Far_branches':[mutation_count_in_context_FB],},
                   index = [c])
                df = pd.concat([df,temp_df])




def statistical_tests():
    g, p, dof, expctd = chi2_contingency(df.transpose(), lambda_="log-likelihood")

    df = df.transpose()

    percentage = df

    percentage = percentage.groupby(level=0).apply(lambda x:
                                                     100 * x / float(x.sum(axis=1)))
    #print(basename, p)


    #if p < 0.05:
    #    print(percentage)

    """
    print(basename, p)
    print(df)
    print(df.sum())
    """
    #df["CB_proportion"] = df.Close_branches / (df.Close_branches + df.Far_branches)
    #df["FB_proportion"] = df.Far_branches / (df.Close_branches + df.Far_branches)

    #if p < 0.05:
    #    print(percentage)

    APOBEC_columns = [col for col in list(percentage) if col.endswith(('A', "G"))]
    NOT_APOBEC_columns = [col for col in list(percentage) if col.endswith(('C', "T"))]


    apobec_df = pd.DataFrame()
    apobec_df["APOBEC_associated"] = df[APOBEC_columns].sum(axis=1)
    apobec_df["not_APOBEC_associated"] = df[NOT_APOBEC_columns].sum(axis=1)
    #print(apobec_df)
    g, p, dof, expctd = chi2_contingency(apobec_df, lambda_="log-likelihood")

    oddsratio, pvalueisher_exact = fisher_exact(apobec_df)
    percentage = apobec_df.groupby(level=0).apply(lambda x:
                                                   100 * x / float(x.sum(axis=1)))
    if pvalueisher_exact <= 0.05:
        print(basename, pvalueisher_exact)
        print(percentage)
    return

    percentage = apobec_df.groupby(level=0).apply(lambda x:
                                                   100 * x / float(x.sum(axis=1)))
    print(basename, p)
    if p < 0.05:
        print(percentage)



    """
    for m in mutations:
        for c in mutations[m]:
            if mutations[m][c] == 0:
                continue
            df = df.append({"family":family, "group":group, "mutation":m, "context":c,
                            "mutation_count_in_context":mutations[m][c], "context_count_overall":context_count[c],
                            "mutation_count_overall":mutation_count[m], "close_mutation":False},
                                   ignore_index = True)
    for m in mutations_close:
        for c in mutations_close[m]:
            if mutations_close[m][c] == 0:
                continue
            df = df.append({"family":family, "group":group, "mutation":m, "context":c,
                            "mutation_count_in_context":mutations_close[m][c],
                            "context_count_overall":context_count_close[c],
                            "mutation_count_overall":mutation_count_close[m], "close_mutation":True},
                                   ignore_index = True)
    """
    return df

    """
    
        mutations_all = {}
    mutations_close_branches = {}
    mutations_far_branches = {}
    for mutation in mutation_count_all:
        mutations_all[mutation] = {}
        mutations_close_branches[mutation] = {}
        mutations_far_branches[mutation] = {}
        for context in context_count_all:
            mutations_all[mutation][context] = 0
            mutations_close_branches[mutation][context] = 0
            mutations_far_branches[mutation][context] = 0
            """



if __name__ == "__main__":
    main()
