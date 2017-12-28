import os
import shutil
import PAML_utilities
import file_utilities
import pbs_jobs
from optparse import OptionParser


def main():
    parser = OptionParser("usage: %prog[options]")
    parser.add_option("-s", "--site_gene_folder", dest="GENES_FOLDER", help="The folder that contains all candida genes subfolders site results")
    parser.add_option("-b", "--branch", dest="BRANCH", help="The branch site model folder where the results for 6, 13 both and 13up will be created")
    parser.add_option("-c", "--branch_choice", dest="CHOICE", help="Choose the desired clade: 13, 6, 13up or both")

    (options, args) = parser.parse_args()

    GENES_FOLDER = options.GENES_FOLDER
    BRANCH = options.BRANCH
    CHOICE = options.CHOICE

    copy_phylips_to_branch_folder(GENES_FOLDER, BRANCH, CHOICE)
    run_branch_site_codeml(GENES_FOLDER, CHOICE)


def copy_phylips_to_branch_folder(all_genes_folder, branch_folder, branch_choice):
    """
    Creates a new gene folder in both 6 and 13 folders of "branch sites" folder
    :param all_genes_folder: No Ref folder where all the genes are seperated into subfolders
    :param branch_folder: The branch site folder where the gene folders will be created inside 6 or 13 subfolder
    :param branch_choice: Which branch subfolder will be created (13, 6, both or 13up)
    """
    for gene_name in os.listdir(all_genes_folder):
        gene_folder = all_genes_folder + gene_name + "/"
        if branch_choice == "6":
            dst_gene_folder_6 = branch_folder + "6/" + gene_name
            os.mkdir(dst_gene_folder_6)
            shutil.copy(gene_folder + gene_name + ".phy", dst_gene_folder_6 + "/" + gene_name + ".phy")
        elif branch_choice =="13":
            dst_gene_folder_13 = branch_folder + "13/" + gene_name
            os.mkdir(dst_gene_folder_13)
            shutil.copy(gene_folder + gene_name + ".phy", dst_gene_folder_13 + "/" + gene_name + ".phy")
        elif branch_choice == "13up":
            dst_gene_folder_13to17 = branch_folder + "13up/" + gene_name
            os.mkdir(dst_gene_folder_13to17)
            shutil.copy(gene_folder + gene_name + ".phy", dst_gene_folder_13to17 + "/" + gene_name + ".phy")
        else:
            dst_gene_folder_both = branch_folder + "both/" + gene_name
            os.mkdir(dst_gene_folder_both)
            shutil.copy(gene_folder + gene_name + ".phy", dst_gene_folder_both + "/" + gene_name + ".phy")



def run_branch_site_codeml(gene_folder, branch_num):
    """
    Creates CLT files for the null and alternative models for every gene in gene folder, according to the branch num.
    Then it runs them on the cluster.
    :param gene_folder: The Branch-site main folder that contains "6" and "13" subfolders
    :param branch_num: 6, 13 or 13up as a string
    """
    unmarked_tree_path = "/sternadi/home/volume2/ella/Candida/Trees/Modified/SNP_seq_without_ref.phy_phyml_tree.txt"
    if branch_num == "6":
        tree_path = "/sternadi/home/volume2/ella/Candida/Trees/Modified/SNP_seq_without_ref.phy_phyml_tree-6.txt"
    elif branch_num == "13":
        tree_path = "/sternadi/home/volume2/ella/Candida/Trees/Modified/SNP_seq_without_ref.phy_phyml_tree-13.txt"
    elif branch_num == "13up":
        tree_path = "/sternadi/home/volume2/ella/Candida/Trees/Modified/SNP_seq_without_ref.phy_phyml_tree-13to17.txt"
    else:
        tree_path = "/sternadi/home/volume2/ella/Candida/Trees/Modified/SNP_seq_without_ref.phy_phyml_tree-13and6.txt"
    for gene_name in os.listdir(gene_folder + "/" + branch_num):
        path_prefix = gene_folder + "/" + branch_num + "/" + gene_name + "/" + gene_name
        seq_path = path_prefix + ".phy"

        null_clt = path_prefix + "-n" + branch_num + ".CLT"
        null_output = path_prefix + "-n" + branch_num + "-Results.txt"

        alternative_clt = path_prefix + "-a" + branch_num + ".CLT"
        alternative_output = path_prefix + "-a" + branch_num + "-Results.txt"

        PAML_utilities.candida_write_ctl_codeml_clade_model(null_clt, seq_path, unmarked_tree_path, null_output, 1)
        PAML_utilities.candida_write_ctl_codeml_clade_model(alternative_clt, seq_path, tree_path, alternative_output, 0)

        candida_codeml_branch_runner(null_clt, alternative_clt, path_prefix, branch_num)


def candida_codeml_branch_runner(clt1, clt2, result_prefix, branch_num, alias ="cml"):
    """
    run codeml program from PAML on cluster - (runs both alternative and null model in one job).
    :param ctl: ctl file path
    :param alias: job name (default: bml)
    :return: job id
    """
    rst1_name = result_prefix + "-" + branch_num + "-n-Result"
    rst2_name = result_prefix + "-" + branch_num + "-a-Result"
    clt1 = file_utilities.check_filename(clt1)
    clt2 = file_utilities.check_filename(clt2)
    base = os.path.split(clt1)[0]
    cmdfile = "codeml.txt"; tnum = 1; gmem = 2
    cmds = "cd %s\n" \
           "/sternadi/home/volume1/taliakustin/software/paml4.8/bin/codeml %s\n" \
           "mv rst %s\n" \
           "/sternadi/home/volume1/taliakustin/software/paml4.8/bin/codeml %s\n" \
           "mv rst %s\n" % (base, clt1, rst1_name, clt2, rst2_name)
    pbs_jobs.create_pbs_cmd(cmdfile, alias, tnum, gmem, cmds)
    job_id = pbs_jobs.submit(cmdfile)
    return job_id

if __name__ == "__main__":
    main()