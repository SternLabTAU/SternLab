import PAML_utilities
import os
import seqFileTools
import pbs_runners
import pbs_jobs
import file_utilities
from os import path
from optparse import OptionParser


def main():
    parser = OptionParser("usage: %prog[options]")
    parser.add_option("-i", "--gene_folder", dest="GENES_FOLDER", help="The folder that contains all candida genes subfolders")

    (options, args) = parser.parse_args()

    GENES_FOLDER = options.GENES_FOLDER

    create_files_for_genes(GENES_FOLDER)


def create_files_for_genes(folder):
    for filename in os.listdir(folder):
        if filename.endswith(".fasta"):
            file = open("/sternadi/home/volume2/ella/Candida/Genes/" + filename, "r")
            #skip the reference gene
            file.readline()
            file.readline()
            folder = "/sternadi/home/volume2/ella/Candida/Genes/NoRef/" + filename[:-6]
            os.makedirs(folder)
            output = open(folder + "/" + filename, "w")

            # write all sequences with A in the begining of the name, and without the stop codon at the end
            for line in file:
                if ">" in line:
                    output.write(line.replace(">",">A"))
                else:
                    output.write(line[:-4].upper() + "\n")
            output.close()
            file.close()

            #create phylip
            seqFileTools.convert_fasta_to_phylip(
                folder + "/" + filename, outfile=None)

            #create CLT
            M8a_CTL = folder + "/" + filename.replace(".fasta", "-M8a.CLT")
            M8_CTL = folder + "/" + filename.replace(".fasta", "-M8.CLT")
            PAML_utilities.candida_write_ctl_codeml_file(
                M8a_CTL,
                folder + "/" + filename.replace(".fasta", ".phy"),
                "/sternadi/home/volume2/ella/Candida/Trees/SNP_seq_without_ref.phy_phyml_tree.txt",
                folder + "/" + filename.replace(".fasta", "-M8a-Results.txt"),
                1)
            PAML_utilities.candida_write_ctl_codeml_file(
                M8_CTL,
                folder + "/" + filename.replace(".fasta", ".phy"),
                "/sternadi/home/volume2/ella/Candida/Trees/SNP_seq_without_ref.phy_phyml_tree.txt",
                folder + "/" + filename.replace(".fasta", "-M8-Results.txt"),
                0)

            #run codeml
            candida_codeml_runner(M8a_CTL, M8_CTL, filename[:-6])


def candida_codeml_runner(ctl1, ctl2, gene_name, alias = "cml"):
    """
    run baseml program from PAML on cluster
    :param ctl: ctl file path
    :param alias: job name (default: bml)
    :return: job id
    """
    rst1_name = gene_name + "-M8a-Result"
    rst2_name = gene_name + "-M8-Result"
    ctl1 = file_utilities.check_filename(ctl1)
    ctl2 = file_utilities.check_filename(ctl2)
    base = os.path.split(ctl1)[0]
    cmdfile = "codeml.txt"; tnum = 1; gmem = 2
    cmds = "cd %s\n" \
           "/sternadi/home/volume1/taliakustin/software/paml4.8/bin/codeml %s\n" \
           "mv rst %s\n" \
           "/sternadi/home/volume1/taliakustin/software/paml4.8/bin/codeml %s\n" \
           "mv rst %s\n" % (base, ctl1, rst1_name, ctl2, rst2_name)
    pbs_jobs.create_pbs_cmd(cmdfile, alias, tnum, gmem, cmds)
    job_id = pbs_jobs.submit(cmdfile)
    return job_id





def run_shorts(folder):
    fixed = 0
    problematic = ["Ca21chr1-orf19.4923.1","Ca21chr3-orf19.6828.1","Ca21chr5-orf19.961.2","Ca21chrR-orf19.1738.1","Ca21chrR-orf19.4380.1"]
    for filename in os.listdir(folder):
        if filename in problematic:
            new = "/sternadi/home/volume2/ella/"+ filename + "/" + filename + ".fasta"
            print(new)
            os.makedirs("/sternadi/home/volume2/ella/"+ filename)
            os.rename("/sternadi/home/volume2/ella/Candida/Genes/NoRef/"+filename + "/" + filename + ".fasta",
                     new)

            seqFileTools.convert_fasta_to_phylip(new, outfile=None)

            M8a_CTL = new.replace(".fasta", "-M8a.CLT")
            M8_CTL = new.replace(".fasta", "-M8.CLT")
            PAML_utilities.candida_write_ctl_codeml_file(
                M8a_CTL,
                new.replace(".fasta", ".phy"),
                "/sternadi/home/volume2/ella/Candida/Trees/SNP_seq_without_ref.phy_phyml_tree.txt",
                new.replace(".fasta", "-M8a-Results.txt"),
                1)
            PAML_utilities.candida_write_ctl_codeml_file(
                M8_CTL,
                new.replace(".fasta", ".phy"),
                "/sternadi/home/volume2/ella/Candida/Trees/SNP_seq_without_ref.phy_phyml_tree.txt",
                new.replace(".fasta", "-M8-Results.txt"),
                0)

            candida_codeml_runner(M8a_CTL, M8_CTL, filename)

            fixed += 1
    print("created shorts for: " + str(fixed))


if __name__ == "__main__":
    main()
