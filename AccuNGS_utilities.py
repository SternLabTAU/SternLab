#! /usr/local/python/anaconda_python-3.6.1

"""
@Author: odedkushnir

"""

import os
import pbs_runners
import glob
import pandas as pd
from cirseq_utilities import *



def main():
    # parser = OptionParser("usage: %prog [options]")
    # parser.add_option("-i", "--input_file", dest="input_file", help="fastq file")
    # parser.add_option("-o", "--output_dir", dest="output_dir", help="output dir")
    # (options, args) = parser.parse_args()
    # file = options.input_file
    # output_dir = options.output_dir

    # 1st thing to do is to index the output files from the Nextseq
    # csv_file = "/Volumes/STERNADILABHOME$/volume3/okushnir/RNAseq/180725_M04473_0026_000000000-BT9GF/index.csv"
    # fastq_path = "/Volumes/STERNADILABHOME$/volume3/okushnir/RNAseq/180725_M04473_0026_000000000-BT9GF/"
    # output_dir = "/Volumes/STERNADILABHOME$/volume3/okushnir/RNAseq/180725_M04473_0026_000000000-BT9GF/indexed"
    # index(csv_file, fastq_path, output_dir)

    # 2nd is to clean the files from --
    # trim(file, output_dir)

    # input_dir = ("/sternadi/datasets/volume1/180503_OST_FINAL_03052018/indexed/")
    # output_dir = ("/sternadi/home/volume3/okushnir/AccuNGS/180503_OST_FINAL_03052018/sheri_clean/")
    # files = glob.glob(input_dir + "*.fastq")
    # for f in files:
    #     trim(f, output_dir)

    # # 3rd merge the files
    # merge_paired_files("/sternadi/home/volume3/okushnir/AccuNGS/180503_OST_FINAL_03052018/clean/",
    #                    "/sternadi/home/volume3/okushnir/AccuNGS/180503_OST_FINAL_03052018/merged")
    #

    # # 4th run pipeline:
    # folders = glob.glob("/sternadi/home/volume3/okushnir/AccuNGS/180503_OST_FINAL_03052018/merged/*")
    #
    # for d in folders:
    #     output_dir = d + "/q30_3UTR"
    #     cmd = "python /sternadi/home/volume1/shared/SternLab/pipeline_runner.py -i %s -o %s -r /sternadi/home/volume3/okushnir/" \
    #           "AccuNGS/180503_OST_FINAL_03052018/merged/HRVB14_from_pWR3.26_1-7212.fasta -NGS_or_Cirseq 2 -rep 2  -q 30" % (d, output_dir)
    #     pbs_runners.script_runner(cmd, alias="pipeline_d")

    #5th analyze the freqs
        # add mutation types
    sample = "P2.SRR1036477.V3.2"
    suffix = "%s.freqs" % sample
    freqs_file = "/volumes/STERNADILABHOME$/volume3/okushnir/Cirseq/PV/Mahoney/P2/P2.SRR1036477.V3.2.freqs"
    # freqs_file = "/volumes/STERNADILABHOME$/volume3/okushnir/AccuNGS/RV-P7_L001-ds.32944248b7874527aa7daeed6203d1da/merged/%s/q30/%s" % \
    #              (sample, suffix)
    virus = "PV"
    seq_meth = "CirSeq"

    if virus == "CVB3":
        ncbi_id ="M16572"
    if virus == "RVB14":
        ncbi_id = "NC_001490"
    if virus == "PV":
        ncbi_id ="V01149"
    if not os.path.isfile(freqs_file[0:-5] + "with.mutation.type.freqs"):
         append_mutation = find_mutation_type(freqs_file, ncbi_id)
    freqs_file_mutations = freqs_file[0:-5] + "with.mutation.type.freqs"

    #freqs_mutation = find_mutation_type("/Volumes/STERNADILABHOME$/volume3/okushnir/AccuNGS/RV-P7_L001-ds.32944248b7874527aa7daeed6203d1da/merged/RV-p71/q30/RV-p71.freqs", "NC_001490")
    #data_mutation.csv file just for RVBp7
    #transition_mutation(freqs_file_mutations, "/Volumes/STERNADILABHOME$/volume3/okushnir/AccuNGS/RV-P7_L001-ds.32944248b7874527aa7daeed6203d1da/merged/RV-p71/q30/")

    # data_mutation.csv file for p1, p7 ant control
    sample_file1 = "/Volumes/STERNADILABHOME$/volume3/okushnir/AccuNGS/RV-P7_L001-ds.32944248b7874527aa7daeed6203d1da/merged/RV-p71/q30/RV-p71.with.mutation.type.freqs"
    sample_file2 = "/Volumes/STERNADILABHOME$/volume3/okushnir/AccuNGS/180503_OST_FINAL_03052018/merged/RV-p11/q30_3UTR_new/RV-p11.with.mutation.type.freqs"
    sample_file3 = "/Volumes/STERNADILABHOME$/volume3/okushnir/AccuNGS/180503_OST_FINAL_03052018/merged/RV-p12/q30_3UTR_new/RV-p12.with.mutation.type.freqs"
    sample_file4 = "/Volumes/STERNADILABHOME$/volume3/okushnir/Cirseq/PV/Mahoney/P3/20170907_q23r2_blastn/PV-p3.1036617.with.mutation.type.freqs"
    control_file = "/Volumes/STERNADILABHOME$/volume3/okushnir/AccuNGS/180503_OST_FINAL_03052018/merged/RV-IVT/q30_3UTR_new/RV-IVT.with.mutation.type.freqs"

    label_control = "RNA Control"
    pass_sample0 = 0
    rep_sample0 = 1
    label_sample1 = "p7 Replica #1"
    pass_sample1 = 7
    rep_sample1 = 1
    label_sample2 = "p1 Replica #1"
    pass_sample2 = 1
    rep_sample2 = 1
    label_sample3 = "p1 Replica #2"
    pass_sample3 = 1
    rep_sample3 = 2
    label_sample4 = "PV Mahoney p3"
    pass_sample4 = 3
    rep_sample4 = 1


    print ("loading " + sample_file1 + " as sample")
    data_mutations1 = pd.read_table(sample_file1)
    data_mutations1["source"] = label_sample1

    print ("loading " + sample_file2 + " as sample")
    data_mutations2 = pd.read_table(sample_file2)
    data_mutations2["source"] = label_sample2

    print("loading " + sample_file3 + " as sample")
    data_mutations3 = pd.read_table(sample_file3)
    data_mutations3["source"] = label_sample3

    print("loading " + sample_file4 + " as sample")
    data_mutations4= pd.read_table(sample_file4)
    data_mutations4["source"] = label_sample4

    print("loading " + control_file + " as homogeneous control")
    data_control = pd.read_table(control_file)
    data_control["source"] = label_control

    data = pd.concat([data_control, data_mutations1, data_mutations2, data_mutations3, data_mutations4])

    transition_mutation(data,
                        "/Volumes/STERNADILABHOME$/volume3/okushnir/Cirseq/PV/Mahoney/P3/20170907_q23r2_blastn/")

    #6th run variant_caller localy to check context mutations


def index(csv_file, fastq_path, output_dir):
    indexing = pd.read_csv(csv_file)  # your_csv_dir
    sample_id = list(indexing.SampleID)
    index_1 = list(indexing.Index1Sequence)
    index_2 = list(indexing.Index2Sequence)

    index_dic = {}
    for i in range(len(sample_id)):
        index_dic[index_1[i] + "+" + index_2[i]] = sample_id[i]

    files = glob.glob(fastq_path)

    for f in files:
        for key in index_dic.keys():
            output = index_dic[key] + f.split("/")[-1].split(".fastq")[0].split("S")[1] + ".fastq"
            output_file = output_dir + output
            pbs_runners.script_runner("grep '%s' %s -A 3 > %s" % (key, f, output_file), alias="OST_Sample")


def trim(input_file, output_dir):
    new_f = output_dir + os.path.basename(input_file)
    new_text = ""
    with open(input_file) as infile:
        for l in infile:
            if l != "--\n":
                new_text += l
    out = open(new_f, "w")
    out.write(new_text)
    out.close()


def merge_paired_files(input_dir, out_dir):
    """
    :param input_dir: tmp directory path of the cirseq pipeline analysis
    :param out_dir: the output directory path
    :return: merged files
    """


    files1 = glob.glob(input_dir + "*_*_R1_001.fastq")
    lst_files = []
    for fastq1 in files1:
        filename = os.path.basename(fastq1)
        lane = filename.split("_")[1]
        sample = filename.split("_")[0]
        fastq2 = ("%s%s_%s_R2_001.fastq") % (input_dir, sample, lane)
        if not os.path.exists(out_dir):
            out_dir = os.system("mkdir %s" % out_dir)
        output_file = "%s/%s_%s_merged.fastq" % (out_dir, sample, lane)

        pbs_runners.script_runner(
            "python /sternadi/home/volume3/okushnir/SternLab/scripts/merge_fastq_files.py -f %s -e %s -o %s -r 60" % (
            fastq1, fastq2, output_file), alias="merge_RV")


if __name__ == "__main__":
    main()

