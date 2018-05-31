

"""
@Author: odedkushnir

"""

import os
import pbs_runners
import glob
import pandas as pd



def main():
    # parser = OptionParser("usage: %prog [options]")
    # parser.add_option("-i", "--input_file", dest="input_file", help="fastq file")
    # parser.add_option("-o", "--output_dir", dest="output_dir", help="output dir")
    # (options, args) = parser.parse_args()
    # file = options.input_file
    # output_dir = options.output_dir

    # 1st thing to do is to index the output files from the Nextseq
    #index(csv_file, fastq_path, output_dir)

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
    folders = glob.glob("/sternadi/home/volume3/okushnir/AccuNGS/180503_OST_FINAL_03052018/merged/*")

    for d in folders:
        output_dir = d + "/q30_3UTR"
        cmd = "python /sternadi/home/volume1/shared/SternLab/pipeline_runner.py -i %s -o %s -r /sternadi/home/volume3/okushnir/" \
              "AccuNGS/180503_OST_FINAL_03052018/merged/HRVB14_from_pWR3.26_1-7212.fasta -NGS_or_Cirseq 2 -rep 2  -q 30" % (d, output_dir)
        pbs_runners.script_runner(cmd, alias="pipeline_d")


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
            "python /sternadi/home/volume1/sheri/SternLab/scripts/merge_fastq_files.py -f %s -e %s -o %s -r 60" % (
            fastq1, fastq2, output_file), alias="merge_RV")


if __name__ == "__main__":
    main()

