#! /usr/local/python_anaconda/bin/python3.4


from file_utilities import check_filename, check_dirname
import argparse
import os
import datetime
import subprocess
import matplotlib.pyplot as plt
plt.switch_backend('agg')
import seaborn as sns
import pandas as pd

sns.set_context('talk')


def main(args):


    pipeline_path = "/sternadi/home/volume1/shared/SternLab/pipeline/runner.pl"

    NGS_or_Cirseq = args.NGS_or_Cirseq


    print("Pipeline to run: %s" % pipeline_path)

    input_dir = args.input_dir
    input_dir = check_dirname(input_dir)
    output = args.output
    output = check_dirname(output, Truedir=False)
    reference = args.ref
    reference = check_filename(reference)


    start = args.start
    end = args.end
    if start not in [1, 2, 3]:
        raise Exception("Not a valid start step - needs to be between 1:3")
    if end not in [2, 3, 4]:
        raise Exception("Not a valid end step - needs to be between 2:4")

    type_of_input_file = args.type_of_input_file
    gaps = args.gaps
    if gaps not in ["Y", "N"]:
        raise Exception("Not a valid gap - must be Y or N")
    q_score = args.q_score

    if q_score == None:
        if NGS_or_Cirseq == 1:
            q_score = 30
        else:
            q_score =23

    blast_id = args.blast

    path_to_save_pipeline_summary = output + "/pipeline_summary.txt"
    print(start, end, q_score, blast_id, NGS_or_Cirseq)

    cmd = "perl {} {} {} {} {} {} {} {} {} {} {}".format(pipeline_path, input_dir, output, reference,
                                                            start, end, type_of_input_file, gaps, NGS_or_Cirseq,
                                                            q_score, blast_id)

    print("running this pipeline command:")
    print(cmd)
    os.system(cmd)


    # get additional statistics about this running
    os.chdir(os.path.join(output, "tmp"))
    os.system("pwd")

    # number of reads that were mapped only once
    only_once_reads = subprocess.getoutput("grep -P '^1\t' *stats  -h | awk '{sum+=$2}END{print sum}'")
    #number of reads that are "contributing to frequency counts”
    freq_contr = subprocess.getoutput("grep 'reads contributing to frequency counts' -h *stats | awk '{sum+=$1}END{print sum}'")
    #number of bases called
    num_based_called = subprocess.getoutput(
        "grep 'num bases called' *stats | awk -F = '{sum+=$2}END{print sum}'")
    #number of reads that were mapped to reference
    num_reads_mapped = subprocess.getoutput(
        "cat *blast | awk '{print $1}' | sort | uniq | wc -l")



    with open(path_to_save_pipeline_summary, "w") as o:
        o.write("---- Pipeline running -----\n")
        o.write("{}\n\n".format(datetime.datetime.now()))
        o.write("Pipeline command used: {}\n\n".format(cmd))
        o.write("Number of reads that were mapped only once: {}\n".format(int(only_once_reads)))
        o.write("Number of reads that are contributing to frequency count: {}\n".format(int(freq_contr)))
        o.write("Number of bases called: {}\n".format(int(num_based_called)))
        o.write("Number of reads mapped to reference: {}\n".format(int(num_reads_mapped)))

    #get back to the freq file directory

    os.chdir(output)
    # create a simple coverage plot
    freq_file_path = os.path.join(output, [f for f in os.listdir(output) if ".freq" in f][0])
    freq_file_path = check_filename(freq_file_path)
    label = os.path.basename(freq_file_path).split('.')[0]

    df = pd.read_csv(freq_file_path, sep='\t')

    df = df.drop_duplicates("Pos")
    plt.plot(df['Pos'].values, df['Read_count'].values, label=label, color='darkorange')
    plt.title("Coverage {}".format(label), fontsize=16)
    plt.xlabel("Position in the genome (bp)")
    plt.ylabel("Read count")
    plt.savefig(os.path.join(output, 'coverage.png'), format='png')



    print("Ran pipeline")


if __name__ == "__main__":
    parser = argparse.ArgumentParser()
    parser.add_argument("-i", "--input_dir", type=str, help="input directory with fastq.gz files",
                        required=True)
    parser.add_argument("-o", "--output", type=str, help="a path to an output directory", required=True)
    parser.add_argument("-r", "--ref", type=str, help="a path to a genome reference seq file (fasta)",
                        required=True)
    parser.add_argument("-s", "--start", type=int, help="start step number. default=1", default=1)
    parser.add_argument("-e", "--end", type=int, help="end step number. default=4", default=4)

    parser.add_argument("-t", "--type_of_input_file", type=str, help="the type of the input file, default=fastq",
                        default='f', required=False)
    parser.add_argument("-g", "--gaps", type=str, help="refer to gaps? Y/N default Y", default='Y',
                        required=False)
    parser.add_argument("-NGS_or_Cirseq", "--NGS_or_Cirseq", type=int, help="NGS/Cirseq? type 2 for Cirseq or 1 for NGS",
                        required=True)
    parser.add_argument("-q", "--q_score", type=int, help="Q-score cutoff, default =23 for cirseq and 30 for NGS")
    parser.add_argument("-b", "--blast", type=int, help="% blast id, default=85",
                        default=85)
    args = parser.parse_args()
    if not vars(args):
        parser.print_help()
        parser.exit(1)

    main(args)
