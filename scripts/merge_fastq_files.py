#! /usr/local/python_anaconda/bin/python3.4


import sys
sys.path.insert(0, "/sternadi/home/volume1/taliakustin/SternLab/")
from optparse import OptionParser
from file_utilities import check_dirname



def main():
    parser = OptionParser("usage: %prog [options]")
    parser.add_option("-f", "--file1", dest="file1", help="fastq1 file")
    parser.add_option("-e", "--file2", dest="file2", help="fastq2 file")
    parser.add_option("-o", "--output_file", dest="output_file", help="output fastq merged file")
    parser.add_option("-n", "--number_of_ns", dest="number_of_Ns",  default = 1, help="number of N's to use as seperators")
    parser.add_option("-s", "--sequencer", dest="sequencer",  default = "M", help="Did the fastq came from Miseq (M) or Nextseq (N)")

    (options, args) = parser.parse_args()
    file1 = options.file1
    file2 = options.file2
    output_file = options.output_file
    number_of_Ns = options.number_of_Ns
    sequencer = options.sequencer

    concat_fastq(file1, file2, output_file, sequencer, number_of_Ns)


def concat_fastq(file1, file2, output_file, sequencer, number_of_Ns = 1):
    """
    merge 2 files of fatsq to 1 file base on the id
    :param file1: fastq file 1 you want to merge
    :param file2: fastq file 2 you want to merge
    :param output_file: merged file
    :param number_of_Ns: how many Ns to put between the reads
    :param sequencer
    :return:
    """
    number_of_Ns = int(number_of_Ns)
    separator = "N" * number_of_Ns
    quality_separator = "!" * number_of_Ns
    if sequencer == "M":
        prefix = "@M"
    elif sequencer == "N":
        prefix =="@N"
    with open(file1, 'r') as fastq1, open(file2, 'r') as fastq2, open(output_file, 'w') as output_file:
        f1 = fastq1.readlines()
        f2 = fastq2.readlines()
        for i in range(0, len(f1)):
            if prefix not in f1[i]:
                continue
            if prefix not in f2[i]:
                continue
            f1[i] = f1[i].strip()
            f2[i] = f2[i].strip()
            if f1[i].split(" ")[0] == f2[i].split(" ")[0]:
                merge_id = f1[i].split(" ")[0]
                merge_seq = f1[i + 1].strip() + separator + f2[i + 1].strip()
                merge_dir = f1[i + 2].strip()
                merge_qual = f1[i + 3].strip() + quality_separator + f2[i + 3].strip()
                output_file.write(merge_id + '\n' + merge_seq + '\n' + merge_dir + '\n' + merge_qual + '\n')
            else:
                print("Problem:" + str(i))

if __name__ == "__main__":
    main()

