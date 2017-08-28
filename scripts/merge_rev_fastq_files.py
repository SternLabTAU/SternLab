import sys
from optparse import OptionParser
from Bio.Seq import Seq
from file_utilities import check_filename, check_dirname
import os



def main():
    parser = OptionParser("usage: %prog [options]")
    parser.add_option("-f", "--file1", dest="file1", help="fastq1 file")
    parser.add_option("-e", "--file2", dest="file2", help="fastq2 file")
    parser.add_option("-o", "--output_file", dest="output_file", help="output fastq merged file")
    (options, args) = parser.parse_args()
    file1 = options.file1
    file2 = options.file2
    output_file = options.output_file

    concat_fastq(file1, file2, output_file)


def concat_fastq(file1, file2, output_file):
    """
    merge 2 files of fatsq to 1 file base on the id
    :param file1: fastq file 1 you want to merge
    :param file2: fastq file 2 you want to merge
    :param output_file: merged file
    :return:
    """
    with open(file1, 'r') as fastq1, open(file2, 'r') as fastq2, open(output_file, 'w') as output_file:
        f1 = fastq1.readlines()
        f2 = fastq2.readlines()
        for i in range(0, len(f1), 4):
            f1[i].strip()
            f2[i].strip()
            if f1[i].split(" ")[0] == f2[i].split(" ")[0]:
                merge_id = f1[i].split(" ")[0]
                merge_seq = str(Seq(f1[i + 1].strip()) + 'N' + Seq(f2[i + 1].strip()).reverse_complement())
                merge_dir = f1[i + 2].strip()
                merge_qual = f1[i + 3].strip() + '!' + f2[i + 3].strip()
                output_file.write(merge_id + '\n' + merge_seq + '\n' + merge_dir + '\n' + merge_qual + '\n')
            else:
                print("Problem:" + str(i))

if __name__ == "__main__":
    main()