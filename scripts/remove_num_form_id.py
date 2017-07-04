"""
@Author: odedkushnir

"""

import sys
from optparse import OptionParser
from file_utilities import check_dirname


def main():
    parser = OptionParser("usage: %prog [options]")
    parser.add_option("-f", "--file1", dest="file1", help="fastq1 file")
    parser.add_option("-o", "--output_file", dest="output_file", help="output fastq merged file")
    (options, args) = parser.parse_args()
    file1 = options.file1
    output_file = options.output_file

    rename_id_fastq(file1, output_file)


def rename_id_fastq(file1, output_file):

    with open(file1, 'r') as fastq1, open(output_file, 'w') as out:
        f1 = fastq1.readlines()
        for i in range(0, len(f1)-1, 4):
            f1[i].strip()
            if f1[i].split(" ")[0][0] == '@':
                id = f1[i].split(" ")[0]
                seq = f1[i + 1]
                dir = f1[i + 2]
                qual = f1[i + 3]
                out.write(id+'\n' + seq + dir + qual)
            else:
                print("Problem" + str(i))

if __name__ == "__main__":
    main()
