#! /usr/local/python_anaconda/bin/python3.4

from optparse import OptionParser
from file_utilities import check_dirname, check_filename
import glob
import re
from seqFileAnalyzer import get_consensus_percentage


def main():
    parser = OptionParser("usage: %prog [options]")
    parser.add_option("-d", "--dir", dest="dir", help="dir input fasta files")
    parser.add_option("-s", "--subtypes", dest="subtypes", help="subtypes to keep, seperated by comma. example A,B")
    (options, args) = parser.parse_args()

    dir = options.dir
    dir = check_dirname(dir)

    subtypes = options.subtypes.split(",")
    input_files = glob.glob(dir + "/*.fasta")


    for file in input_files:
        fasta = open(file, "r").read()
        basename = file.split(".fasta")[0]
        for subtype in subtypes:
            pattern = re.compile(">%s[^>]*" % subtype)
            results = pattern.findall(fasta)
            print("%s: %s: %s" % (file, subtype, str(len(results))))
            fasta_out = "".join(results)
            output_file = basename + "_%s.aln" % subtype
            output = open(output_file, "w")
            output.write(fasta_out)
            output.close()
    print("splited HIV files to subtypes")

if __name__ == "__main__":
    main()


