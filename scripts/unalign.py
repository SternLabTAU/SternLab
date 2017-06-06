#! /usr/local/python_anaconda/bin/python3.4

from optparse import OptionParser
from file_utilities import check_filename
from seqFileTools import unalign

def main():
    parser = OptionParser("usage: %prog [options]")
    parser.add_option("-i", "--input", dest="input", help="file that contains the initial alignment/dataset")
    parser.add_option("-o", "--output", dest="output_file", help="output file for unaligned file", default=None)
    parser.add_option("-f", "--format", dest="format", help="format of input file", default="fasta")
    parser.add_option("-g", "--gap", dest="gap", help="gap charachter to unalign", default="-")
    (options, args) = parser.parse_args()

    input = options.input
    output = options.output
    format = options.format
    gap = options.gap


    #check filnames
    if input == None:
        parser.error("you must specify a filname")
    input = check_filename(input)
    if output != None:
        output = check_filename(output)

    output = unalign(input, in_format=format, gap=gap, outfile=output)
    print("unaligned %s to %s" % (input, output))


if __name__ == "__main__":
    main()

