#! /usr/local/python_anaconda/bin/python3.4

from optparse import OptionParser
from seqFileTools import remove_description
from file_utilities import check_filename

def main():
    parser = OptionParser("usage: %prog [options]")
    parser.add_option("-i", "--input", dest="input", help="file that contains the initial alignment/dataset")
    parser.add_option("-o", "--output", dest="output", help="output file - optional", default=None)
    (options, args) = parser.parse_args()

    input = options.input
    input = check_filename(input)
    output = options.output
    if output != None:
        output = check_filename(output)

    output = remove_description(input, output)
    print("saved %s without description in %s" %(input, output))



if __name__ == "__main__":
    main()

