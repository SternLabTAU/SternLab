#! /usr/local/python_anaconda/bin/python3.4

from optparse import OptionParser
from file_utilities import check_dirname, check_filename
import glob
from seqFileAnalyzer import get_consensus_percentage
import pandas as pd


def main():
    parser = OptionParser("usage: %prog [options]")
    parser.add_option("-d", "--dir", dest="dataset_dir", help="dir of dataset files")
    parser.add_option("-o", "--output", dest="output_file", help="output file name")
    (options, args) = parser.parse_args()

    dataset_dir = options.dataset_dir
    dataset_dir = check_dirname(dataset_dir)

    output_file = options.output_file
    output_file = check_filename(output_file, Truefile=False)

    aln_files = glob.glob(dataset_dir + "/*/*best.fas")


    df = pd.DataFrame(columns=["filename", "1", "0.9", "0.8", "0.7", "0.6", "0.5", "0.4", "0.3", "0.2"])

    for aln_file in aln_files:
        consensus, consensus_percentage = get_consensus_percentage(aln_file)
        filename = aln_file.split(dataset_dir)[-1]
        print(aln_file)
        df = df.append({"filename":filename, "1": consensus_percentage[1],
                       "0.9": consensus_percentage[0.9], "0.8":consensus_percentage[0.8],
                       "0.7": consensus_percentage[0.7], "0.6":consensus_percentage[0.6],
                       "0.5": consensus_percentage[0.5], "0.4": consensus_percentage[0.4],
                       "0.3": consensus_percentage[0.3], "0.2": consensus_percentage[0.2]},
                       ignore_index=True)


    df.to_csv(output_file, index=False)

if __name__ == "__main__":
    main()


