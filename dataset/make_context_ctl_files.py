#! /usr/local/python_anaconda/bin/python3.4

from optparse import OptionParser
from file_utilities import check_dirname, check_filename
import glob
from seqFileAnalyzer import get_consensus_percentage
import pandas as pd
from PAML_utilities import write_ctl_file


def main():
    parser = OptionParser("usage: %prog [options]")
    parser.add_option("-d", "--dir", dest="dataset_dir", help="dir of dataset files")
    parser.add_option("-c", "--c", dest="context_dir", help="dir of context alignent files")
    (options, args) = parser.parse_args()

    dataset_dir = options.dataset_dir
    dataset_dir = check_dirname(dataset_dir)

    context_dir = options.context_dir
    context_dir = check_dirname(context_dir)


    aln_files = glob.glob(dataset_dir + "/*/*best.fas")
    for aln_file in aln_files:
        print(aln_file)

        basename_dataset_dir = aln_file.split(".phy")[0].split("_aln")[0].split(".aln")[0]
        tree = glob.glob(basename_dataset_dir + "[_.A-z]*tree*")[0]
        rooted_tree = glob.glob(basename_dataset_dir + ".r.tree")[0]
        basename = aln_file.split(".phy")[0].split("_aln")[0].split(".aln")[0].split(dataset_dir)[-1]
        context_alns = glob.glob(context_dir + "/%s[_.A-z]*.aln" %basename)


        for context_aln in context_alns:
            context_basename = context_aln.split(".aln")[0]
            ctl_gtr = context_basename + "_gtr.ctl"
            ctl_unr = context_basename + "_unr.ctl"
            output_gtr = context_basename + "_gtr.mlb"
            output_unr = context_basename + "_unr.mlb"

            write_ctl_file(ctl_gtr, context_aln, tree, output_gtr, 7, fix_alpha="1", alpha="0")
            write_ctl_file(ctl_unr, context_aln, rooted_tree, output_unr, 8, fix_alpha="1", alpha="0")


        


if __name__ == "__main__":
    main()


