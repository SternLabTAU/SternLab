#! /usr/local/python_anaconda/bin/python3.4

from optparse import OptionParser
from file_utilities import check_dirname, check_filename, change_filename
import glob
from seqFileAnalyzer import get_consensus_percentage
import pandas as pd
from PAML_utilities import write_ctl_file


def change_phy_aln_to_aln():
    # used it already
    files = glob.glob("/sternadi/nobackup/volume1/phyVirus/*.phy_aln*")
    print("going to change %i files" % len(files))
    for f in files:
        out = "_aln".join(f.split(".phy_aln"))
        change_filename(f, out)

    print("Done")

def remove_viridae():
    files = glob.glob("/sternadi/nobackup/volume1/phyVirus/*viridae*")
    print("going to change %i files" % len(files))
    for f in files:
        out = "".join(f.split("viridae"))
        change_filename(f, out)

    print("Done")

def remove_0_():
    files = glob.glob("/sternadi/nobackup/volume1/phyVirus/*/*_0_*")
    print("going to change %i files" % len(files))
    for f in files:
        out = ".".join(f.split("_0_", 1))
        change_filename(f, out)

    print("Done")

def change_aln_best_to_aln():
    files = glob.glob("/sternadi/nobackup/volume1/phyVirus/*_aln.best*")
    print("going to change %i files" % len(files))
    for f in files:
        out = ".aln".join(f.split("_aln.best", 1))
        change_filename(f, out)

    print("Done")

def change_aln_best_to_aln_point():
    files = glob.glob("/sternadi/nobackup/volume1/phyVirus/*.aln.best*")
    print("going to change %i files" % len(files))
    for f in files:
        out = ".aln".join(f.split(".aln.best", 1))
        change_filename(f, out)

    print("Done")


if __name__ == "__main__":
    main()


