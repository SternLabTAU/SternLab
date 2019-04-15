#! /usr/local/python_anaconda/bin/python3.4
import sys, gzip
from Bio import SeqIO,Seq
from Bio.Seq import MutableSeq
from Bio.SeqRecord import SeqRecord
from optparse import OptionParser
import pandas as pd
import numpy as np

def replace_base(row, mut_reference):
    mut_reference[int(row["Pos"])-1]=row["Base"]

def main():
    parser = OptionParser("usage: %prog [options]\nTry running %prog --help for more information")
    parser.add_option("-f", "--fasta", dest="fasta", help="reference fasta file")
    parser.add_option("-p", "--freqs", dest="freqs", help="frequency file")
    parser.add_option("-o", "--output_file", dest="output_file", help="output fasta file")
    parser.add_option("-c", "--min_coverage", dest="coverage", help="minimal coverage of position to be updated (default: 1000)", type="int")
    (options, args) = parser.parse_args()
    fasta = options.fasta
    freqs = options.freqs
    output_file = options.output_file
    min_coverage = 1000
    if options.coverage:
        min_coverage= options.coverage

    if options.fasta is None or options.fasta is None or options.output_file is None:
        parser.error("Missing file input")

    reference=None
    input_record=None

    mutable_reference=None

    with open(fasta, "rU") as handle:
        for record in SeqIO.parse(handle, "fasta"):
            reference=record.seq
            input_record=record

    mutable_reference=MutableSeq(str(reference))

    data_freqs = pd.read_table(freqs)
    data_freqs = data_freqs[data_freqs["Read_count"] >= min_coverage]
    data_freqs = data_freqs[data_freqs["Rank"] == 0]
    #data_freqs = data_freqs[data_freqs["Base"] != data_freqs["Ref"]]
    data_freqs = data_freqs[data_freqs["Pos"] == np.round(data_freqs['Pos'])] #remove insertion
    data_freqs = data_freqs[data_freqs["Base"] != "-"] #remove deletion
    data_freqs.apply(replace_base, args=(mutable_reference,), axis=1) #updates mutable reference to hold correct consensus

    new_sequence=Seq.Seq(str(mutable_reference))
    input_record.seq=new_sequence
    with open(output_file, "w") as output_handle:
        SeqIO.write(input_record, output_handle, "fasta")

if __name__ == "__main__":
    main()