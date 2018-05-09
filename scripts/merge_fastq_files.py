#! /usr/local/python_anaconda/bin/python3.4
from Bio import SeqIO,Seq
from Bio.SeqRecord import SeqRecord
import sys, gzip
from optparse import OptionParser

def merger_generator(forward_handle,reverse_handle, rep_length):
    for a, b in zip (SeqIO.parse(forward_handle, "fastq"),SeqIO.parse(reverse_handle, "fastq")):
        if (a.id.split(" ")[0] !=  b.id.split(" ")[0]):
                    print("Problem:" + str(i))
        new_seq_id=a.id.split(" ")[0]
        new_seq_str = str(a.seq) + ("N"*rep_length) + str(b.seq)
        a_quals=a.letter_annotations["phred_quality"]
        b_quals=b.letter_annotations["phred_quality"]
        new_seq_qual= a_quals+[1.0 for a in range(rep_length)]+b_quals

        new_seq=SeqRecord(new_seq_str,id=new_seq_id,description="",letter_annotations={"phred_quality":new_seq_qual})
        yield new_seq

def main():
    parser = OptionParser("usage: %prog [options]\nTry running %prog --help for more information")
    parser.add_option("-f", "--file1", dest="file1", help="fastq1 file")
    parser.add_option("-e", "--file2", dest="file2", help="fastq2 file")
    parser.add_option("-o", "--output_file", dest="output_file", help="output fastq merged file")
    parser.add_option("-r", "--rep_length", dest="rep_length", help="amount of N bases to repeat (default: 1)", type="int")
    (options, args) = parser.parse_args()
    file1 = options.file1
    file2 = options.file2
    output_file = options.output_file
    rep_length = 1
    if options.rep_length:
        rep_length = options.rep_length

    if options.file1 is None or options.file2 is None or options.output_file is None:
        parser.error("Missing file input")
    if file1.endswith(".gz"):
        concat_fastq_gz(file1, file2, output_file, rep_length)
    else:
        concat_fastq(file1, file2, output_file, rep_length)

def concat_fastq(file1, file2, output_file, rep_length):
    with open(file1, "r") as forward_handle:
        with open(file2, "r") as reverse_handle :
            with open (output_file,"w") as merged_handle:
                SeqIO.write(merger_generator(forward_handle,reverse_handle, rep_length), merged_handle, "fastq")

def concat_fastq_gz(file1, file2, output_file, rep_length):
    with gzip.open(file1, "rt") as forward_handle:
        with gzip.open(file2, "rt") as reverse_handle :
            with gzip.open (output_file,"wt") as merged_handle:
                SeqIO.write(merger_generator(forward_handle,reverse_handle, rep_length), merged_handle, "fastq")
                
if __name__ == "__main__":
    main()

