#! /usr/local/python_anaconda/bin/python3.4
from __future__ import print_function
import re, sys
from Bio import pairwise2, SeqIO, Seq
from Bio.pairwise2 import format_alignment

def align(seq,target):
    return pairwise2.align.localms(seq, target,10,-5,-10,-3)
    
def extract_primer_id(alignment, primer_id_start_idx, primer_id_end_idx):
    #print alignment[2:]
    #print "Score=" + str(alignment[2])
    #3print(format_alignment(*alignment))
    return alignment[0][alignment[3]+primer_id_start_idx:alignment[3]+primer_id_end_idx]

def primer_id_generator(fastq_handle, fasta, primer_id_start_idx, primer_id_end_idx):
    for seqRecord in SeqIO.parse(fastq_handle, "fastq"):
        toAlign=seqRecord.seq
        alignments = align(toAlign,fasta)
        if alignments[0][2] > 230:    
            yield extract_primer_id(alignments[0], primer_id_start_idx, primer_id_end_idx)
        alignments = align(toAlign.reverse_complement(),fasta)
        if alignments[0][2] > 230:        
            yield extract_primer_id(alignments[0], primer_id_start_idx, primer_id_end_idx)    

def main(args):


    primer_id_template=r'(N{8,})'
    input_fasta= args[0]
    input_fastq= args[1]
    output_file= args[2]
    record = SeqIO.read(input_fasta, "fasta")
    
    m = re.search(primer_id_template, str(record.seq))
    recZone = record.seq[max(0,m.start(0)-20):min(len(record.seq)-1,m.end(0)+20)]
    m = re.search(primer_id_template, str(recZone))
    
    with open(input_fastq, "r") as fastq_handle:
        primer_ids = primer_id_generator(fastq_handle, recZone, m.start(0), m.end(0))
        with open(output_file, "w") as output_handle:
            for primer_id in primer_ids:
                print (primer_id, file=output_handle)
        
if __name__ == "__main__":
    main(sys.argv[1:])