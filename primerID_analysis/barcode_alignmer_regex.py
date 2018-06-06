#! /usr/local/python_anaconda/bin/python3.4
from __future__ import print_function
import re, sys, io
from Bio import pairwise2, SeqIO, Seq
from Bio.pairwise2 import format_alignment

def align(seq,target):
    return pairwise2.align.localms(seq, target,10,-5,-10,-3)
    
def extract_primer_id(alignment, primer_id_start_idx, primer_id_end_idx):
    #print alignment[2:]
    #print "Score=" + str(alignment[2])
    print(format_alignment(*alignment))
    return alignment[0][alignment[3]+primer_id_start_idx:alignment[3]+primer_id_end_idx]

def primer_id_alignment_generator(fastq_handle, fasta, primer_id_start_idx, primer_id_end_idx):
    for seqRecord in SeqIO.parse(fastq_handle, "fastq"):
        toAlign=seqRecord.seq
        alignments = align(toAlign,fasta)
        if alignments[0][2] > 230:    
            yield extract_primer_id(alignments[0], primer_id_start_idx, primer_id_end_idx)
        alignments = align(toAlign.reverse_complement(),fasta)
        if alignments[0][2] > 230:        
            yield extract_primer_id(alignments[0], primer_id_start_idx, primer_id_end_idx)    

def primer_id_regex_generator(fastq_handle, sense_regex, antisense_regex):
    for seqRecord in SeqIO.parse(fastq_handle, "fastq"):
        toMatch=str(seqRecord.seq)
        sense_match = re.search(sense_regex, toMatch)
        if sense_match:
            yield sense_match.group(1)
            #yield "matched"
        else:
            antisense_match = re.search(antisense_regex, toMatch)
            if antisense_match:
                yield str(Seq.Seq(antisense_match.group(1)).reverse_complement())
                #yield antisense_match.group(1)#revcomp


def main(args):


    primer_id_template=r'(N{8,})'
    input_fasta= args[0]
    input_fastq= args[1]
    output_file= args[2]
    record = SeqIO.read(input_fasta, "fasta")
    
    m = re.search(primer_id_template, str(record.seq))
    recZone = record.seq[max(0,m.start(0)-20):min(len(record.seq)-1,m.end(0)+20)]
    surround=15
    primerIdZone=record.seq[max(0,m.start(0)-surround):min(len(record.seq)-1,m.end(0)+surround)]
    #print ("primer-ID area", primerIdZone)
    #print ("primer-ID length", len(primerIdZone)-2*surround)
    sense_template = str(primerIdZone[:15]) + r'(\w{15})' +str(primerIdZone[-15:])
    antisense_template= str(primerIdZone[-15:].reverse_complement()) + r'(\w{15})' +str(primerIdZone[:15].reverse_complement())
    #print (sense_template)
    #print (antisense_template)
    
    m = re.search(primer_id_template, str(recZone))
    
    with open(input_fastq, "r") as fastq_handle:
        #primer_ids = primer_id_alignment_generator(fastq_handle, recZone, m.start(0), m.end(0))
        primer_ids=primer_id_regex_generator(fastq_handle, sense_template, antisense_template)
        with io.open(output_file, "wb") as output_handle:
            for primer_id in primer_ids:
                print (primer_id, file=output_handle)
        
if __name__ == "__main__":
    main(sys.argv[1:])