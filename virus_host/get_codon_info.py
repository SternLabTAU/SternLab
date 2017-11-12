#! /usr/local/python_anaconda/bin/python3.4

from seqFileAnalyzer import get_codon_freqs

codon_data = get_codon_freqs("/sternadi/nobackup/volume1/talia_temp/virushostdb.cds.fna")
codon_info.to_csv("/sternadi/nobackup/volume1/talia_temp/codon_info_script.csv")
print("done")