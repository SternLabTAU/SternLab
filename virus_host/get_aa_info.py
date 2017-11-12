#! /usr/local/python_anaconda/bin/python3.4

from seqFileAnalyzer import get_amino_acid_freqs

aa_data = get_amino_acid_freqs("/sternadi/nobackup/volume1/talia_temp/virushostdb.cds.faa")
aa_data.to_csv("/sternadi/nobackup/volume1/talia_temp/aa_info_script.csv")
print("done")