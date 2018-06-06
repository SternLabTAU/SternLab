#! /usr/local/python_anaconda/bin/python3.4

# all primers from upac wildcards
# By Tal
# updated 2018-05-10

import sys
import argparse


def PropagatePrimers( original_sequence ):

	Regular_letters = ["A", "C", "G", "T"]

	iupac_dict = { 
	"R" : ["A", "G"], 
	"Y" : ["C", "T"],
	"S" : ["G", "C"],
	"W" : ["A", "T"],
	"K" : ["G", "T"],
	"M" : ["A", "C"],
	"B" : ["C", "G", "T"],
	"D" : ["A", "G", "T"],
	"H" : ["A", "C", "T"],
	"V" : ["A", "C", "G"],
	"N" : ["A", "C", "G", "T"]
	}

	current_seq_list = list()
	current_seq_list.append(original_sequence)



	# go through the orginal sequence, replace nucleotide codes with 
	# the nucleotides themselves
	for i in range(0, len(original_sequence)):

		if original_sequence[i] in Regular_letters:
			continue

		if original_sequence[i] not in iupac_dict.keys():
			raise Exception("Unkown character in primer (" + original_sequence + "): " + original_sequence[i])

		# code encountered - make replacements
		new_seq_list = list()

		to_replace_nucleotides = iupac_dict[ original_sequence[i] ]


		for current_seq in current_seq_list:

			for current_nucleotide in to_replace_nucleotides:

				tmp_seq = list(current_seq)
				tmp_seq[i] = current_nucleotide
				new_seq_list.append( "".join(tmp_seq) )

			current_seq_list = new_seq_list

	return current_seq_list
	

def PropagatePrimerPair( primer_F, primer_R, use_semicolon=False ):

	F_primer_list = PropagatePrimers(primer_F)
	R_primer_list = PropagatePrimers(primer_R)

	all_combinations = list()

	for current_f_primer_idx in range(0, len(F_primer_list)):
		for current_r_primer_idx in range(0, len(R_primer_list)):


			tmp_str = "primer_F" + str(current_f_primer_idx) + "\t" + F_primer_list[current_f_primer_idx] 
			if ( use_semicolon ):
				tmp_str = tmp_str + ";\tprimer_R" + str(current_r_primer_idx) + "\t" + R_primer_list[current_r_primer_idx]
			else:
				tmp_str = tmp_str + "\tprimer_R" + str(current_r_primer_idx) + "\t" + R_primer_list[current_r_primer_idx]
			all_combinations.append( tmp_str )

	return all_combinations


def main():

	if len(sys.argv) < 4:
		print( "Usage: ")

	parser = argparse.ArgumentParser(description="Populates all actual primers resulting from IUPAC wildcards")
	parser.add_argument( "-f", help="Forward primer" )
	parser.add_argument( "-r", help="Reverse primer" )
	parser.add_argument( "-target", help="formatting (thermo) default is none", nargs="?" )
	args = parser.parse_args()
	

	f_primer = args.f
	r_primer = args.r
	target = args.target


	if f_primer == None:
		raise Exception("No forward primer given")

	if r_primer == None:
		raise Exception("No reverse primer given")


	print( "Forward primer: " + f_primer )
	print( "Reverse primer: " + r_primer + "\n" )

	
	if target == "thermo":
		final_seq_list = PropagatePrimerPair( f_primer, r_primer, True )	
	else:
		final_seq_list = PropagatePrimerPair( f_primer, r_primer, False )

	for current_seq in final_seq_list:
			print(current_seq)


if __name__ == "__main__":
    main()

