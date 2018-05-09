import pandas as pd
import numpy as np
import argparse
from Bio.Seq import Seq

''' This script adds to an existing data frame a column of Type which indicates the mutation type.
	Required columns:
	* Mutation
	* Pos
	also a reference genome is needed.

	Type column will include the following:
	* synonymous
	* non-synonymous
	* stop
	* read-through
'''


# coding regions

coding = list(range(130,1312)) + list(range(1335,1728))  + list(range(1761,3398))
#coding = list(range(1676,1902)) 
START = [130,1335, 1678, 3398]
#START = [1676]

def main(args):

	df = pd.read_csv(args.in_file)
	df['Mutation'] = df['Ref'] + df['Base']
	df = df[df.Ref != df.Base]
	#df = remove_dot_from_pos(df)
	df = df[df.Base != '-']
	#df["Pos"] = df["Pos"].astype(int)
	positions = df.Pos.values

	# get the reference genome
	with open(args.ref, 'r') as o:
		reference = o.read().replace('\n', '')

	types = []
	default = 'non-coding'
	bases = get_mutation_base_vector(df.Mutation.values)

	for i, pos in enumerate(positions):
		if pos not in coding:
			types.append(default)

		else:
			offset = get_offset(pos)
			#offset = 1676
			mut_type = get_mutation_type(pos, bases[i], reference, offset)
			types.append(mut_type)

	# add the new column
	df["Type"] = types

	# save data frame
	df.to_csv(args.out_file, encoding='utf-8', index = False)


def get_offset(pos):
	''' note: in this case there are overlaps. I will not consider the whole lysis protein'''
	if 130 <= pos <=1311:
		return 130
	elif 1335 <= pos <= 1729:
		return 1335
	elif 1761 <= pos <= 3398:
		return 1761



def get_mutation_base_vector(mutations):
	return [mut[1] for mut in mutations]

def get_mutation_type(pos, base, ref, offset=1):
	''' this method returns the mutation type after changing base in pos'''

	codon = (pos - offset + 1) % 3		# subtract offset from pos to fit it to AUG, translation start. 
	mut_type = ''
	mut_seq = ''
	
	if codon == 1:
		original = ref[pos - 1 : pos + 2]
		prot = Seq(original).transcribe().translate()
		mutant = Seq(base + original[1:]).transcribe().translate()
		mut_seq = Seq(base + original[1:]).transcribe()

	elif codon == 2:
		original = ref[pos - 2 : pos + 1]
		prot = Seq(original).transcribe().translate()
		mutant = Seq(original[0] + base + original[2]).transcribe().translate()
		mut_seq = Seq(original[0] + base + original[2]).transcribe()

	else:
		assert(codon == 0)
		original = ref[pos - 3 : pos]
		prot = Seq(original).transcribe().translate()
		mutant = Seq(original[0:2] + base).transcribe().translate()
		mut_seq = Seq(original[0:2] + base).transcribe()

	# check the type of the mutation
	if mutant == '*' and prot != '*' and mut_seq != 'UGA': 	# stop to stop will be a synonymous mutation, and so as read through
		mut_type = "stop"
	elif prot != '*' and mut_seq == 'UGA':
		mut_type = "read-through"
	elif prot == mutant:
		mut_type = "synonymous"
	elif prot != mutant:
		mut_type = "non-synonymous"
	return mut_type

def remove_dot_from_pos(data_frame):


    is_dot = (data_frame.Pos).astype(int) == (data_frame.Pos)
    data_frame["TMP"] = is_dot
    data_frame = data_frame[data_frame.TMP != False]
    data_frame.drop("TMP",axis=1, inplace=True)

    return data_frame



if __name__ == "__main__":
	parser = argparse.ArgumentParser()
	parser.add_argument("-i", "--in_file", type=str, help="a path to a data frame containing 'Mutation', 'Pos' columns", required=True)
	parser.add_argument("-r", "--ref", type=str, help="a path to a reference file", required=True)
	parser.add_argument("-o", "--out_file", type=str, help="a path to output file in which the result will be saved", required=True)
	parser.add_argument("-l", "--overlap", type=str, help="Y to include overlaps", required=False, default='N')
	args = parser.parse_args()
	main(args)