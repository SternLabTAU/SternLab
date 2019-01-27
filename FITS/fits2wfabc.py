# fits2wfabc.py
# converts fits data file to wfabc data file
# 2019-01-03 first version created - output to stdout

import sys
import pandas as pd


# load and pre-processing
def LoadFITSFile( my_filename, population_size ):

	# maybe define column data types with dtype?
	freqs_df = pd.read_table( my_filename )

	# take only the mutant
	freqs_df = freqs_df[ freqs_df["allele"] == 1 ];

	# we need copy number, not frequencies
	freqs_df["copy_number"] = freqs_df["freq"] * population_size
	freqs_df["copy_number"] = freqs_df["copy_number"].astype(int)

	# do we have a position column?
	if not "pos" in freqs_df.columns:
		freqs_df["pos"] = -1

	return(freqs_df)



def FITSData2WFABCData( fits_df, population_size ):

	gen_list = fits_df["gen"].unique()
	pos_list = fits_df["pos"].unique()

	num_generations = len(gen_list)
	num_positions = len(pos_list)

	sample_size_list = [str(population_size)] * num_generations

	print( str(num_positions) + " " + str(num_generations) )

	print( ",".join( map(str, gen_list) ) )

	for current_position in pos_list:
		position_df = fits_df[ fits_df["pos"] == current_position ]

		copy_number_list = position_df["copy_number"].tolist()

		if len(copy_number_list) != num_generations:
			print("Error: number of generations is " + str(num_generations) + " but copy number list size is " + str(len(copy_number_list )) )
			return

		print( ",".join( map(str, sample_size_list) ) )
		print( ",".join( map(str, copy_number_list) ) )


if len(sys.argv) < 3:
	print( "Syntax: fits2wfabc.py fits_data population_size" )
	quit(1)


data_filename = sys.argv[1]
pop_size = int(sys.argv[2])

fits_df = LoadFITSFile( data_filename, pop_size )

FITSData2WFABCData( fits_df, pop_size )
