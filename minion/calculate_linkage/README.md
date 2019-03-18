This script calculates the conditional probabilities for every mutation in list
to appear with any other mutation in the given list, including those probabilities 
for the WT variants for those positions. The script gets the blast dataframe path, 
and the mutations dataframe path created by parse_blasts.py (see association_tests directory), 
a csv path with the mutations to check and a directory to save the results to.
   

usage: calculate_linkage.py [-h] -b INPUT_BLAST_DF -m INPUT_MUTATION_DF -p
                            INPUT_CHOSEN_MUTATION -o OUTPUT_FILE

optional arguments:
  -h, --help            show this help message and exit
  -b INPUT_BLAST_DF, --input_blast_df INPUT_BLAST_DF
                        path to blasts df csv
  -m INPUT_MUTATION_DF, --input_mutation_df INPUT_MUTATION_DF
                        path to mutations df csv
  -p INPUT_CHOSEN_MUTATION, --input_chosen_mutation INPUT_CHOSEN_MUTATION
                        path to csv file with mutations to check linkage of.
                        Every mutation should have its own row, no header row,
                        and the mutations should be written in the following
                        format: "A1664.0G"
  -o OUTPUT_FILE, --output_file OUTPUT_FILE
                        a path to an output file

						
						
Format for the csv with the mutations to check linkage for: every mutation should 
have its own row, the csv should have no header row, and the mutations should be written 
in the following format: "A1664.0G". For example:
"A1664.0G
A535.0G
T1764.0-"