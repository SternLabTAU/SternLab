We provide two ways to analyze the linkage between variants, or to identify the different strains in a population:
1. Split the population up into strains, where we define a strain as a combination of the occurence or lack of occurence of the varaints provided. 
2. For a provided list of variants, determine the linkage between every pair of variants. 

Both scripts require a list of variants to focus on. This list can either be created by AssociVar or provided by the user.

**Format for the csv with the mutations to check linkage for:** every mutation should have its own row, the csv should have a header row titled "variant", and the mutations should be written in the following format: "A1664.0G". This file is created automatically with the variant_association_test.py script in AssociVar. For example:

**variant**
------ |
A535.0G |
T1764.0-  |

### Strain Analysis


**usage:** 
strain_analysis.py [-h] -b INPUT_BLAST_DF -m INPUT_MUTATION_DF -p
                          INPUT_CHOSEN_MUTATIONS -o OUTPUT_FILE
                          [-f MINIMAL_MUTATION_FREQUENCY] -d
                          DELETION_ERROR_CUTOFF -s SUBSTITUTION_ERROR_CUTOFF

  -b INPUT_BLAST_DF, --input_blast_df INPUT_BLAST_DF
                        path to blasts df csv
			
  -m INPUT_MUTATION_DF, --input_mutation_df INPUT_MUTATION_DF
                        path to mutations df csv
			
  -p INPUT_CHOSEN_MUTATIONS, --input_chosen_mutations INPUT_CHOSEN_MUTATIONS
                        path to csv file with mutations to separate into
                        strains. Every mutation should have its own row, the
                        header row titled "variant", and the mutations should
                        be written in the following format: "A1664.0G". The
                        output file variants_chosen.csv from
                        association_test_variant.py can be used here.
			
  -o OUTPUT_FILE, --output_file OUTPUT_FILE
                        a path to an output file
			
  -f MINIMAL_MUTATION_FREQUENCY, --minimal_mutation_frequency MINIMAL_MUTATION_FREQUENCY
                        frequency cutoff for a single mutation. Only mutations
                        that are on the list and appear at least in this
                        frequency in the population will be included in the
                        strain analysis.
			
  -d DELETION_ERROR_CUTOFF, --deletion_error_cutoff DELETION_ERROR_CUTOFF
  
  -s SUBSTITUTION_ERROR_CUTOFF, --substitution_error_cutoff SUBSTITUTION_ERROR_CUTOFF


### Calculate Linkage
This script calculates the conditional probabilities for every mutation in list
to appear with any other mutation in the given list, including those probabilities 
for the WT variants for those positions. The script gets the blast dataframe path 
and the mutations dataframe path created by parse_blasts.py (see association_tests directory), 
a csv path with the mutations to check and a directory to save the results to.
   

**usage:** 
calculate_linkage.py [-h] -b INPUT_BLAST_DF -m INPUT_MUTATION_DF -p
                            INPUT_CHOSEN_MUTATION -o OUTPUT_FILE

  
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

						

**Output is interpreted as follows:**

Each cell represents the probability of the column when row is true (P(column|row)).
Row and column names represent either the mutation variant or the WT for that position. For example, "A535.0G" represents the mutation A to G at position 535, and "535" represents the WT base for that position.
It is recommended to open with excel and use Conditional Formatting > Color Scales to view this as presented in the article.
