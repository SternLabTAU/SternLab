import pandas as pd
import argparse
import FITS.fitness_utilities
import scipy.stats as stats
import seaborn as sns
import matplotlib.pyplot as plt

sns.set_context("poster")

primers = list(range(1, 20)) + list(range(1291, 1304)) + list(range(1179, 1200)) + list(range(2270, 2288)) +\
                list(range(2167, 2188)) + list(range(3548, 3570))
problematic = [17, 18, 19, 20, 21, 22, 23, 183, 188, 189, 190, 196, 274, 317, 364, 452, 762, 2719, 3117, 3133, 3139,
                   3143, 3146, 3150, 3401, 3539, 3542]

proteins = {"A_protein":list(range(130,1311)), "CP_protein":list(range(1335,1727)),\
			 "lysis_protein":list(range(1678,1905)),"replicase_protein":list(range(1761,3398))}


def main(args):

	df = pd.read_csv(args.in_file, sep='\t')
	rep1 = 'A'
	rep2 = 'B'
	all_types = ['AG', 'GA', 'TC', 'CT']

	# remove problematic positions
	df = df[~df["Pos"].isin(primers + problematic)]

	# remove the ? from category and ambiguous positions (replicas do not agree on category)
	df = FITS.fitness_utilities.remove_sign(df)
	df = FITS.fitness_utilities.filter_replicas_ambiguity(df)
	df.reset_index(drop=True, inplace=True)

	# plot correlation - arrays should be from the same length
	FITS.fitness_utilities.plot_correlation(df)

	# plot mutation types by degree (syn non-sun and stop)
	FITS.fitness_utilities.plot_mutation_types_diff(df, 37, rep1)
	FITS.fitness_utilities.plot_mutation_types_diff(df, 37, rep2)
	FITS.fitness_utilities.plot_mutation_types_diff(df, 41, rep1)
	FITS.fitness_utilities.plot_mutation_types_diff(df, 41, rep2)

	# plot fitness values according mutation type (transitions)
	FITS.fitness_utilities.plot_mutations_diff(df, 37, rep1)
	FITS.fitness_utilities.plot_mutations_diff(df, 37, rep2)
	FITS.fitness_utilities.plot_mutations_diff(df, 41, rep1)
	FITS.fitness_utilities.plot_mutations_diff(df, 41, rep2)

	# plot fitness values boxplot according degree
	FITS.fitness_utilities.plot_degree_boxplot(df)

	# plot fitness value heat map according to replica
	#FITS.fitness_utilities.plot_heatmap(df, rep1)
	#FITS.fitness_utilities.plot_heatmap(df, rep2)

	# for each protein plot its dfe and discrete dfe

	for protein in proteins:
		mut_df = df[df['Pos'].isin(proteins[protein])]
		FITS.fitness_utilities.plot_dfe(mut_df, protein)
		df_a = mut_df[(mut_df.Replica == 'A')]
		FITS.fitness_utilities.discretize_fitness_dfe(df_a, "{} Replica A".format(protein))
		df_b = mut_df[(mut_df.Replica == 'B')]
		FITS.fitness_utilities.discretize_fitness_dfe(df_b, "{} Replica B".format(protein))

	# for each mutation type plot its dfe and discrete dfe

	for t in all_types:
		mut_df = df[df['Mutation'] == t]
		FITS.fitness_utilities.plot_dfe(mut_df, t)
		df_a = mut_df[(mut_df.Replica == 'A')]
		FITS.fitness_utilities.discretize_fitness_dfe(df_a, "Mutation type {} Replica A".format(t))
		df_b = mut_df[(mut_df.Replica == 'B')]
		FITS.fitness_utilities.discretize_fitness_dfe(df_b, "Mutation type {} Replica B".format(t))

	# discrete DFE for each replica

	df_a = df[(df.Replica == 'A')]
	FITS.fitness_utilities.discretize_fitness_dfe(df_a, "Replica A")

	df_b = df[(df.Replica == 'B')]
	FITS.fitness_utilities.discretize_fitness_dfe(df_b, "Replica B")

	#FITS.fitness_utilities.test_kolmogorov_smirnov(df, rep1, rep2, 37, 41, mut_type)




if __name__ == "__main__":
	parser = argparse.ArgumentParser()
	parser.add_argument("-i", "--in_file", type=str, help="a path to a fits output file containing fitness values ", required=True)
	args = parser.parse_args()
	main(args)