# This script calulates mutation rate by linear regression

import numpy as np
import matplotlib.pyplot as plt
import pandas as pd
import seaborn as sns
import os.path
import statsmodels.api as sm
from mutation_rate_utils import *
import argparse

######################## init key variables ###########################

all_mutations = ['AG', 'AC', 'AT', 'GA', 'GC', 'GT', 'CA', 'CG', 'CT', 'TG', 'TA', 'TC']
transitions = ['AG','GA','TC','CT']
transversions = [mut for mut in all_mutations if mut not in transitions]

limit = 2

THRESHOLD = 0.00001
#######################################################################

def main(args):

	df = pd.read_csv(args.in_file)
	lr_df = filter_freqs_2_regression(df)

	# init a vector to contain all slopes
	slopes = {}
	pvals = {}

	mutations = args.type
	if mutations == 1:
		mutations = transitions
	elif mutations ==2:
		mutations = transversions
	else:
		mutations = all_mutations
	
	# for each mutations type fit a linear regression model for each position in the genome and save the data to a file
	calculate_regression_slopes(df, slopes, p_values, mutation_type=mutations, limit=limit)
	all_slopes = pd.DataFrame.from_dict(slopes)

	# save the median vector to a file.
	if args.out_dir:
		all_slopes_df.to_csv(os.path.join(out, 'regression_slopes.csv'), encoding='utf-8', index = False)

	# plot the outputs
	#plot_all_slopes(averaged_slopes, replica, degree)
	#plot_only_transitions(averaged_slopes, replica, degree)


########################### Plots ##############################################

# TODO - adjust to all types


def plot_all_slopes(averaged_slopes, replica, degree):
	''' this method creats a joint plot of all mutation substitution rates'''

	order = ['AA', 'AC','AG','AT','CA','CC','CG','CT','GA','GC','GG','GT','TA','TC','TG','TT']

	# create a data frame to contain slopes and two intersection points (time 0 and 20)
	all_mut = all_mutations[:12]	# ignore all kinds
	averaged_slopes = averaged_slopes[:12]
	
	df1 = pd.DataFrame({"mut":all_mut, "slopes":averaged_slopes})
	df2 = df1.copy()

	# in each data frame create two columns (x,y) to represent a point in the mutation regression line
	df1["x"] = 0
	df2["x"] = 20 

	df1["y"] = TIME_0_RATE
	df2["y"] = df2.x * df2.slopes + TIME_0_RATE

	# create a third df of non mutations AA CC GG TT
	df3 = pd.DataFrame({"mut":['AA', 'CC','GG', 'TT'], "slopes":[0,0,0,0], "x":[0,0,0,0], "y":[0,0,0,0]})

	# unite all data frams and sort according to mutation type
	df = pd.concat([df1,df2,df3])
	df = df.sort_values(by="mut")
	df = df.reset_index()
	del df["index"]
	
	# create a slopes vector corresponds to mutations in order list
	slopes = [0] * len(order)
	for i, mut in enumerate(order):
		slopes[i] = df[df.mut == mut].slopes.head(n=1).iloc[0]


	# create the plot
	sns.set(style="ticks")
	grid = sns.FacetGrid(df, col="mut", hue="slopes", col_order=order, col_wrap=4, size=2)
	grid.set(xticks=np.arange(21), yticks=[0.000000001, 0.01],xlim=(0,21), ylim=(0.000000001, 0.01))
	grid.map(plt.plot, "x", "y", marker="o", ms=4)
	pylab.get_current_fig_manager().window.showMaximized()
	plt.yscale("log")

	# set titles
	for ax, title in zip(grid.axes.flat, slopes):
		ax.set_title("slope = " + str(title))

	sns.plt.show()



def plot_only_transitions(averaged_slopes, replica, degree):
	''' this method creats a joint plot of all mutation substitution rates.
	'''
	order = ['AG','CT','GA','TC']

	# create a data frame to contain slopes and two intersection points (time 0 and 15)
	all_mut = [all_mutations[0], all_mutations[3], all_mutations[8], all_mutations[11]]	# ignore all kinds
	averaged_slopes = [averaged_slopes[0], averaged_slopes[3], averaged_slopes[8], averaged_slopes[11]]
	
	df1 = pd.DataFrame({"mut":all_mut, "slopes":averaged_slopes})
	df2 = df1.copy()

	# in each data frame create two columns (x,y) to represent a point in the mutation regression line
	df1["passage"] = 0
	df2["passage"] = 20

	df1["frequency"] = TIME_0_RATE
	df2["frequency"] = df2.passage * df2.slopes + TIME_0_RATE


	# unite all data frams and sort according to mutation type
	df = pd.concat([df1,df2])
	df = df.sort_values(by="mut")
	df.reset_index(drop=True, inplace=True)
	
	
	# create a slopes vector corresponds to mutations in order list
	slopes = [0] * len(order)
	for i, mut in enumerate(order):
		slopes[i] = df[df.mut == mut].slopes.head(n=1).iloc[0]


	# create the plot
	sns.set(style="ticks")
	grid = sns.FacetGrid(df, col="mut", hue="slopes", col_order=order, col_wrap=4, size=2)
	grid.set(xticks=np.arange(21), yticks=[0.000001, 0.01],xlim=(0,21), ylim=(0.000001, 0.01))
	grid.map(plt.plot, "passage", "frequency", marker="o", ms=4)
	pylab.get_current_fig_manager().window.showMaximized()
	plt.yscale("log")

	# set titles
	#for ax, title in zip(grid.axes.flat, mut):
	#	ax.set_title(title)

	sns.plt.show()



if __name__ == "__main__":
	parser = argparse.ArgumentParser()
	parser.add_argument("-i", "--in_file", type=str, help="a path to an input file containing mutation data in a data frame format", required=True)
	parser.add_argument("-o", "--out_dir", type=str, help="a path to a directory in which the output files will be saved", default=None)
	parser.add_argument("-t", "--type", type=int, help="the type of mutations 1 for transitions, 2 for transversions, 3 for all ", default=1)



	args = parser.parse_args()

	main(args)
