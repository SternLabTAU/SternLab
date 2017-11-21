

import os.path
import sys
import os
import pbs_jobs
import argparse


def main(args):
	dest_dir_path = args.out_dir
	base_path = args.in_dir
	files = os.listdir(base_path)
	filenames = [f for f in files if "R1" in f and ".fastq" in f]

	full_file_paths = [os.path.join(base_path, filename) for filename in filenames]
	umerge_wrapper(full_file_paths, dest_dir_path)



def umerge_wrapper(file_lst, dest_dir):
	''' this method recieves a list of fastq files of the forward strand and for each file merges it with the corresponding reverse strand
	the merge is done by maoz_umerge_runner
	input - 
			file_lst - a list of fastq files, all R1 (forward)
			dest_dir - full path of the directory in which the merged files will be written to
	'''
	for filename in file_lst:
		maoz_umerge_runner(filename, os.path.join(dest_dir, os.path.basename(filename).split('.')[0]) + "_merge.fastq")


''' This is the merge command in python recieved by Maoz'''

def maoz_umerge_runner(input_file, output_file, alias = "umerge"):
	''' this method merges forward and reverse fastq files and assumes that R1 and R2 files exist
	input - input_file of the R1 (forward)
	output - merged file of R1 and R2 
	'''
	cmddfile = "umerge"; alias; tnum = 1; gmem = 2;
	cmds = "/usr/local/bin/usearch -fastq_nostagger -fastq_qmaxout 80 -fastq_qmax 80 -fastq_mergepairs %s -fastqout %s -report %s.report" % (input_file, output_file, output_file)
	pbs_jobs.create_pbs_cmd(cmddfile, alias, tnum, gmem, cmds)
	pbs_jobs.submit(cmddfile)


if __name__ == "__main__":
	parser = argparse.ArgumentParser()
	parser.add_argument("-i", "--in_dir", type=str, help="a path to an input directory in which R1 fastq files are saved (forward strand)", required=True)
	parser.add_argument("-o", "--out_dir", type=str, help="a path to an output directory in which the results will be saved", required=True)
	args = parser.parse_args()

	main(args)




