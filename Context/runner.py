import pandas
from pbs_runners import *
import os

betas_to_test = [0, 0.5, 1, 3, -1]
non_neutral_sites = [0, 10, 20, 50]
numSim = 100

running_folder = r'/sternadi/home/volume1/daniellem1/Context/paper/noise'

# motif - XXCGX --> XXTGX
# make directories for each one of possible combinations, if not exists

base_motif = 'P4__G_CT'

for beta in betas_to_test:
    beta_dir = os.path.join(running_folder, 'beta_{}'.format(beta))
    os.makedirs(beta_dir, exist_ok=True)
    for num in non_neutral_sites:
        curr_dir = os.path.join(beta_dir, '{}_{}'.format(base_motif, num))
        os.makedirs(curr_dir, exist_ok=True)

        # run the array
        cmd_noise = 'module load python/anaconda25_python-2.7.14\npython /sternadi/home/volume1/daniellem1/Context/code/main_noise.py\
         -f /sternadi/home/volume1/shared/data/ref_genomes/sabin2.full_genome.U882C.A2973G.C4905U.C5526U.fasta\
          -fr {} -c {} -n {} -b {} -n 1'.format(os.path.join(curr_dir, 'sim_{}_$PBS_ARRAY_INDEX.freq'.format(num)),
                               os.path.join(curr_dir, 'sim_{}_regression_$PBS_ARRAY_INDEX.csv'.format(num)),
                                          float(num/100), beta)

        array_script_runner(cmd_noise, numSim, alias='FullAnalysis_beta_{}_p', load_python=False)

running_folder = r'/sternadi/home/volume1/daniellem1/Context/paper/not_noise'

for beta in betas_to_test:
    beta_dir = os.path.join(running_folder, 'beta_{}'.format(beta))
    os.makedirs(beta_dir, exist_ok=True)
    for num in non_neutral_sites:
        curr_dir = os.path.join(beta_dir, '{}_{}'.format(base_motif, num))
        os.makedirs(curr_dir, exist_ok=True)

        cmd_regular = 'module load python/anaconda25_python-2.7.14\npython /sternadi/home/volume1/daniellem1/Context/code/main_noise.py\
                 -f /sternadi/home/volume1/shared/data/ref_genomes/sabin2.full_genome.U882C.A2973G.C4905U.C5526U.fasta\
                  -fr {} -c {} -n {} -b {}'.format(os.path.join(curr_dir, 'sim_{}_$PBS_ARRAY_INDEX.freq'.format(num)),
                                                   os.path.join(curr_dir,
                                                                'sim_{}_regression_$PBS_ARRAY_INDEX.csv'.format(num)),
                                                   float(num / 100), beta)
        array_script_runner(cmd_noise, numSim, alias='FullAnalysis_beta_{}_p', load_python=False)










