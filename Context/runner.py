import pandas
from pbs_runners import *
import os
from tqdm import tqdm
import matplotlib.pyplot as plt
import seaborn as sns
from sklearn.metrics import roc_curve, auc

betas_to_test = [0, 0.5, 1, 3, -1]
non_neutral_sites = [0, 10, 20, 50]
numSim = 100

running_folder = r'/sternadi/home/volume1/daniellem1/Context/paper/noise'

# motif - XXCGX --> XXTGX
# make directories for each one of possible combinations, if not exists

base_motif = 'P4__G_CT'

for beta in tqdm(betas_to_test):
    beta_dir = os.path.join(running_folder, 'beta_{}'.format(beta))
    os.makedirs(beta_dir, exist_ok=True)
    for num in tqdm(non_neutral_sites):
        curr_dir = os.path.join(beta_dir, '{}_{}'.format(base_motif, num))
        os.makedirs(curr_dir, exist_ok=True)

        # run the array
        cmd_noise = 'module load python/anaconda25_python-2.7.14\npython /sternadi/home/volume1/daniellem1/Context/code/main_noise.py\
         -f /sternadi/home/volume1/shared/data/ref_genomes/sabin2.full_genome.U882C.A2973G.C4905U.C5526U.fasta\
          -fr {} -c {} -n {} -b {} -ns 1'.format(os.path.join(curr_dir, 'sim_{}_$PBS_ARRAY_INDEX.freq'.format(num)),
                               os.path.join(curr_dir, 'sim_{}_regression_$PBS_ARRAY_INDEX.csv'.format(num)),
                                          float(num/100), beta)

        array_script_runner(cmd_noise, numSim, alias='FullAnalysis_beta_{}_p'.format(beta), load_python=False)

running_folder = r'/sternadi/home/volume1/daniellem1/Context/paper/not_noise'

for beta in tqdm(betas_to_test):
    beta_dir = os.path.join(running_folder, 'beta_{}'.format(beta))
    os.makedirs(beta_dir, exist_ok=True)
    for num in tqdm(non_neutral_sites):
        curr_dir = os.path.join(beta_dir, '{}_{}'.format(base_motif, num))
        os.makedirs(curr_dir, exist_ok=True)

        cmd_regular = 'module load python/anaconda25_python-2.7.14\npython /sternadi/home/volume1/daniellem1/Context/code/main_noise.py\
                 -f /sternadi/home/volume1/shared/data/ref_genomes/sabin2.full_genome.U882C.A2973G.C4905U.C5526U.fasta\
                  -fr {} -c {} -n {} -b {}'.format(os.path.join(curr_dir, 'sim_{}_$PBS_ARRAY_INDEX.freq'.format(num)),
                                                   os.path.join(curr_dir,
                                                                'sim_{}_regression_$PBS_ARRAY_INDEX.csv'.format(num)),
                                                   float(num / 100), beta)
        array_script_runner(cmd_regular, numSim, alias='FullAnalysis_beta_{}_p'.format(beta), load_python=False)


def create_roc_curve(xs, ys):
    colors = ['aqua', 'darkorange', 'cornflowerblue', 'deeppink']
    i = 0
    for key in xs:
        if 'beta_-1/' not in key:
            continue
        x = xs[key]
        y = ys[key]
        curr_auc = auc(x, y)
        plt.plot(x, y, lw=2,
                      label=key.split('CT_')[-1].split('/')[0] + '% non-neutral sites (area = {:.2f})'.format(curr_auc),
                 color=colors[i])
        plt.plot([0, 1], [0, 1], linestyle='--', lw=2, color='gray',
                                   label = 'Luck')
        plt.xlim([-0.05, 1.05])
        plt.ylim([-0.05, 1.05])
        plt.xlabel('False Positive Rate', fontsize=18)
        plt.ylabel('True Positive Rate', fontsize=18)
        plt.title('ROC curve', fontsize=20)
        sns.despine()
        plt.legend()
        plt.show()



def plot_pr_curve(xs, ys):

    colors = ['aqua', 'darkorange', 'cornflowerblue', 'deeppink']
    i = 0

    for key in xs:
        if 'beta_3/' not in key:
            continue
        precision = xs[key]
        recall = ys[key]
        ap = get_AP(recall, precision)
        plt.step(recall, precision, color=colors[i], alpha=0.2,
                          where = 'post', label = key.split('CT_')[-1].split('/')[
                                                            0] + '% non-neutral sites (AP={:.2f})'.format(ap))
        plt.fill_between(recall, precision, step='post', alpha=0.2,
                                        color = colors[i])
        i += 1

    plt.xlabel('Recall', fontsize=18)
    plt.ylabel('Precision', fontsize=18)
    plt.title('ROC curve ' + r'$\beta$ = -1', fontsize=20)
    sns.despine()
    plt.legend(loc="upper left")
    plt.show()
    plt.savefig(r'/Volumes/STERNADILABHOME$/volume1/daniellem1/Context/paper/tmp/precision-recall_B-1.png',
                     dpi=400, bbox_inches='tight')



