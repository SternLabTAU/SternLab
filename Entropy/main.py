

from OU_model_constructor import *


def main():


# run the MC simulation and wrute results to log.

    super_folder = r'/Volumes/STERNADILABHOME$/volume1/daniellem1/Entropy/data/Phylogeny/family/Asfarviridae'
    features = ['k5']
    out= r'/Volumes/STERNADILABHOME$/volume1/daniellem1/Entropy/data/OU_model/MonteCarlo'

    BM_OU_runner(super_folder, features, out, MC=True, log=True)





if __name__ == "__main__":
    main()