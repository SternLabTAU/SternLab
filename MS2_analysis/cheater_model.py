import numpy as np
from itertools import permutations
import math
import pandas as pd
import matplotlib
matplotlib.use('TkAgg')
import matplotlib.pyplot as plt
import seaborn as sns
import math


GEN = 30     #number of generations
N = 10**9 # population size
n1 = 10**9     # WT population size (using plaque assay will be constant
n2_init = 100 # a start population size for the first iteration
BASE_WT_P = 2
EPS = 10**(-4)

def main():
    model = create_all_generations_data(n1, n2_init, N, 10)

    model['f1'] = model['n1']/(model['n1'] + model['n2'])
    model['f2'] = model['n2'] / (model['n1'] + model['n2'])
    model.to_csv(r'/Users/daniellemiller/Google Drive/Lab/Analysis/cheaters model/model_probs.csv')
    print(model)
    sns.pointplot(x='time', y='f2', data=model)
    plt.show()




def get_all_pairs(k):
    """
    generats a list with pirs (i, j) which sum to k for all 0<=q<=k
    :param k: an integer
    :return: a list of tuples (i,j)
    """
    options = list(range(k+1))
    all_pairs = list(permutations(options, 2)) + [(x,x) for x in options]
    pairs = [pair for pair in all_pairs if sum(pair) in options]

    return pairs



def get_probabilities(n1,n2,k1,k2):
    """
    calculates the probability of an infection as a condition of the MOI
    :param n1: WT population size
    :param n2: Parasite population size
    :param k1: number of WT
    :param k2: number of Parasite
    :return: probability according to poisson distribution
    """

    return ((n1/N)**k1)*((n2/N)**k2)*(math.exp(-(n1+n2)/N))/(math.factorial(k1)*math.factorial(k2))

def update(n1, n2, k, pairs, N):
    # updates a row in the data frame

    # create a matrix of zeros to be filled with the probabilities
    P = np.zeros(shape=(k+1,k+1))
    W = np.zeros(shape=(k+1,k+1))
    W_wt = np.zeros(shape=(k + 1, k + 1))


    for pair in pairs:
        i = pair[0]
        j = pair[1]

        # update poisson probabilities
        P[i,j] = get_probabilities(n1,n2,i,j)
        if (i == 0 and j == 0) or (i == 0 and j != 0) or (i != 0 and j/i >= 1.5):    # this one is undefined, there is no fitness for no infection at all
            W[i, j] = 0 + EPS
        elif j == 0 and i != 0:
            W[i, j] = 0 + EPS
        else:
            W[i,j] = BASE_WT_P*(math.exp(j**1.4))



    # get n2 value
    i = 0
    j = 0
    new_n2 = 0
    for i in range(k+1):
        for j in range(1,k+1):  # need Parasite to be != 0
            new_n2 += W[i,j]*P[i,j]

    return (new_n2 * N, P, W)



def create_all_generations_data(n1,n2_init, N, k):

    # create the basic data frame
    names = ['time', 'N', 'n1', 'n2'] + ['P{}-{}'.format(i,j) for i in range(k+1) for j in range(k+1)]
    model = {k:[] for k in names}
    n2 = n2_init

    pairs = get_all_pairs(k)

    for i in range(1, GEN+1):
        model['time'].append(i)
        model['N'].append(N)
        model['n1'].append(n1)
        model['n2'].append(n2)

        result = update(n1,n2, k, pairs, N)
        n2 = result[0]
        P = result[1]



        for j in range(k+1):
            for q in range(k+1):
                model['P{}-{}'.format(j,q)].append(P[j,q])


    # convert into a data frame
    print(model)
    model_df = pd.DataFrame.from_dict(model)
    return model_df




def get_pij_prob(n1,n2,max_moi):

    # num WT > parasite
    p10 = sum([get_probabilities(n1,n2,i,j) for i in range(1,max_moi+1) for j in range(max_moi+1) if i > j])

    # num parasite > num WT
    p01 = sum([get_probabilities(n1,n2,i,j) for i in range(max_moi+1) for j in range(1,max_moi+1) if i < j])

    # num WT = num parasite
    p11 = sum([get_probabilities(n1,n2,i,i) for i in range(1,max_moi+1)])

    # no infection
    p00 = get_probabilities(n1,n2,0,0)

    return {'p00': p00, 'p10':p10,'p01':p01, 'p11':p11}


if __name__ == "__main__":
    main()


