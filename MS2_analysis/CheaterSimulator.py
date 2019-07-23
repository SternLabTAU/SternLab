import numpy as np
from itertools import permutations
import math
import pandas as pd
import matplotlib.pyplot as plt
import seaborn as sns
import math
from tqdm import tqdm
from scipy import stats
from scipy.stats import hypergeom
import os
import argparse



THRESHOLD = 10**-5
CT_BASE_FITNESS = 0
SYN_BASE_FITNESS = 0.1
WT_BASE_FITNESS = 1


def nCr(n,r):
    f = math.factorial
    return f(n) / f(r) / f(n-r)

class ChtrModel(object):
    """
    This is an implementation of the cheater dynamics model of bacteriophage MS2
    class ChtrModel holds all the initial parameters of the model.
    we add to this multiplayer model an additional player, which is a synonymous mutation with cheater like behaviour.
    """
    def __init__(self, N, n1, G, r, B, MOI=1):
        self.N = N
        self.n1 = n1
        self.genome = G
        self.r = r
        self.b = B
        self.moi = MOI
        self.cheater = []
        self.wt = []
        self.syn = []
        self.ecoli = []
        self.time = []



    def resistantN(self):
        return self.N * (1- self.r) # r is the factor of resistant bacteria



def sum_2_k_triplets(k):
    #generats a list with pairs (i, j) which sum to k for all 0<=q<=k

    options = list(range(k+1))
    all_triplets = list(permutations(options, 3)) + [(x,x,x) for x in options]
    triplets = [tri for tri in all_triplets if sum(tri) in options]

    return triplets


class Cycle(object):

    def __init__(self, N, n1, n2, moi, k):
        self.N = N
        self.n1 = n1
        self.n2 = n2
        self.moi = moi
        self.k = k
        self.payoffs = np.matrix([[1,1], [4,0]])

    def setN(self, N):
        self.N = N

    def set_n1(self, n1):
        self.n1 = n1

    def set_n2(self, n2):
        self.n2 = n2

    def set_n3(self, n3):
        self.n3 = n3

    def set_moi(self):
        self.moi = (self.n1 + self.n2) / self.N


    def poisson_prob(self,k1,k2):
        # according the poisson distribution calculate the probabilities
        return ((self.n1/self.N)**k1)*((self.n2/self.N)**k2)*\
               (math.exp(-(self.n1+self.n2)/self.N))/\
               (math.factorial(k1)*math.factorial(k2))


    def binomial_prob(self, k1, k2):
        # according the poisson distribution calculate the probabilities
        return stats.binom.pmf(k1, k1+k2, self.n1/(self.n1+self.n2)) * \
               stats.binom.pmf(k2, k1+k2, self.n2/(self.n1+self.n2))




    def normalize_payoff_rows(self):
        # this method normalizes the payoff matrix rows to sum to 1
        payoffs = self.payoffs
        W = np.zeros(shape=(2, 2))

        for i in range(2):
            sum_i = payoffs[i,].sum()
            for j in range(2):
                W[i,j] = payoffs[i,j] / sum_i

        return W


    def infection_proba(self, k, b):

        # init the fractions for each n_i
        n1 = 0
        n2 = 0

        pairs = sum_2_k_triplets(k)
        P = np.zeros(shape=(k+1, k+1))

        p_no_infection = np.exp(-self.moi)

        for pair in pairs:
            i = pair[0]
            j = pair[1]


            # add the value to the matrix

            P[i, j] = self.poisson_prob(i, j)

        # psum is the probability of infection. sum over all options and subtract the probability of no infection.
        psum = sum([P[i,j] for i in range(k+1) for j in range(k+1)]) - P[0,0]

        # if we have a 10^3 more viruses then cells the poission probability == zero!!
        if psum == 0:
            for pair in pairs:
                i = pair[0]
                j = pair[1]
                P[i, j] = self.binomial_prob(i, j)

            # probabilities should sum to 1.
            psum = sum([P[i, j] for i in range(k + 1) for j in range(k + 1) ]) - p_no_infection

        # get the relative fitness matrix for the next calculations
        # W = self.update_payoff()
        W = self.payoffs
        # W = self.normalize_payoff_rows()
        # W = self.payoffs_by_freq()

        # start filling the probabilities.
        for pair in pairs:
            i = pair[0]
            j = pair[1]

            # only cheater
            if i ==0 and j != 0:
                n2 += (P[i,j]/ psum) * W[1,1] * b * j


            # only wt infection, update n1 only
            if i != 0 and j == 0:
                n1 += (P[i, j] / psum) * W[0,0] * b


            # wt and cheater infection. wt and cheaters are coming out
            if i != 0 and j != 0:
                n1 += (P[i, j] / psum) * W[0, 1] * b * j/(j+i)
                n2 += (P[i, j] / psum) * W[1, 0] * b * i/(j+i)



        updated_n1 = n1 * self.N + self.n1 * (sum([P[0, j] / psum for j in range(k + 1) for q in range(k + 1)]) + p_no_infection)
        updated_n2 = n2 * self.N + self.n2 * (sum([P[i, 0] / psum for i in range(k + 1) for q in range(k + 1)]) + p_no_infection)


        return updated_n1, updated_n2



    def simulate_cycle(self, model, passage):

        # init first round
        self.setN(model.N * (1 - model.r))
        self.set_n1(model.n1)

        if model.cheater == []:
            model.cheater.append(self.n2)   # first iteration' self.n2 is defined

        else:
            dilution_factor = model.N / model.wt[-1]# in each passage we have a dilution factor.
            self.set_n2(model.cheater[-1] * dilution_factor + self.n1 *10 ** -5)
            model.cheater.append(self.n2)


        model.time.append(passage)
        model.ecoli.append(self.N)
        model.wt.append(self.n1)


        self.set_moi()

        print(self.N, self.n1, self.n2, self.moi)


        #while there are still cells to infect
        while self.N >= THRESHOLD:
            #update cheater, wt, ecoli and moi
            updatedN = (1 - model.r) * self.N * math.exp(-self.moi) * (2 ** 6)

            updated_n1, updated_n2 = self.infection_proba(self.k, model.b)

            if updatedN < THRESHOLD:
                model.ecoli.append(updatedN)
                model.wt.append(updated_n1)
                model.cheater.append(updated_n2)
                model.time.append(passage)

                break

            self.set_n2(updated_n2)
            self.setN(updatedN)
            self.set_n1(updated_n1)
            self.set_moi()

            # update models progression
            model.ecoli.append(self.N)
            model.wt.append(self.n1)
            model.cheater.append(self.n2)
            model.time.append(passage)



        print(self.N, self.n1, self.n2, self.moi)

def main(args):


    #create the model
    mdl = ChtrModel(N=args.N, n1=args.n1, G=args.genome, r=args.res, B=args.burst, MOI=args.expmoi)
    c = Cycle(N=mdl.N, n1=mdl.n1, n2=args.n2, moi=mdl.moi, k=args.k)

    # simulate cycles
    num_passages = args.passage
    for p in tqdm(range(1, num_passages + 1)):
        print("passage is : {}".format(p))
        c.simulate_cycle(mdl, p)



    df = pd.DataFrame({'Passage':mdl.time, 'N': mdl.ecoli, 'n1':mdl.wt, 'n2':mdl.cheater,})
    df['f1'] = df.apply(lambda row: row['n1'] / (row['n1'] + row['n2']), axis=1)
    df['f2'] = df.apply(lambda row: row['n2'] / (row['n1'] + row['n2']), axis=1)
    df.to_csv(args.out, index=False)

    sns.set_style('white')
    sns.pointplot(x='Passage', y='f1', data=df.drop_duplicates('Passage'), color='#4696AA')
    sns.pointplot(x='Passage', y='f2', data=df.drop_duplicates('Passage'), color='#DF1111', labels='T1764-')
    sns.despine(offset=10)


    plt.xlabel('Passage', fontsize=18)
    plt.ylabel('Cheater frequencies', fontsize=18)
    plt.ylim(0,1)
    plt.savefig(os.path.join(args.out, 'singlePlayer.png'), format='png', dpi=400,
                bbox_inches='tight')
    # plt.show()
    plt.gcf().clear()

if __name__ == '__main__':
    parser = argparse.ArgumentParser()
    parser.add_argument("-N", "--N", type=int, help="num of cells",default=10 ** 10)
    parser.add_argument("-n1", "--n1", type=int, help="wt particles", default=10 ** 10)
    parser.add_argument("-g", "--genome", type=int, help="genome lentgh", default=3569)
    parser.add_argument("-r", "--res", type=float, help="resistence fraction", default=0.3)
    parser.add_argument("-b", "--burst", type=int, help="burst size", default=10 ** 5)
    parser.add_argument("-m", "--expmoi", type=float, help="moi of the experiment", default=1)
    parser.add_argument("-n2", "--n2", type=int, help="initial number of cheaters", default=10**5)
    parser.add_argument("-k", "--k", type=int, help="maximal cell infected particles", default=15)
    parser.add_argument("-p", "--passage", type=int, help="num of passages to simulate", default=20)
    parser.add_argument("-o", "--out", type=str, help="path to save the output figure", required=True)

    args = parser.parse_args()
    main(args)





