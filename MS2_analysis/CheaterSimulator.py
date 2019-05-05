import numpy as np
from itertools import permutations
import pandas as pd
import matplotlib.pyplot as plt
import seaborn as sns
import math
from tqdm import tqdm
import argparse
import os


THRESHOLD = 10**-3
BASE_FITNESS = 2

class ChtrModel(object):
    """
    This is an implementation of the cheater dynamics model of bacteriophage MS2
    class ChtrModel holds all the initial parameters of the model.
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
        self.ecoli = []
        self.time = []


    def resistantN(self):
        return self.N * (1- self.r) # r is the factor of resistant bacteria



def sum_2_k_pairs(k):
    #generats a list with pairs (i, j) which sum to k for all 0<=q<=k

    options = list(range(k+1))
    all_pairs = list(permutations(options, 2)) + [(x,x) for x in options]
    pairs = [pair for pair in all_pairs if sum(pair) in options]

    return pairs


class Cycle(object):

    def __init__(self, N, n1, n2, moi, k):
        self.N = N
        self.n1 = n1
        self.n2 = n2
        self.moi = moi
        self.k = k

    def setN(self, N):
        self.N = N

    def set_n1(self, n1):
        self.n1 = n1

    def set_n2(self, n2):
        self.n2 = n2

    def set_moi(self):
        self.moi = (self.n1 + self.n2) / self.N
        #self.moi = (self.n1) / self.N


    def poisson_prob(self,k1,k2):
        # according the poisson distribution calculate the probabilities 
        return ((self.n1/self.N)**k1)*((self.n2/self.N)**k2)*(math.exp(-(self.n1+self.n2)/self.N))/\
               (math.factorial(k1)*math.factorial(k2))

    def infection_proba(self, k, b):

        # init the fraction n2
        n2 = 0
        pairs = sum_2_k_pairs(k)
        P = np.zeros(shape=(k+1, k+1))

        for pair in pairs:
            i = pair[0]
            j = pair[1]

            # add the value to the matrix
            P[i, j] = self.poisson_prob(i, j)

            # co-infection, i wt's j cheaters, 3.5 arbitrarily chosen to limit cheater fitness when there are much less wt
            if i != 0 and j != 0 :
                n2 += self.poisson_prob(i, j) * BASE_FITNESS  * max(1, (j+1/i)) * b * j #we have j cheaters, not one,
                # each has a burst size b
            # only cheater infection - there are cheaters coming out
            if i == 0 and j != 0:
                n2 += self.poisson_prob(i, j) * BASE_FITNESS * j
        #print("only cheater")
        #print(sum([P[i,j] for i in range(k+1) for j in range(k+1) if i==0 and j!=0]))

        return n2 * self.N




    def simulate_cycle(self, model, passage):

        # init first round
        self.setN(model.N * (1 - model.r))
        self.set_n1(model.n1)

        if model.cheater == []:
            model.cheater.append(self.n2)   # first iteration' self.n2 is defined
        else:
            dilution_factor = model.N / model.wt[-1]# in each passage we have a dilution factor.
            print(dilution_factor)
            if dilution_factor == 1:
                dilution_factor = 0.001
            self.set_n2(model.cheater[-1] * dilution_factor + self.n1 *10**-6)
            model.cheater.append(self.n2)

        model.time.append(passage)
        model.ecoli.append(self.N)
        model.wt.append(self.n1)


        self.set_moi()

        print(self.N, self.n1, self.n2, self.moi)

        while self.N >= THRESHOLD:
            #update cheater, wt, ecoli and moi
            updatedN = (1 - model.r) * self.N * math.exp(-self.moi) * (2 ** 6)
            if updatedN < THRESHOLD:
                break
            updated_n2 = self.infection_proba(self.k, model.b)
            # we have a wt population which do not infect. the breast size is for particles that infected an ecoli.
            updated_n1 = self.n1 * (1-math.exp(-self.n1/self.N)) * model.b + self.n1 * math.exp(-self.n1/self.N)

            self.set_n2(updated_n2)
            self.setN(updatedN)
            self.set_n1(updated_n1)
            self.set_moi()

            # update models progression
            model.ecoli.append(self.N)
            model.wt.append(self.n1)
            model.cheater.append(self.n2)
            model.time.append(passage)


        #
        # print("before while - N is {}".format(self.N))
        # updatedN = self.N
        # while updatedN >= THRESHOLD:
        #     print("updatedN {}".format(updatedN))
        #     updatedN = (1 - model.r) * self.N * math.exp(-self.moi) * (2**6)
        #     updated_n2 = self.infection_proba(self.k)
        #     print("self.N {}".format(self.N))
        #     self.setN(updatedN)
        #     print("self.N {}".format(self.N))
        #     self.set_n1(self.n1 * model.b)
        #     self.set_n2(updated_n2)
        #     print("moi now is : {}, n1={} N={}".format(self.moi, self.n1, self.N))
        #     self.set_moi()
        #
        #     # update model
        #     model.ecoli.append(self.N)
        #     model.wt.append(self.n1)
        #     model.cheater.append(self.n2)
        #     model.time.append(passage)
        # print("after while - N is {}".format(self.N))
        # # remove the last almost zero result from calculation
        # model.ecoli.pop()
        # model.wt.pop()
        # model.cheater.pop()
        # model.time.pop()



def main(args):


    #create the model
    mdl = ChtrModel(N=args.N, n1=args.n1, G=args.genome, r=args.res, B=args.burst, MOI=args.expmoi)
    c = Cycle(N=mdl.N, n1=mdl.n1, n2=args.n2, moi=mdl.moi, k=args.k)

    # simulate cycles
    num_passages = args.passage
    for p in tqdm(range(1, num_passages + 1)):
        print("passage is : {}".format(p))
        c.simulate_cycle(mdl, p)



    df = pd.DataFrame({'Passage':mdl.time, 'N': mdl.ecoli, 'n1':mdl.wt, 'n2':mdl.cheater})
    df['f1'] = df.apply(lambda row: row['n1'] / (row['n1'] + row['n2']), axis=1)
    df['f2'] = df.apply(lambda row: row['n2'] / (row['n1'] + row['n2']), axis=1)
    df.to_csv(os.path.join(args.out, 'dynamics_info.csv'), index=False)

    sns.set_style('white')
    sns.pointplot(x='Passage', y='f2', data=df.drop_duplicates('Passage'), color='#3093C4')
    sns.despine(offset=10)
    # plt.title('r = {}; k = {}; w = {}'.format(mdl.r, c.k, BASE_FITNESS), fontsize=18, x=0.1)
    plt.xlabel('Passage', fontsize=18)
    plt.ylabel('Cheater frequencies', fontsize=18)
    plt.ylim(0,1)
    plt.savefig(os.path.join(args.out, 'dynamics.png'), format='png', dpi=400,
                bbox_inches='tight')
    # plt.show()
    plt.gcf().clear()

if __name__ == '__main__':
    parser = argparse.ArgumentParser()
    parser.add_argument("-N", "--N", type=int, help="num of cells",default=10**9)
    parser.add_argument("-n1", "--n1", type=int, help="wt particles", default=10 ** 9)
    parser.add_argument("-g", "--genome", type=int, help="genome lentgh", default=3569)
    parser.add_argument("-r", "--res", type=float, help="resistence fraction", default=0.3)
    parser.add_argument("-b", "--burst", type=int, help="burst size", default=10 ** 5)
    parser.add_argument("-m", "--expmoi", type=float, help="moi of the experiment", default=1)
    parser.add_argument("-n2", "--n2", type=int, help="initial number of cheaters", default=1000)
    parser.add_argument("-k", "--k", type=int, help="maximal cell infected particles", default=10)
    parser.add_argument("-p", "--passage", type=int, help="num of passages to simulate", default=20)
    parser.add_argument("-o", "--out", type=int, help="path to save the output figure", default=20)

    args = parser.parse_args()
    main(args)





