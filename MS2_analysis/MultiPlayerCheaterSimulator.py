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



THRESHOLD = 10**-5
CT_BASE_FITNESS = 0
SYN_BASE_FITNESS = 0
WT_BASE_FITNESS = 1
payoffs = np.matrix([[1,0.8,0.8], [2,0,0.2], [2,0.2,0.2]])

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

    def __init__(self, N, n1, n2, n3, moi, k):
        self.N = N
        self.n1 = n1
        self.n2 = n2
        self.n3 = n3
        self.moi = moi
        self.k = k

    def setN(self, N):
        self.N = N

    def set_n1(self, n1):
        self.n1 = n1

    def set_n2(self, n2):
        self.n2 = n2

    def set_n3(self, n3):
        self.n3 = n3

    def set_moi(self):
        self.moi = (self.n1 + self.n2 + self.n3) / self.N


    def poisson_prob(self,k1,k2,k3):
        # according the poisson distribution calculate the probabilities
        return ((self.n1/self.N)**k1)*((self.n2/self.N)**k2)*((self.n3/self.N)**k3)*\
               (math.exp(-(self.n1+self.n2+self.n3)/self.N))/\
               (math.factorial(k1)*math.factorial(k2)*math.factorial(k3))


    def binomial_prob(self, k1, k2, k3):
        # according the poisson distribution calculate the probabilities
        return stats.binom.pmf(k1, k1+k2+k3, self.n1/(self.n1+self.n2+self.n3)) * \
               stats.binom.pmf(k2, k1+k2+k3, self.n2/(self.n1+self.n2+self.n3)) * \
               stats.binom.pmf(k3, k1+k2+k3, self.n3/(self.n1+self.n2+self.n3)) * nCr(self.k, k1+k2+k3)


    def infection_proba(self, k, b):

        # init the fractions for each n_i
        n1 = 0
        n2 = 0
        n3 = 0
        triplets = sum_2_k_triplets(k)
        P = np.zeros(shape=(k+1, k+1, k+1))

        for triplet in triplets:
            i = triplet[0]
            j = triplet[1]
            q = triplet[2]

            # add the value to the matrix

            P[i, j, q] = self.poisson_prob(i, j, q)


        psum = sum([P[i,j,q] for i in range(k+1) for j in range(k+1) for q in range(k+1)])

        # if we have a 10^3 more viruses then cells the poission probability == zero!!
        if psum == 0:
            for triplet in triplets:
                i = triplet[0]
                j = triplet[1]
                q = triplet[2]
                P[i, j, q] = self.binomial_prob(i, j, q)

            # probabilities should sum to 1.
            psum = sum([P[i, j, q] for i in range(k + 1) for j in range(k + 1) for q in range(k + 1)])

        # start filling the probabilities.
        for triplet in triplets:
            i = triplet[0]
            j = triplet[1]
            q = triplet[2]


            # assuming syn has a working replicase both n2 and n3 are coming out after the infection
            if i ==0 and j != 0 and q != 0:
                n2 += (P[i,j,q]/ psum ) * q*payoffs[1,2] * b * j * (q/j)
                n3 += (P[i,j,q]/ psum ) * j*payoffs[2,1] * b * q

            # we have in infection of all three. we assume that this model is aditive
            if i != 0 and j != 0 and q != 0:
                n1 += (P[i, j, q] / psum) * ((j*payoffs[0, 1] + q*payoffs[0, 2]) / (j+q)) * b * i
                n2 += (P[i, j, q] / psum) * ((i*payoffs[1, 0] + q*payoffs[1, 2]) / (i+q)) * b * j * ((i+q)/j)
                n3 += (P[i, j, q] / psum) * ((i*payoffs[2, 0] + j*payoffs[2, 1]) / (i+j)) * b * q * (i/q)

            # only wt infection, update n1 only
            if i != 0 and j == 0 and q == 0:
                n1 += (P[i, j, q] / psum) * WT_BASE_FITNESS * b * i

            # wt and cheater infection. wt and cheaters are coming out
            if i != 0 and j != 0 and q == 0:
                n1 += (P[i, j, q] / psum) * j*payoffs[0, 1] * b * i
                n2 += (P[i, j, q] / psum) * i*payoffs[1, 0] * b * j * (i/j)

            # infection of wt and syn. this case we have wt coming out as well as syns.
            if i != 0 and j == 0 and q != 0:
                n1 += (P[i, j, q] / psum) * q*payoffs[0, 2] * b * i
                n3 += (P[i, j, q] / psum) * i*payoffs[2, 0] * b * q * (i/q)

        updated_n1 = n1 * self.N + self.n1 * (sum([P[0, j, q] / psum for j in range(k + 1) for q in range(k + 1)]) + P[0,0,0]/psum)
        updated_n2 = n2 * self.N + self.n2 * (sum([P[i, 0, q] / psum for i in range(k + 1) for q in range(k + 1)]) + P[0,0,0]/psum)
        updated_n3 = n3 * self.N + self.n3 * (sum([P[i, j, 0] / psum for j in range(k + 1) for i in range(k + 1)]) + P[0,0,0]/psum)

        addition = np.log((self.n1+self.n3)/self.n2)
        if addition < 0:
            addition = -1/addition
        updated_n2 = updated_n2 * addition

        return updated_n1, updated_n2, updated_n3



    def simulate_cycle(self, model, passage):

        # init first round
        self.setN(model.N * (1 - model.r))
        self.set_n1(model.n1)

        if model.cheater == []:
            model.cheater.append(self.n2)   # first iteration' self.n2 is defined

        if model.syn == []:
            model.syn.append(self.n3)    # first iteration' self.n3 is defined

        else:
            dilution_factor = model.N / model.wt[-1]# in each passage we have a dilution factor.

            self.set_n2(model.cheater[-1] * dilution_factor + self.n1 *10**-5)
            self.set_n3(model.syn[-1] * dilution_factor + self.n1 * 10 **-5)
            model.cheater.append(self.n2)
            model.syn.append(self.n3)

        model.time.append(passage)
        model.ecoli.append(self.N)
        model.wt.append(self.n1)


        self.set_moi()

        print(self.N, self.n1, self.n2, self.n3, self.moi)

        #while there are still cells to infect
        while self.N >= THRESHOLD:
            #update cheater, wt, ecoli and moi
            updatedN = (1 - model.r) * self.N * math.exp(-self.moi) * (2 ** 6)

            updated_n1, updated_n2, updated_n3 = self.infection_proba(self.k, model.b)

            if updatedN < THRESHOLD:
                model.ecoli.append(updatedN)
                model.wt.append(updated_n1)
                model.cheater.append(updated_n2)
                model.syn.append(updated_n3)
                model.time.append(passage)

                break

            self.set_n2(updated_n2)
            self.set_n3(updated_n3)
            self.setN(updatedN)
            self.set_n1(updated_n1)
            self.set_moi()

            # update models progression
            model.ecoli.append(self.N)
            model.wt.append(self.n1)
            model.cheater.append(self.n2)
            model.syn.append(self.n3)
            model.time.append(passage)



        print(self.N, self.n1, self.n2, self.n3, self.moi)


def main():

    out = r'/Users/daniellemiller/Google Drive/Lab/Analysis/cheaters model/multiPlayer_mdl_passage_cycle.csv'
    #create the model
    mdl = ChtrModel(N=10**9, n1=10**9, G=3569, r=0.2, B=10000, MOI=1)
    c = Cycle(N=mdl.N, n1=mdl.n1, n2=100, n3=100, moi=mdl.moi, k=10)

    num_passages = 20
    for p in tqdm(range(1, num_passages + 1)):
        print("passage is : {}".format(p))
        c.simulate_cycle(mdl, p)



    df = pd.DataFrame({'Passage':mdl.time, 'N': mdl.ecoli, 'n1':mdl.wt, 'n2':mdl.cheater, 'n3':mdl.syn})
    df['f1'] = df.apply(lambda row: row['n1'] / (row['n1'] + row['n2'] + row['n3']), axis=1)
    df['f2'] = df.apply(lambda row: row['n2'] / (row['n1'] + row['n2'] + row['n3']), axis=1)
    df['f3'] = df.apply(lambda row: row['n3'] / (row['n1'] + row['n2'] + row['n3']), axis=1)
    df.to_csv(out, index=False)

    sns.set_style('white')
    sns.pointplot(x='Passage', y='f2', data=df.drop_duplicates('Passage'), color='#3093C4')
    sns.despine(offset=10)
    # plt.title('r = {}; k = {}; w = {}'.format(mdl.r, c.k, BASE_FITNESS), fontsize=18, x=0.1)
    plt.xlabel('Passage', fontsize=18)
    plt.ylabel('Cheater frequencies', fontsize=18)
    plt.ylim(0,1)
    plt.savefig(r'/Users/daniellemiller/Google Drive/Lab/Analysis/cheaters model/multiPlayer.png', format='png', dpi=400,
                bbox_inches='tight')
    # plt.show()
    plt.gcf().clear()

if __name__ == '__main__':
    main()





