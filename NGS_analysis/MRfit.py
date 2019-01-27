import pandas as pd
import numpy as np
import os
import matplotlib.pyplot as plt
import seaborn as sns
from sklearn.model_selection import train_test_split
from sklearn.linear_model import LinearRegression, Lasso, LassoCV
from sklearn.linear_model import lasso_path, enet_path
from sklearn.metrics import mean_squared_error as mse, roc_curve, auc, accuracy_score, explained_variance_score, r2_score
from itertools import cycle

sns.set_style('white')


'''
Mutation rate fit by Lasso regression model.
we will use FITS outputs as labels and add important features'''

# this is fits runs for Mahoney WT1 results.
data = r'/Users/daniellemiller/Google Drive/FITS/quarter_allelic_summary.csv'
reference = r'/Users/daniellemiller/Google Drive/FITS/Mahoney.WT1.fasta'

# secondary structures by shape
shape = set(list(range(135,474)) + list(range(482,537)) + list(range(544,688))+\
        list(range(694,802)) + list(range(1779,1872)) + list(range(2327,2470))+\
        list(range(5276,5397)) + list(range(5715,5829)) + list(range(5900,5993))+\
        list(range(6180,6341)) + list(range(6354,6477)) + list(range(6929,7060))+\
        list(range(4443,4505)) + list(range(5742,5969)) + list(range(6902,7118)) +list(range(7132,7284)))

def create_input_matrix(data, reference):
    """
    This method create the input data matrix using an existing data frame, reference genome and additional information
    :param data: a data frame with mutation rtes estimations by FITS
    :param reference: reference file path
    :return: a matrix X and vector Y
    """

    df = pd.read_csv(data)
    df = df.dropna()
    df = df[df['Prob'] >= 0.95]

    with open(reference, 'r') as o:
        ref = o.read().replace('\n', '').split('genome')[-1]

    # add features
    df['Mutation_count'] = df['Freq'] * df['Read_count']
    df = df.drop(columns=['significance', 'Freq', 'Rank', 'Prob'])

    df['shape'] = df['Pos'].apply(lambda x: 0 if x in shape else 1)

    # add previous nucleotide context XY. notice that the first position is ref[-1] as it is not defined
    df['prev_nuc'] = df['Pos'].apply(lambda x: ref[x - 2])
    df['2prev_nuc'] = df['Pos'].apply(lambda x: ref[x - 3])
    df['next_nuc'] = df['Pos'].apply(lambda x: ref[x] if x != 7440 else 'g')
    df['2next_nuc'] = df['Pos'].apply(lambda x: ref[x + 1] if x != 7440 and x!= 7439 else 'g')


    # transform all columns which are not numbers into categorical data
    for col in df.select_dtypes(exclude=['int', 'float']).columns:
        df[col] = pd.Categorical(df[col], categories=df[col].unique()).codes


    return df



def fitModel(df):
    """
    fit the model to the requested matrix
    :param df: a data frame with all features
    :return:
    """

    train, test = train_test_split(df, test_size=0.3)
    train_x = train.drop(['MR'], axis=1)
    train_y = train['MR']




def coefs_paths(X,y):
    """
    plot the coefficients path by lasso regression
    :param X: train data in a data frame format with columns names
    :param y: labels vector data frame format
    :return: plots the lasso path
    """
    columns = list(X.columlns())
    X = np.array(X)
    y = np.array(y)

    X /= X.std(axis=0)  # Standardize data (easier to set the l1_ratio parameter)

    # Compute paths

    eps = 5e-3  # the smaller it is the longer is the path

    print("Computing regularization path using the lasso...")
    alphas_lasso, coefs_lasso, _ = lasso_path(X, y, eps, fit_intercept=False)

    print("Computing regularization path using the positive lasso...")
    alphas_positive_lasso, coefs_positive_lasso, _ = lasso_path(
        X, y, eps, positive=True, fit_intercept=False)
    print("Computing regularization path using the elastic net...")
    alphas_enet, coefs_enet, _ = enet_path(
        X, y, eps=eps, l1_ratio=0.8, fit_intercept=False)

    print("Computing regularization path using the positive elastic net...")
    alphas_positive_enet, coefs_positive_enet, _ = enet_path(
        X, y, eps=eps, l1_ratio=0.8, positive=True, fit_intercept=False)

    # Display results

    plt.figure(1)
    neg_log_alphas_lasso = -np.log10(alphas_lasso)
    neg_log_alphas_enet = -np.log10(alphas_enet)
    for coef_l, coef_e, c in zip(coefs_lasso, coefs_enet, columns):
        l1 = plt.plot(neg_log_alphas_lasso, coef_l)
        l2 = plt.plot(neg_log_alphas_enet, coef_e, linestyle='--')
        plt.text(neg_log_alphas_lasso[-1]+0.1, coef_l[-1],c, fontsize=10,color=l1[0].get_color())

    sns.despine(offset=10)
    plt.xlabel('-Log(alpha)')
    plt.ylabel('coefficients')
    plt.title('Lasso and Elastic-Net Paths')
    plt.legend((l1[-1], l2[-1]), ('Lasso', 'Elastic-Net'), loc='lower left')
    plt.axis('tight')

    plt.figure(2)
    neg_log_alphas_positive_lasso = -np.log10(alphas_positive_lasso)
    for coef_l, coef_pl, c in zip(coefs_lasso, coefs_positive_lasso, columns):
        l1 = plt.plot(neg_log_alphas_lasso, coef_l)
        l2 = plt.plot(neg_log_alphas_positive_lasso, coef_pl, linestyle='--')
        plt.text(neg_log_alphas_positive_lasso[-1] + 0.1, coef_l[-1], c, fontsize=10, color=l1[0].get_color())

    sns.despine(offset=10)
    plt.xlabel('-Log(alpha)')
    plt.ylabel('coefficients')
    plt.title('Lasso and positive Lasso')
    plt.legend((l1[-1], l2[-1]), ('Lasso', 'positive Lasso'), loc='lower left')
    plt.axis('tight')

    plt.figure(3)
    neg_log_alphas_positive_enet = -np.log10(alphas_positive_enet)
    for (coef_e, coef_pe, c) in zip(coefs_enet, coefs_positive_enet, columns):
        l1 = plt.plot(neg_log_alphas_enet, coef_e)
        l2 = plt.plot(neg_log_alphas_positive_enet, coef_pe, linestyle='--')
        plt.text(neg_log_alphas_positive_enet[-1] + 0.1, coef_e[-1], c, fontsize=10, color=l1[0].get_color())

    sns.despine(offset=10)
    plt.xlabel('-Log(alpha)')
    plt.ylabel('coefficients')
    plt.title('Elastic-Net and positive Elastic-Net')
    plt.legend((l1[-1], l2[-1]), ('Elastic-Net', 'positive Elastic-Net'),
               loc='lower left')
    plt.axis('tight')
    plt.show()