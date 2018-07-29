import pandas as pd
from sklearn.model_selection import train_test_split
from sklearn.linear_model import LinearRegression
from sklearn.metrics import mean_squared_error as mse, roc_curve, auc, accuracy_score, explained_variance_score, r2_score
import matplotlib
matplotlib.use('TkAgg')
import matplotlib.pyplot as plt
from sklearn.decomposition import PCA
import seaborn as sns
from tqdm import tqdm
from sklearn.feature_selection import RFE
import numpy as np




fits_path = r'/Volumes/STERNADILABHOME$/volume1/daniellem1/MS2/FITS/Mutation_rate_inference/mutation_rate_with_type'
ref = r'/Volumes/STERNADILABHOME$/volume1/shared/data/MS2/ms2_atcc_moran.fasta'
data_path = r'/Users/daniellemiller/Google Drive/Lab/Analysis/ms2_20170420/Mutation analysis/Di-nucleotide/freq_with_aa_170708_alternated_stops_with_type_dinucleotide_prev'


fits = pd.read_csv(fits_path, sep='\t')
data = pd.read_csv(data_path, sep='\t')
data = data[data['Mutation'].isin(['AG','GA','CT','TC'])]

with open(ref, 'r') as o:
    o.readline()    # read first ">..." line
    reference = o.read().replace('\n', '')

fits['ID'] = fits['Pos'].astype(str) + "_" + fits['Degree'].astype(str) + "_" +  fits['Replica'] + "_" +fits['Mutation']
data['ID'] = data['Pos'].astype(str) + "_" +  data['Degree'].astype(str) + "_" +  data['Replica'] + "_" + data['Mutation']


fits = fits[['ID', 'Mutation_rate']]

result = pd.merge(data, fits, on='ID')
result.drop(['ID', 'Sample', 'wtAA', 'mutAA', 'wtAA2', 'mutAA2', 'Replica'], axis=1, inplace=True)
#result.drop_duplicates(['ID'], inplace=True)


# apply ML model

for col in result.select_dtypes(exclude=['int', 'float']).columns:
    print(col)
    result[col] = pd.Categorical(result[col], categories=result[col].unique()).codes

# remove duplicates



explaind_var = []
r2 = []
test_mse = []
top_f = []
sim = []
for i in tqdm(range(1000)):

    train, test = train_test_split(result, test_size = 0.3)
    train_x = train.drop(['Mutation_rate'], axis=1)
    train_y = train['Mutation_rate']

    mdl = LinearRegression()
    mdl.fit(train_x, train_y)
    m = mdl.coef_

    rfe = RFE(mdl, n_features_to_select=2)
    x = result.drop(['Mutation_rate'], axis=1)
    r = rfe.fit(x, result['Mutation_rate'])
    top_f.extend(x.columns[r.support_].values)
    sim.extend([i] * len(x.columns[r.support_].values))


    test_x = test.drop(['Mutation_rate'], axis=1)
    test_y = test['Mutation_rate']

    pred = mdl.predict(test_x)
    explaind_var.append(explained_variance_score(test_y, pred))
    test_mse.append(mse(pred, test_y))
    r2.append(r2_score(test_y, pred))




# plot the results
df = pd.DataFrame({'MSE': test_mse, 'R_squared':r2, 'Explained_variance':explaind_var})
df_coeffs = pd.DataFrame({'Sim':sim, 'TOP':top_f})
print(df_coeffs)
df.boxplot()
plt.title("Regression results for mutation rate model")
plt.show()



# plot ROC curve and get AUC for classification problems

# get predicted values using predict_proba instead of just predict to get non linear fit of the roc curve
# probas_ = mdl.predict_proba(test_x)
#
#
# # Compute ROC curve and area the curve
# fpr, tpr, thresholds = roc_curve(test_y, probas_[:, 1])
# roc_auc = auc(fpr, tpr)
# plt.plot(fpr, tpr,lw=2, color='darkorange',
#          label='ROC fold (area = %0.2f)' % (roc_auc))
#
#
# plt.plot([0, 1], [0, 1], linestyle='--', lw=2, color='k',
#          label='Luck')
#
# plt.xlim([-0.05, 1.05])
# plt.ylim([-0.05, 1.05])
# plt.xlabel('False Positive Rate')
# plt.ylabel('True Positive Rate')
# plt.title('ROC curve Mutation rate prediction')
# #plt.savefig(r'/Users/daniellemiller/Dropbox/סדנא/GO/Ontotype/Validation/up_2_100_genes/with_all/ontotype_outputs/Figures/{} ROC curve Nonsense'.format(CANCER))
# plt.legend(loc="lower right")
# plt.show()