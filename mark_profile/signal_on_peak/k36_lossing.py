import os,sys
import numpy as np
import pandas as pd
from itertools import product
import seaborn as sns


sigmat = pd.read_csv('mark_profile.sig.mat.txt', header=0, index_col=0, sep='\t')

## nsd2i devide dmso
marks = sigmat.columns.str.split('_').str[3].unique()
comb = list(product(['NSD2i'], ['D1','D5','D9'], ['Rep1', 'Rep2'], marks))


## log2 fold change
# log2((x+1)/(y+1))

normed = pd.DataFrame()
for c in comb:
    # print(c[0])
    # print(c[1])

    normed[f'{c[1]}_{c[2]}_{c[3]}'] = np.log2((sigmat[f'{c[0]}_{c[1]}_{c[2]}_{c[3]}']+1)/(sigmat[f'Vehicle_{c[1]}_{c[2]}_{c[3]}']+1))



## k36me2 lossing
log2k36me2 = normed.filter(like='K36me2')


# filter reps with diff. signs
log2k36me2 = log2k36me2[(log2k36me2.groupby(log2k36me2.columns.str.split('_').str[0], axis=1).prod() > 0).sum(axis=1) == 3]


# merge reps
log2k36me2_merged  = log2k36me2.groupby(log2k36me2.columns.str.split('_').str[0], axis=1).mean()


## export
log2k36me2_merged.to_csv('k36me2_changes.txt', sep='\t')
