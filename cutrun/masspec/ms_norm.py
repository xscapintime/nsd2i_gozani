import os, sys, glob
import pandas as pd
from itertools import product
from matplotlib import pyplot as plt
import matplotlib
import seaborn as sns

matplotlib.rc('pdf', fonttype=42)


## read ms table
msdat = pd.read_csv('miapaca2_ms.csv', sep=',', header=0, index_col=0)

## nsd2i devide dmso
comb = list(product(['DMSO','NSD2i'], ['D1','D5','D9'], ['Rep1', 'Rep2','Rep3']))

normed = pd.DataFrame()
for c in comb[9:]:
    # print(c[0])
    # print(c[1])

    normed[f'{c[1]}_{c[2]}'] = msdat[f'{c[0]}_{c[1]}_{c[2]}']/msdat[f'DMSO_{c[1]}_{c[2]}']


# normed.index = normed.index.str.replace('H33','H3')

## export 
normed.to_csv('absabun_fc.csv')


## heatmap of intereseted marks
# marks have C&R data
marks = ['H3K27ac', 'H3K27me3', 'H3K36me1', 'H3K36me2', 'H3K36me3', 'H3K4me1',
         'H3K4me3', 'H3K9me3']

mark_incr = normed.loc[[row in marks for row in normed.index],:]
mark_incr = pd.concat([mark_incr, normed.loc[['H33K27me3','H33K36me1','H33K36me2','H33K36me3'],:]])


plt.figure(figsize=(6,4))
ax = sns.heatmap(mark_incr, cmap="RdBu", mask=mark_incr.isnull(),
                 annot=True, fmt=".2f", linewidth=.5,
                 cbar_kws={'label': 'NSD2i/DMSO'},
                 center=0)
ax.set(xlabel="", ylabel="")
ax.set_facecolor('xkcd:gray')
plt.xticks(rotation=330,ha='left', va='top', rotation_mode='anchor')

plt.tight_layout()
plt.savefig('ms_fc.pdf')
plt.close()
