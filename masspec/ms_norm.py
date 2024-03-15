######### nomalizing to the least ##########


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


## calculate scaling facotr for cut & run
sf_df = msdat.loc[[row in marks for row in normed.index],~msdat.columns.str.contains('Rep3')]

ms_long =  sf_df.reset_index().melt(id_vars='index', var_name='Sample', value_name='ms')

sf_df = sf_df.apply(lambda x: x/x.min(), axis=1)
sf_df_long = sf_df.reset_index().melt(id_vars='index', var_name='Sample', value_name='Scaling Factor')
sf_df_long['sample'] = sf_df_long['Sample'].str.replace('DMSO','Vehicle') + \
    '_' + sf_df_long['index'].str.replace('ac','Ac')

sf_df_long = sf_df_long.merge(ms_long, on=['index','Sample'], how='inner')

sf_df_long.to_csv('ms_scaling_factor.txt', index=False, sep='\t', header=False)


## calculate subsampling factor for cut & run
## unique reads number
# uniq_stats = pd.read_csv('miapaca2_cr_bt2.clean.stats.txt',sep='\t',index_col=0)

# subfactor_df = pd.merge(sf_df_long, uniq_stats[['uniq_hg']], left_on='sample', right_index=True, how='inner')

# subfactor_df['expcted_reads'] = subfactor_df.apply(lambda row: row['uniq_hg'] if row['Scaling Factor'] == 1 \
#                                                    else row['Scaling Factor'] * subfactor_df.loc[(subfactor_df['index'] == row['index']) & (subfactor_df['Scaling Factor'] == 1), 'uniq_hg'].values[0], axis=1)

# subfactor_df['subsample_fraction'] = subfactor_df['expcted_reads'] / subfactor_df['uniq_hg']

# subfactor_df.to_csv('ms_subsample_factor.txt', index=False, sep='\t', header=False)