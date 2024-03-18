import os 
import pandas as pd
from itertools import product
from matplotlib import pyplot as plt
import matplotlib
import seaborn as sns

matplotlib.rc('pdf', fonttype=42)


## read barcode stats
# Read the kmet
kmet = pd.read_csv('kmet_nuccount.txt', sep='\t', index_col=0)

# Read the kac file
kac = pd.read_csv('kac_nuccount.txt', sep='\t', index_col=0)

# Merge 
nuc_stats = pd.concat([kmet,kac]).dropna(axis=1)

# total barcode number
tot_barcode = pd.DataFrame(nuc_stats.sum())
tot_barcode['samples'] = tot_barcode.index.str.split('.').str[0]
tot_barcode = tot_barcode.groupby('samples').sum()

## unique reads number
uniq_stats = pd.read_csv('miapaca2_cr_bt2.clean.stats.txt',sep='\t',index_col=0)


## merge stats
stats_df = uniq_stats.copy()
stats_df['tot_barcode'] = tot_barcode


## scaling factor via nucloesome barcodes
# normalize by unique reads
# barcodes scaling factor
nuc_norm_hg = tot_barcode.apply(lambda x : uniq_stats['uniq_hg'][x.index] / x)
nuc_norm_hgec = tot_barcode.apply(lambda x : uniq_stats['uniq_all'][x.index] / x)

stats_df['nuc_norm_hg'] = nuc_norm_hg
stats_df['nuc_norm_hgec'] = nuc_norm_hgec


# export
nuc_norm_hg.to_csv('nuc_norm_hg.txt', sep='\t', header=False)
nuc_norm_hgec.to_csv('nuc_norm_hgec.txt', sep='\t', header=False)


## on-targt normalization
sum_ab = pd.DataFrame()

for i in range(0, len(nuc_stats.columns), 2):
    # Check if there's a next column to avoid IndexError
    if i+1 < len(nuc_stats.columns):
        nuc_stats.columns = nuc_stats.columns.str.replace('Ac','ac')
        nuc_stats.index = nuc_stats.index.str.replace('-','_')
        
        col_name =  nuc_stats.columns[i].split('.')[0]
        sum_ab[col_name] = nuc_stats.iloc[:, i] + nuc_stats.iloc[:, i+1]        
        
# Sum up every two rows
sum_ab = sum_ab.groupby(sum_ab.index.str.replace('_A$','').str.replace('_B$','')).sum()


## check antobody specificity
def on_target_ratio(df):
    new_df = df.copy()  # Create a copy of the input dataframe
    for column in new_df.columns:
        parts = column.split('_')
        if len(parts) >= 4 and parts[3] == 'IgG':
            column_sum = new_df[column].sum()
            new_df[column] = new_df[column] / column_sum
        else:
            targetsum = new_df[new_df.index.str.split('_').str[0] == parts[3]][column].sum()
            new_df[column] = new_df[column] / targetsum
    return new_df

ontar_norm = on_target_ratio(sum_ab)


# export
ontar_norm.to_csv('ontar_norm.txt', sep='\t', header=True)


## scale on target nucleosome barcodes
targetsum = {}
for column in sum_ab.columns:
    parts = column.split('_')
    targetsum['samples'] = targetsum.get('samples', []) + [column]

    if len(parts) >= 4 and parts[3] == 'IgG':
        targetsum['sum'] = targetsum.get('sum', []) + [sum_ab[column].sum()]
    else:
        targetsum['sum'] = targetsum.get('sum', []) + [sum_ab[sum_ab.index.str.split('_').str[0] == parts[3]][column].sum()]


targetsum_df = pd.DataFrame(targetsum, index=targetsum['samples'])
targetsum_df = targetsum_df[targetsum_df.columns.drop('samples')]

# targetsum_norm = targetsum_df.apply(lambda x : uniq_stats['uniq_hg'][x.index] / x[2])

stats_df.index = stats_df.index.str.replace('Ac','ac')

stats_df['target_barcode'] = targetsum_df
stats_df['target_norm_hg'] = stats_df['uniq_hg'] /  stats_df['target_barcode'] 
stats_df['target_norm_hgec'] = stats_df['uniq_all']/ stats_df['target_barcode']

stats_df = pd.concat([stats_df, sum_ab.T], axis=1)

stats_df.to_csv('barcode_stats_df.txt', sep='\t', header=True)

targetsum_norm = stats_df[['target_norm_hg']]
targetsum_norm = targetsum_norm.dropna()
# targetsum_norm.index = targetsum_norm.index.str.replace('Ac','ac')

targetsum_norm.to_csv('tarnuc_norm_hg.txt', sep='\t', header=False)


#### NSD2i vs Control
## human vs all barcodes
# long to wide
nuc_norm_hg.index = nuc_norm_hg.index.str.replace('Ac','ac')
nuc_norm_hg['mark'] = nuc_norm_hg.index.str.split('_').str[3]
nuc_norm_hg['timerep'] = [ '_'.join(n) for n in  nuc_norm_hg.index.str.split('_').str[0:3] ]
nuc_norm_hg_wide = nuc_norm_hg.pivot_table(index='mark', columns='timerep', values=0)

# norm
comb = list(product(['D1','D5','D9'], ['Rep1', 'Rep2']))

allbar_normed = pd.DataFrame()
for c in comb:
    # print(c[0])
    # print(c[1])

    allbar_normed[f'{c[0]}_{c[1]}'] = nuc_norm_hg_wide[f'NSD2i_{c[0]}_{c[1]}']/nuc_norm_hg_wide[f'Vehicle_{c[0]}_{c[1]}']


plt.figure(figsize=(6,4))
ax = sns.heatmap(allbar_normed, cmap="RdBu", mask=allbar_normed.isnull(),
                 annot=True, fmt=".2f", linewidth=.5,
                 cbar_kws={'label': 'NSD2i/DMSO'},
                 center=0)
ax.set(xlabel="", ylabel="")
ax.set_facecolor('xkcd:gray')
plt.xticks(rotation=330,ha='left', va='top', rotation_mode='anchor')

plt.tight_layout()
plt.savefig('allbar_normed.pdf')
plt.close()


# also norm to IgG
allbar_normed_igg = allbar_normed.iloc[:-1,:]/allbar_normed.iloc[-1,:]

plt.figure(figsize=(6,4))
ax = sns.heatmap(allbar_normed_igg, cmap="RdBu", mask=allbar_normed_igg.isnull(),
                 annot=True, fmt=".2f", linewidth=.5,
                 cbar_kws={'label': 'NSD2i/DMSO'},
                 center=0)
ax.set(xlabel="", ylabel="")
ax.set_facecolor('xkcd:gray')
plt.xticks(rotation=330,ha='left', va='top', rotation_mode='anchor')

plt.tight_layout()
plt.savefig('allbar_normed_igg.pdf')
plt.close()



## human vs target barcodes
# long to wide
targetsum_norm['mark'] = targetsum_norm.index.str.split('_').str[3]
targetsum_norm['timerep'] = [ '_'.join(n) for n in targetsum_norm.index.str.split('_').str[0:3] ]
targetsum_df_wide = targetsum_norm.pivot_table(index='mark', columns='timerep', values='target_norm_hg')

# norm
targetbar_normed = pd.DataFrame()
for c in comb:
    # print(c[0])
    # print(c[1])

    targetbar_normed[f'{c[0]}_{c[1]}'] = targetsum_df_wide[f'NSD2i_{c[0]}_{c[1]}']/targetsum_df_wide[f'Vehicle_{c[0]}_{c[1]}']


plt.figure(figsize=(6,4))
ax = sns.heatmap(targetbar_normed, cmap="RdBu", mask=targetbar_normed.isnull(),
                 annot=True, fmt=".2f", linewidth=.5,
                 cbar_kws={'label': 'NSD2i/DMSO'},
                 center=0)
ax.set(xlabel="", ylabel="")
ax.set_facecolor('xkcd:gray')
plt.xticks(rotation=330,ha='left', va='top', rotation_mode='anchor')

plt.tight_layout()
plt.savefig('tarbar_normed.pdf')
plt.close()



# also norm to IgG
targetbar_normed_igg = targetbar_normed.iloc[:-1,:]/targetbar_normed.iloc[-1,:]

plt.figure(figsize=(6,4))
ax = sns.heatmap(targetbar_normed_igg, cmap="RdBu", mask=targetbar_normed_igg.isnull(),
                 annot=True, fmt=".2f", linewidth=.5,
                 cbar_kws={'label': 'NSD2i/DMSO'},
                 center=0)
ax.set(xlabel="", ylabel="")
ax.set_facecolor('xkcd:gray')
plt.xticks(rotation=330,ha='left', va='top', rotation_mode='anchor')

plt.tight_layout()
plt.savefig('tarbar_normed_igg.pdf')
plt.close()



