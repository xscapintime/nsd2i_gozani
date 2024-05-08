import os,glob
import pandas as pd
import numpy as np
from scipy import stats
from umap.umap_ import UMAP
import seaborn as sns
import matplotlib
from matplotlib import pyplot as plt
import matplotlib.font_manager
import matplotlib.cm as cm
import matplotlib.colors as mcolors
from Ensembl_converter import EnsemblConverter



plt.style.use('seaborn-poster')
matplotlib.rcParams['pdf.fonttype'] = 42
matplotlib.rcParams['ps.fonttype'] = 42
matplotlib.rcParams['font.family'] = 'sans-serif'
matplotlib.rcParams['font.sans-serif'] = 'Arial'
matplotlib.rcParams.update({'font.size': 6})

## load ensemble id to symbol matching table
ensembl_syb = pd.read_csv('human_ensembl_syb.tsv', header=0, index_col=None, sep='\t')
ensembl_syb = dict(zip(ensembl_syb['ensembl_gene_id'], ensembl_syb['hgnc_symbol']))


## load tpm table
tpm = pd.read_csv('../tpm/tximport-tpm.csv', header=0, index_col=0)

## change index to gene symbol
tpm = tpm.set_index(tpm.index.map(ensembl_syb))
tpm = tpm[tpm.index.notnull()]

## filter lowly expressed transcripts
# half of the samples should have TPM > 1
tpm = tpm.loc[(tpm >= 1).sum(axis=1) > tpm.shape[1]/2]


## mean acorss replicates
tpm_ave = tpm.groupby(tpm.columns.str[:-1], axis=1).mean()


## zscore
tpm_zscore = stats.zscore(tpm_ave, axis=1,ddof=1)
tpm_zscore = tpm_zscore[['Day1_N', 'Day5_N', 'Day9_N', 'Day1_C', 'Day5_C', 'Day9_C']]

## normalize
# tpm_N = tpm_ave.filter(like='N') + 0.1
# tpm_C = tpm_ave.filter(like='C') + 0.1
# normed =  np.log2(np.divide(*tpm_N.align(tpm_C, axis=0)))

# normed_z = stats.zscore(tpm_ave.filter(like='N') - tpm_ave.filter(like='C').values, axis=1, ddof=1)


#### UMAP ####
# seed
np.random.seed(1024)

# UMAP demension redunction on tpm of all samples
umap = UMAP(n_neighbors=20, min_dist=0.25)


# embedding
embedding = pd.DataFrame(umap.fit_transform(tpm), columns = ['UMAP1','UMAP2'])
embedding.index = tpm.index



########## viz #######
# 0 as center
vcenter = 2
# vmin, vmax = fcstat.iloc[:,:171].min().min(), fcstat.iloc[:,:171].max().max()
vmin, vmax = (0,4)
normalize = mcolors.TwoSlopeNorm(vcenter=vcenter, vmin=vmin, vmax=vmax)

groups = tpm_ave.columns

for g in groups:
    plt.figure(figsize=(5.5,4))

    ax = sns.scatterplot(x='UMAP1', y='UMAP2', data=embedding,
                    hue=np.log10(tpm_ave[g]+1).to_list(),
                    linewidth=0, s=3, alpha=.7, edgecolor=None,
                    palette="viridis",
                    hue_norm=normalize)
    
    
    # norm = plt.Normalize(fcstat.iloc[:,:171].min().min(), fcstat.iloc[:,:171].max().max())
    sm = plt.cm.ScalarMappable(cmap="viridis", norm=normalize)
    sm.set_array([])
    
    # Remove the legend and add a colorbar
    ax.get_legend().remove()
    ax.figure.colorbar(sm)
    
    # plt.show()
    plt.tight_layout()
    plt.savefig(f'umap_{g}_log10tpm.new.pdf')
    plt.close()




## color by zscore
# 0 as center
vcenter = 0
# vmin, vmax = fcstat.iloc[:,:171].min().min(), fcstat.iloc[:,:171].max().max()
vmin, vmax = (-2,2)
normalize = mcolors.TwoSlopeNorm(vcenter=vcenter, vmin=vmin, vmax=vmax)

groups = tpm_zscore.columns

for g in groups:
    plt.figure(figsize=(5.5,4))

    ax = sns.scatterplot(x='UMAP1', y='UMAP2', data=embedding,
                    hue=tpm_zscore[g].to_list(),
                    linewidth=0, s=3, alpha=.7, edgecolor=None,
                    palette="viridis",
                    hue_norm=normalize)
    
    
    # norm = plt.Normalize(fcstat.iloc[:,:171].min().min(), fcstat.iloc[:,:171].max().max())
    sm = plt.cm.ScalarMappable(cmap="viridis", norm=normalize)
    sm.set_array([])
    
    # Remove the legend and add a colorbar
    ax.get_legend().remove()
    ax.figure.colorbar(sm)
    
    # plt.show()
    plt.tight_layout()
    plt.savefig(f'umap_{g}_zscore.new.pdf')
    plt.close()



## convert ensemble id to gene symbol
#### this is way too slow
#converter = EnsemblConverter()
#symbolids = converter.convert_ids(embedding.index)
####


## color with K37ac levels of enahcners of these genes
ac_hancer = pd.read_csv('../../mark_profile/signal_on_enhancer/msnormedct_on_enhancer.txt', header=0, index_col=None, sep='\t')

# group by gene
ave_bygene = ac_hancer.iloc[:,4:].groupby('connected_gene').mean()

# mean of replicates
ave_bygene_mean = ave_bygene.groupby(ave_bygene.columns.str.replace('_Rep1|_Rep2','', regex=True), axis=1).mean()


# keep only genes in the embedding
plot_dat = embedding.merge(ave_bygene_mean, left_index=True, right_index=True, how='inner')



## unfinished
sns.scatterplot(x='UMAP1', y='UMAP2', data=plot_dat,
                    hue=np.log2(plot_dat['NSD2i_D5_H3K27Ac']+0.1).to_list(),
                    linewidth=0, s=3, alpha=.7, edgecolor=None,
                    palette="viridis")

### did not show patterns'





# ### if embed base on N/C
# ## does not work well
# tpm_normed = np.log2(np.divide(*(tpm.iloc[:,np.where(tpm.columns.str.contains('_N'))[0]]+1).\
#                    align((tpm.iloc[:,np.where(tpm.columns.str.contains('_C'))[0]]+1), axis=0)))

# embedding_norm = pd.DataFrame(umap.fit_transform(tpm_normed), columns = ['UMAP1','UMAP2'])
# embedding_norm.index = tpm.index

# umap = UMAP(n_neighbors=10, min_dist=0)

# embedding_norm = pd.DataFrame(umap.fit_transform(tpm_normed), columns = ['UMAP1','UMAP2'])
# embedding_norm.index = tpm.index


# vcenter = 0
# # vmin, vmax = fcstat.iloc[:,:171].min().min(), fcstat.iloc[:,:171].max().max()
# vmin, vmax = (-2,2)
# normalize = mcolors.TwoSlopeNorm(vcenter=vcenter, vmin=vmin, vmax=vmax)

# groups = tpm_zscore.columns

# for g in groups:
#     plt.figure(figsize=(5.5,4))

#     ax = sns.scatterplot(x='UMAP1', y='UMAP2', data=embedding_norm,
#                     hue=tpm_zscore[g].to_list(),
#                     linewidth=0, s=3, alpha=.7, edgecolor=None,
#                     palette="viridis",
#                     hue_norm=normalize)
    
    
#     # norm = plt.Normalize(fcstat.iloc[:,:171].min().min(), fcstat.iloc[:,:171].max().max())
#     sm = plt.cm.ScalarMappable(cmap="viridis", norm=normalize)
#     sm.set_array([])
    
#     # Remove the legend and add a colorbar
#     ax.get_legend().remove()
#     ax.figure.colorbar(sm)
    
#     # plt.show()
#     plt.tight_layout()