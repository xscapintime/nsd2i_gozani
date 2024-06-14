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


# merge transcripts, only ~100 genes with more than 1 transcript
tpm = tpm.groupby(tpm.index).mean()

## log2tpm
tpm_normed = np.log2(np.divide(*(tpm.iloc[:,np.where(tpm.columns.str.contains('_N'))[0]]+1).\
                    align((tpm.iloc[:,np.where(tpm.columns.str.contains('_C'))[0]]+1), axis=0)))

tpm_normed.columns = tpm_normed.columns.str.replace('N', 'Rep').str.replace('Day', 'D') + '_log2tpm'


## load log2fc cut&run
ct_enhancer_normed = pd.read_csv('../../mark_profile/signal_on_enhancer/ct_enhancer_normed.txt', header=0, index_col=0, sep='\t')
ct_enhancer_normed.columns = ct_enhancer_normed.columns.str.replace('NSD2i_', '')

# merge enhancers??
# ct_enhancer_normed = ct_enhancer_normed.groupby(ct_enhancer_normed.index).mean()

## merge tpm and ct_enhancer
tpm_enhancersig = tpm_normed.merge(ct_enhancer_normed, left_index=True, right_index=True, how='inner')


## merge replicates for coloring
rep_merged = tpm_enhancersig.groupby(tpm_enhancersig.columns.str.replace('_Rep1|_Rep2|_Rep3|_Rep4','', regex=True), axis=1).mean()

# distribution

rep_merged_long = rep_merged.stack().reset_index()
rep_merged_long['day'] = rep_merged_long['level_1'].str.split('_').str[0]
rep_merged_long['mark'] = rep_merged_long['level_1'].str.split('_').str[1]

g = sns.FacetGrid(rep_merged_long, col="day", row="mark", margin_titles=True)
g.map(sns.histplot, 0)
plt.savefig('histplot_alllog2.pdf')



#### UMAP ####
# seed
np.random.seed(1024)

# UMAP demension redunction on tpm of all samples
umap = UMAP(n_neighbors=25, min_dist=0.15, metric='manhattan')


# embedding
embedding = pd.DataFrame(umap.fit_transform(tpm_enhancersig), columns = ['UMAP1','UMAP2'])
embedding.index = tpm_enhancersig.index


########## viz #######
# 0 as center
vcenter = 0
# vmin, vmax = fcstat.iloc[:,:171].min().min(), fcstat.iloc[:,:171].max().max()
vmin, vmax = (-2,2)
normalize = mcolors.TwoSlopeNorm(vcenter=vcenter, vmin=vmin, vmax=vmax)

for g in rep_merged.columns:
    plt.figure(figsize=(7,5))

    ax = sns.scatterplot(x='UMAP1', y='UMAP2', data=embedding,
                    hue=rep_merged.loc[:,g].to_list(),
                    linewidth=0, s=1, alpha=.5, edgecolor=None,
                    palette="viridis",
                    hue_norm=normalize)
    
    
    # norm = plt.Normalize(fcstat.iloc[:,:171].min().min(), fcstat.iloc[:,:171].max().max())
    sm = plt.cm.ScalarMappable(cmap="viridis", norm=normalize)
    sm.set_array([])
    
    # Remove the legend and add a colorbar
    ax.get_legend().remove()
    ax.figure.colorbar(sm)
    
    ax.set_title(g)

    # plt.show()
    plt.tight_layout()
    plt.savefig(f'umap_{g}.png')
    plt.close()




## only color kras genes
kras_dn = ['ABCB11','ABCG4','ACTC1','ADRA2C','AKR1B10','ALOX12B','AMBN','ARHGDIG','ARPP21','ASB7','ATP4A','ATP6V1B1','BARD1',
            'BMPR1B','BRDT','BTG2','C5','CACNA1F','CACNG1','CALCB','CALML5','CAMK1D','CAPN9','CCDC106','CCNA1','CCR8','CD207',
            'CD40LG','CD80','CDH16','CDKAL1','CELSR2','CHRNG','CHST2','CKM','CLDN16','CLDN8','CLPS','CLSTN3','CNTFR','COL2A1',
            'COPZ2','COQ8A','CPA2','CPB1','CPEB3','CYP11B2','CYP39A1','DCC','DLK2','DTNB','EDAR','EDN1','EDN2','EFHD1','EGF',
            'ENTPD7','EPHA5','FGF16','FGF22','FGFR3','FGGY','FSHB','GAMT','GDNF','GP1BA','GP2','GPR19','GPR3','GPRC5C','GRID2',
            'GTF3C5','HNF1A','HSD11B2','HTR1B','HTR1D','IDUA','IFI44L','IFNG','IGFBP2','IL12B','IL5','INSL5','IRS4','ITGB1BP2',
            'ITIH3','KCND1','KCNE2','KCNMB1','KCNN1','KCNQ2','KLHDC8A','KLK7','KLK8','KMT2D','KRT1','KRT13','KRT15','KRT4','KRT5',
            'LFNG','LGALS7','LYPD3','MACROH2A2','MAGIX','MAST3','MEFV','MFSD6','MSH5','MTHFR','MX1','MYH7','MYO15A','MYOT','NGB',
            'NOS1','NPHS1','NPY4R','NR4A2','NR6A1','NRIP2','NTF3','NUDT11','OXT','P2RX6','P2RY4','PAX3','PAX4','PCDHB1','PDCD1',
            'PDE6B','PDK2','PKP1','PLAG1','PNMT','PRKN','PRODH','PROP1','PTGFR','PTPRJ','RGS11','RIBC2','RSAD2','RYR1','RYR2',
            'SCGB1A1','SCN10A','SELENOP','SERPINA10','SERPINB2','SGK1','SHOX2','SIDT1','SKIL','SLC12A3','SLC16A7','SLC25A23',
            'SLC29A3','SLC30A3','SLC38A3','SLC5A5','SLC6A14','SLC6A3','SMPX','SNCB','SNN','SOX10','SPHK2','SPRR3','SPTBN2',
            'SSTR4','STAG3','SYNPO','TAS2R4','TCF7L1','TCL1A','TENM2','TENT5C','TEX15','TFAP2B','TFCP2L1','TFF2','TG','TGFB2',
            'TGM1','THNSL2','THRB','TLX1','TNNI3','TSHB','UGT2B17','UPK3B','VPREB1','VPS50','WNT16','YBX2','YPEL1','ZBTB16',
            'ZC2HC1C','ZNF112']


for g in rep_merged.columns:
    plt.figure(figsize=(7,5))

    ax = sns.scatterplot(x='UMAP1', y='UMAP2',
                    data=embedding[~embedding.index.isin(kras_dn)],
                    linewidth=0, s=1, alpha=.3, edgecolor=None,
                    color="#8B8989")

    ax = sns.scatterplot(x='UMAP1', y='UMAP2',
                    data=embedding[embedding.index.isin(kras_dn)],
                        hue=rep_merged.loc[rep_merged.index.isin(kras_dn),g].to_list(),
                        linewidth=0, s=3, alpha=.7, edgecolor=None,
                        palette="viridis",
                        hue_norm=normalize)

    # norm = plt.Normalize(fcstat.iloc[:,:171].min().min(), fcstat.iloc[:,:171].max().max())
    sm = plt.cm.ScalarMappable(cmap="viridis", norm=normalize)
    sm.set_array([])
    
    # Remove the legend and add a colorbar
    ax.get_legend().remove()
    ax.figure.colorbar(sm)
    
    ax.set_title(g)

    # plt.show()
    plt.tight_layout()
    plt.savefig(f'umap_{g}_kras.png')
    plt.close()




# ## color by zscore
# # 0 as center
# vcenter = 0
# # vmin, vmax = fcstat.iloc[:,:171].min().min(), fcstat.iloc[:,:171].max().max()
# vmin, vmax = (-2,2)
# normalize = mcolors.TwoSlopeNorm(vcenter=vcenter, vmin=vmin, vmax=vmax)

# groups = tpm_zscore.columns

# for g in groups:
#     plt.figure(figsize=(5.5,4))

#     ax = sns.scatterplot(x='UMAP1', y='UMAP2', data=embedding,
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
#     plt.savefig(f'umap_{g}_zscore.new.pdf')
#     plt.close()



