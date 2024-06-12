import os,glob
import pandas as pd
import numpy as np
import seaborn as sns
import matplotlib
from matplotlib import pyplot as plt
import matplotlib.font_manager

plt.style.use('seaborn-poster')
matplotlib.rcParams['pdf.fonttype'] = 42
matplotlib.rcParams['ps.fonttype'] = 42

## load ensemble id to symbol matching table
ensembl_syb = pd.read_csv('../rnaseq/umap/human_ensembl_syb.tsv', header=0, index_col=None, sep='\t')
ensembl_syb = dict(zip(ensembl_syb['ensembl_gene_id'], ensembl_syb['hgnc_symbol']))

## load tpm table
tpm = pd.read_csv('../rnaseq/tpm/tximport-tpm.csv', header=0, index_col=0)

## change index to gene symbol
tpm = tpm.set_index(tpm.index.map(ensembl_syb))
tpm = tpm[tpm.index.notnull()]

# merge transcripts, only ~100 genes with more than 1 transcript
tpm = tpm.groupby(tpm.index).mean()

## log2tpm
tpm_normed = np.log2(np.divide(*(tpm.iloc[:,np.where(tpm.columns.str.contains('_N'))[0]]+1).\
                    align((tpm.iloc[:,np.where(tpm.columns.str.contains('_C'))[0]]+1), axis=0)))

tpm_normed.columns = tpm_normed.columns.str.replace('N', 'Rep').str.replace('Day', 'D') + '_log2tpm'


## load enhancer annotaion to get the genes that eventually will be using
ct_enhancer_normed = pd.read_csv('../mark_profile/signal_on_enhancer/ct_enhancer_normed.txt', header=0, index_col=0, sep='\t')
ct_enhancer_normed.columns = ct_enhancer_normed.columns.str.replace('NSD2i_', '')


## merged gene set
# gene_uni = set(ct_enhancer_normed.index) & set(tpm_normed.index)
gene_uni = set(tpm_normed.index)


# kras down
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


kras_dn = [x for x in kras_dn if x in gene_uni]


### selecting random genes from conrol d1
krasdn_d1ctrl = tpm.loc[kras_dn,tpm.columns.str.contains('Day1_C')]
all_d1ctrl = tpm.loc[gene_uni,tpm.columns.str.contains('Day1_C')]

def euclidean_distance(row1, row2):
    return np.sqrt(np.sum((row1 - row2) ** 2))

# Calculate distances
ctrl_sets = []

for i in range(krasdn_d1ctrl.shape[0]):
    distances = all_d1ctrl.apply(lambda row: euclidean_distance(row.values, krasdn_d1ctrl.iloc[[i]].values), axis=1)
    sorted_indices = np.argsort(distances)

    # keep the 5 rows with the smallest distances
    closest_indices = [idx for idx in sorted_indices if all_d1ctrl.index[idx] in gene_uni][1:6]
    ctrl_sets.append(closest_indices)



## convert index to genes
ctrl_set_genes = [all_d1ctrl.index[idx].to_list() for idx in ctrl_sets]

ctrl_genes = pd.DataFrame(ctrl_set_genes)
# ctrl_genes = ctrl_genes.drop_duplicates()

ctrl_genes.columns = ['group1','group2','group3','group4','group5']
ctrl_genes['krasdn'] = kras_dn



## box plot of kras dn vs 5 control sets

kras_dnset = tpm_normed.loc[kras_dn].stack().to_frame().assign(geneset='krasdn')
ctrl_set1 = tpm_normed.loc[ctrl_genes.iloc[:,0]].stack().to_frame().assign(geneset='ctrl1')
ctrl_set2 = tpm_normed.loc[ctrl_genes.iloc[:,1]].stack().to_frame().assign(geneset='ctrl2')
ctrl_set3 = tpm_normed.loc[ctrl_genes.iloc[:,2]].stack().to_frame().assign(geneset='ctrl3')
ctrl_set4 = tpm_normed.loc[ctrl_genes.iloc[:,3]].stack().to_frame().assign(geneset='ctrl4')
ctrl_set5 = tpm_normed.loc[ctrl_genes.iloc[:,4]].stack().to_frame().assign(geneset='ctrl5')


## concat all sets
all_sets = pd.concat([kras_dnset,ctrl_set1,ctrl_set2,ctrl_set3,ctrl_set4,ctrl_set5]).reset_index()
all_sets['day'] = all_sets['level_1'].str.split('_').str[0]
all_sets['rep'] = all_sets['level_1'].str.split('_').str[1]



## plot the box
g = sns.catplot(data=all_sets, kind="violin",
            palette=["#4876FF", "#CDC9C9", "#CDC9C9", "#CDC9C9","#CDC9C9","#CDC9C9"],
            x="geneset", y=0, col="day", aspect=1, 
            linewidth=1.2)
g.set_axis_labels("", "log2(NSD2i/Vehicle) TPM")
plt.savefig('ctrl_from_D1.pdf', bbox_inches='tight')
plt.close()



####### selecting random genes log2 normed tpm considering all time points ######
krasdn_normed = tpm_normed.loc[kras_dn,:]
other_normed = tpm_normed.loc[gene_uni,:]

# Calculate distances
ctrl_sets_normed = []

for i in range(krasdn_normed.shape[0]):
    distances = other_normed.apply(lambda row: euclidean_distance(row.values, krasdn_normed.iloc[[i]].values), axis=1)
    sorted_indices = np.argsort(distances)

    # keep the 5 rows with the smallest distances
    closest_indices = [idx for idx in sorted_indices if all_d1ctrl.index[idx] in gene_uni][1:6]
    ctrl_sets_normed.append(closest_indices)

## convert index to genes
ctrl_sets_normed = [other_normed.index[idx].to_list() for idx in ctrl_sets_normed]


ctrl_genes_normed = pd.DataFrame(ctrl_sets_normed)
# ctrl_genes_normed = ctrl_genes_normed.drop_duplicates()

ctrl_genes_normed.columns = ['group1','group2','group3','group4','group5']
ctrl_genes_normed['krasdn'] = kras_dn


## box plot of kras dn vs 5 control sets
ctrl_set1_normed = tpm_normed.loc[ctrl_genes_normed.iloc[:,0]].stack().to_frame().assign(geneset='ctrl1')
ctrl_set2_normed = tpm_normed.loc[ctrl_genes_normed.iloc[:,1]].stack().to_frame().assign(geneset='ctrl2')
ctrl_set3_normed = tpm_normed.loc[ctrl_genes_normed.iloc[:,2]].stack().to_frame().assign(geneset='ctrl3')
ctrl_set4_normed = tpm_normed.loc[ctrl_genes_normed.iloc[:,3]].stack().to_frame().assign(geneset='ctrl4')
ctrl_set5_normed = tpm_normed.loc[ctrl_genes_normed.iloc[:,4]].stack().to_frame().assign(geneset='ctrl5')


## concat all sets
all_sets_normed = pd.concat([kras_dnset,ctrl_set1_normed,ctrl_set2_normed,ctrl_set3_normed,ctrl_set4_normed,ctrl_set5_normed]).reset_index()
all_sets_normed['day'] = all_sets_normed['level_1'].str.split('_').str[0]
all_sets_normed['rep'] = all_sets_normed['level_1'].str.split('_').str[1]



## plot the box
g = sns.catplot(data=all_sets_normed, kind="violin",
            palette=["#4876FF", "#CDC9C9", "#CDC9C9", "#CDC9C9","#CDC9C9","#CDC9C9"],
            x="geneset", y=0, col="day", aspect=1, 
            linewidth=1.2)
g.set_axis_labels("", "log2(NSD2i/Vehicle) TPM")
plt.savefig('ctrl_from_allnormedtpm.pdf', bbox_inches='tight')
plt.close()


# export
ctrl_genes.to_csv('ctrlsets_from_D1.csv', index=False)
ctrl_genes_normed.to_csv('ctrlsets_from_allnormedtpm.csv', index=False)
