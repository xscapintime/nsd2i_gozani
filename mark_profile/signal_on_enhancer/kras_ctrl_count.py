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


## load log2fc cut&run
ct_enhancer_normed = pd.read_csv('ct_enhancer_normed.txt', header=0, index_col=0, sep='\t')
ct_enhancer_normed.columns = ct_enhancer_normed.columns.str.replace('NSD2i_', '')


## kras down
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

## ctrl sets
kras_ctrl_sets = pd.read_csv('../../kras_vsctrl/ctrlsets_from_allnormedtpm.csv')



## count enhancer number
kras_df = ct_enhancer_normed[ct_enhancer_normed.index.isin(kras_dn)]

ctrl_sets = {}

# Loop to create each DataFrame
for i in range(1, 6):
    key = f'ctrl{i}'
    
    # Use .isin to filter and stack, reset index for the desired operation
    df =  ct_enhancer_normed[ct_enhancer_normed.index.isin(kras_ctrl_sets.iloc[:, i-1])]

    # Assign the DataFrame to the dictionary
    ctrl_sets[key] = df


ctrl1 = ctrl_sets['ctrl1']
ctrl2 = ctrl_sets['ctrl2']
ctrl3 = ctrl_sets['ctrl3']
ctrl4 = ctrl_sets['ctrl4']
ctrl5 = ctrl_sets['ctrl5']


## stats df
stas_df = pd.DataFrame(columns=['Gene Number', 'Enhancer Number', 'Average Enhancer Number'],
                       index=['kras_df', 'ctrl1', 'ctrl2', 'ctrl3', 'ctrl4', 'ctrl5'])


for n in ['kras_df', 'ctrl1', 'ctrl2', 'ctrl3', 'ctrl4', 'ctrl5']:
    df = eval(n)
    gene_number = len(df.index.unique())
    enhancer_number = df.shape[0]
    avg_enhancer_number = enhancer_number / gene_number

    stas_df.loc[n] = [gene_number, enhancer_number, avg_enhancer_number]


sns.barplot(data=stas_df, x=stas_df.index, y='Average Enhancer Number')



##########------------------------------------##########
# split by k27ac
new_dat = kras_ctrl_ct_tmp[['connected_gene', 'geneset','mark_x','value_x','day']].drop_duplicates()
new_ac = new_dat[new_dat['mark_x'] == 'H3K27Ac']
new_ac['group'] = np.where(new_ac['value_x'] == 0, 'no change', np.where(new_ac['value_x'] > 0, 'gain', 'loss'))


g = sns.catplot(data=new_ac, x='geneset', kind='count',
            palette=["#4876FF", "#CDC9C9", "#CDC9C9", "#CDC9C9","#CDC9C9","#CDC9C9"],
            col='group', aspect=1.4)


sns.countplot(data=kras_ctrl_ct_tmp[['connected_gene', 'geneset']].drop_duplicates(), x='geneset')