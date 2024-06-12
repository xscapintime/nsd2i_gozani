import os,glob
import pandas as pd
import numpy as np
import seaborn as sns
import matplotlib
from matplotlib import pyplot as plt
import matplotlib.font_manager
from itertools import product
from statannotations.Annotator import Annotator


plt.style.use('seaborn-poster')
matplotlib.rcParams['pdf.fonttype'] = 42
matplotlib.rcParams['ps.fonttype'] = 42


## load log2fc cut&run
ct_enhancer_normed = pd.read_csv('ct_enhancer_normed.txt', header=0, index_col=0, sep='\t')
ct_enhancer_normed.columns = ct_enhancer_normed.columns.str.replace('NSD2i_', '')


# stats of enhancer count
enhancer_count = ct_enhancer_normed.index.value_counts().to_dict()

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
kras_ctrl_sets = pd.read_csv('../../kras_vsctrl/ctrlsets_from_D1.csv')

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


## stats data

def count_signs(column):
    positive = (column > 0).sum()
    negative = (column < 0).sum()
    nochange  = (column == 0).sum()
    return pd.Series({'gain': positive, 'loss': negative, 'nochange': nochange})

stas = []
for n in ['kras_df', 'ctrl1', 'ctrl2', 'ctrl3', 'ctrl4', 'ctrl5']:
    dd = eval(n)
    countstat = dd.loc[:,dd.columns.str.contains('H3K27Ac')].apply(count_signs).stack().reset_index()
    # countstat['rep'] = countstat['level_1'].str.split('_').str[1]
    countstat['day'] = countstat['level_1'].str.split('_').str[0]
    countstat = countstat.groupby(['level_0', 'day']).mean().reset_index()
    countstat['pergene'] = countstat[0] / len(kras_dn)
    countstat['set'] = n
    stas.append(countstat)

stas_df = pd.concat(stas)


# plot enhancer per gene
plt.figure(figsize=(6, 5))
g = sns.barplot(data=stas_df.groupby(['day','set']).sum().loc[['D1']].reset_index(),
            x='set', y='pergene',
            palette=["#4876FF", "#CDC9C9", "#CDC9C9", "#CDC9C9","#CDC9C9","#CDC9C9"],
            order=['kras_df', 'ctrl1', 'ctrl2', 'ctrl3', 'ctrl4', 'ctrl5'])
g.set(xlabel="Genesets", ylabel="Ehancer per gene") 
plt.savefig('enhancer_per_gene.pdf', bbox_inches='tight')


# plot by lose or gain k27ac
g = sns.catplot(data=stas_df, x='set', y='pergene', hue='level_0', kind='bar', palette='Set2',col='day')
g.set_axis_labels("Genesets", "Ehancer per gene")
plt.savefig('enhancer_k27actype_per_gene.pdf', bbox_inches='tight')



### enhancer count violin dot plot + t-test
count_df = pd.DataFrame()

for column in kras_ctrl_sets.columns:
    count_df[column] = kras_ctrl_sets[column].map(enhancer_count).fillna(0)

count_df.index = kras_ctrl_sets.krasdn
count_df = count_df[['krasdn','group1', 'group2', 'group3', 'group4', 'group5']]

plot_df = count_df.stack().reset_index()


plt.figure(figsize=(7, 6))
g = sns.violinplot(data=plot_df, x='level_1', y=0,
               palette = ["#4876FF", "#CDC9C9", "#CDC9C9", "#CDC9C9","#CDC9C9","#CDC9C9"],
               saturation=0.7, linewidth=1, inner='box')

g = sns.stripplot(data=plot_df, x="level_1", y=0, jitter=True,
              color='black', alpha=0.3, dodge=True, size=2)

plt.ylabel('Enhancer number per gene')
plt.xlabel('')

# stats
comb = list(product([count_df.columns[0]],count_df.columns[1:]))

annotator = Annotator(g, comb, data=plot_df, x='level_1', y=0)
annotator.configure(test="t-test_paired", text_format="full",loc='inside')
annotator.apply_and_annotate()

plt.savefig('enhancer_count_violin.pdf', bbox_inches='tight')




## if remove SGK1
count_df_new = count_df.drop('SGK1', inplace=False)

plot_df_new = count_df_new.stack().reset_index()

plt.figure(figsize=(7, 6))
g = sns.violinplot(data=plot_df_new, x='level_1', y=0,
               palette = ["#4876FF", "#CDC9C9", "#CDC9C9", "#CDC9C9","#CDC9C9","#CDC9C9"],
               saturation=0.7, linewidth=1, inner='box')

g = sns.stripplot(data=plot_df_new, x="level_1", y=0, jitter=True,
              color='black', alpha=0.3, dodge=True, size=2)

plt.ylabel('Enhancer number per gene')
plt.xlabel('')

# stats
comb = list(product([count_df_new.columns[0]],count_df_new.columns[1:]))

annotator = Annotator(g, comb, data=plot_df_new, x='level_1', y=0)
annotator.configure(test="t-test_paired", text_format="full",loc='inside')
annotator.apply_and_annotate()

plt.savefig('enhancer_count_violin_drop.pdf', bbox_inches='tight')


