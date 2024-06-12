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

# stats of enhancer count
enhancer_count = ct_enhancer_normed.index.value_counts().to_dict()

## load ensemble id to symbol matching table
ensembl_syb = pd.read_csv('../../rnaseq/umap/human_ensembl_syb.tsv', header=0, index_col=None, sep='\t')
ensembl_syb = dict(zip(ensembl_syb['ensembl_gene_id'], ensembl_syb['hgnc_symbol']))

## load tpm table
tpm = pd.read_csv('../../rnaseq/tpm/tximport-tpm.csv', header=0, index_col=0)

## change index to gene symbol
tpm = tpm.set_index(tpm.index.map(ensembl_syb))
tpm = tpm[tpm.index.notnull()]

# merge transcripts, only ~100 genes with more than 1 transcript
tpm = tpm.groupby(tpm.index).mean()

## log2tpm
tpm_normed = np.log2(np.divide(*(tpm.iloc[:,np.where(tpm.columns.str.contains('_N'))[0]]+1).\
                    align((tpm.iloc[:,np.where(tpm.columns.str.contains('_C'))[0]]+1), axis=0)))

tpm_normed.columns = tpm_normed.columns.str.replace('N', 'Rep').str.replace('Day', 'D') + '_log2tpm'


## load control gene sets
kras_ctrl_sets = pd.read_csv('../../kras_vsctrl/ctrlsets_from_D1.csv')
kras_ctrl_sets = kras_ctrl_sets[['krasdn','group1','group2','group3','group4','group5']]
kras_ctrl_sets.index = kras_ctrl_sets['krasdn']

## order of base expression (day1 control)
base_expr = tpm.loc[kras_ctrl_sets['krasdn'],tpm.columns.str.contains('Day1_C')].\
    apply(lambda x: np.mean(x), axis=1).sort_values(ascending=False).index

kras_ctrl_sets = kras_ctrl_sets.loc[base_expr]


## map enhancer count
kras_ctrl_sets_enhancern= pd.DataFrame()

for column in kras_ctrl_sets.columns:
    kras_ctrl_sets_enhancern[column] = kras_ctrl_sets[column].map(enhancer_count)

g = sns.relplot(
    data=kras_ctrl_sets_enhancern.stack().reset_index(), x="gene", 
    y=kras_ctrl_sets_enhancern.stack().reset_index()[0],
    col="level_1",kind="line",linewidth=2, hue="level_1",
    col_wrap=3, height=3, aspect=1.5,
    palette=["#4876FF", "#CDC9C9", "#CDC9C9", "#CDC9C9","#CDC9C9","#CDC9C9"],
)
g.set(xticks=[])

plt.savefig('kras_ctrl_sets_enhancern_baseexprorder.pdf', bbox_inches='tight')

## map log2tpm
krasdn_tpmnormed = tpm_normed.loc[kras_ctrl_sets.index].stack().reset_index()
krasdn_tpmnormed.columns = ['gene','sample','log2tpm']
krasdn_tpmnormed['set'] = 'krasdn'
# krasdn_tpmnormed['krasrefgene'] = krasdn_tpmnormed['gene']

for n in ['group1','group2','group3','group4','group5']:
    ctrl_df = tpm_normed.loc[kras_ctrl_sets[n]].stack().reset_index()
    ctrl_df.columns = ['gene','sample', 'log2tpm']
    ctrl_df['set'] = n
    # ctrl_df['krasrefgene'] =  krasdn_tpmnormed['gene']

    krasdn_tpmnormed = pd.concat([krasdn_tpmnormed, ctrl_df], axis=0)

krasdn_tpmnormed['refkrasgene'] = pd.concat([kras_ctrl_sets['krasdn'].repeat(12)]*6, ignore_index=True)
krasdn_tpmnormed['day'] = krasdn_tpmnormed['sample'].str.split('_').str[0]
krasdn_tpmnormed['rep'] = krasdn_tpmnormed['sample'].str.split('_').str[1]
krasdn_tpmnormed['enhancer_n'] = krasdn_tpmnormed['gene'].map(enhancer_count)



## 
g = sns.relplot(
    data=krasdn_tpmnormed, x="refkrasgene", 
    y='log2tpm',col="set", row='day',
    kind="line",linewidth=2, hue="set",
    height=3, aspect=1.5,
    palette=["#4876FF", "#CDC9C9", "#CDC9C9", "#CDC9C9","#CDC9C9","#CDC9C9"],
)
g.set(xticks=[])
plt.savefig('kras_ctrl_sets_log2tpm_baseexprorder.pdf', bbox_inches='tight')




g = sns.relplot(
    data=krasdn_tpmnormed, y="enhancer_n", 
    x='log2tpm',col="set", row='day',
    kind="scatter", hue="set",s=10,
    height=3, aspect=1.5, edgecolor=None,
    palette=["#4876FF", "#CDC9C9", "#CDC9C9", "#CDC9C9","#CDC9C9","#CDC9C9"],
)
g.set(xticks=[])
plt.savefig('kras_ctrl_sets_log2tpm_enhancern_baseexprorder.pdf', bbox_inches='tight')




## heatmap
avrage_d1 = krasdn_tpmnormed[krasdn_tpmnormed['day'] =='D1'].groupby(['gene','day']).mean().reset_index()

d1df = krasdn_tpmnormed[krasdn_tpmnormed['day'] =='D1'].groupby


sns.heatmap(data=krasdn_tpmnormed[krasdn_tpmnormed['day'] =='D1'].pivot_table(index='refkrasgene', columns='set', values='log2tpm'),
            cmap='coolwarm', center=0)