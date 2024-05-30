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


## load enhancer count
ac_hancer = pd.read_csv('../../mark_profile/signal_on_enhancer/msnormedct_on_enhancer.txt', header=0, index_col=None, sep='\t')
ac_hancer.groupby('connected_gene')['genehancer_id'].count().sort_values(ascending=False).head(10)


# NSD2i/CK
nsd2i =  ac_hancer[sorted(ac_hancer.columns[ac_hancer.columns.str.contains('NSD2i')])]
ck =  ac_hancer[sorted(ac_hancer.columns[ac_hancer.columns.str.contains('Vehicle')])]
ac_normed = np.log2((nsd2i+1)/(ck+1).values)
ac_normed = ac_normed[np.sum(ac_normed,axis=1) > 0]

# add gene symbol and enhancer id
ac_normed = pd.merge(ac_hancer[['genehancer_id','connected_gene']], ac_normed, left_index=True, right_index=True)


# find genes loses k36me2 at enhancers
consisrep_idx = ac_normed.loc[:,ac_normed.columns.str.contains('H3K36me2')].\
    apply(lambda x: x[0] * x[1], axis=1).loc[lambda x: x > 0].index

least_k36 = ac_normed.loc[consisrep_idx, ac_normed.columns.str.contains('H3K36me2')].\
        apply(lambda x: np.mean(x), axis=1).sort_values().index

# export, ordered by increasing k36me2
ac_normed.loc[least_k36].to_csv('ct_enhancer_normed_leastk36me2.txt', header=True, index=False, sep='\t')


# filter RNA gene, lncRNA, and pseudogene
ac_normed_norna = ac_normed.loc[least_k36].loc[~ac_normed['connected_gene'].str.contains('lnc|HSALNG|piR|LOC|ENSG|RF00|LINC0|CM0|MIR|RNU|NONHSA|AS1|LINC')]

ac_normed_norna.to_csv('ct_enhancer_normed_leastk36me2_norna.txt', header=True, index=False, sep='\t')


## top genes losing k36me2 at enhancers
topgenes = ac_normed_norna.loc[ac_normed_norna.loc[:, ac_normed_norna.columns.str.contains('H3K36me2')].\
            apply(lambda x: np.mean(x), axis=1).sort_values(ascending=True).index].head(100)['connected_gene']

# count enhancer number
topg_n_enhancer = ac_normed_norna.loc[ac_normed_norna['connected_gene'].isin(topgenes)].groupby('connected_gene')['genehancer_id'].nunique().sort_values(ascending=False).to_dict()


# plot data

plot_dat = ac_normed_norna.loc[ac_normed_norna['connected_gene'].isin(list(topg_n_enhancer.keys())[0:20])].melt(id_vars=["genehancer_id", "connected_gene"], var_name="original_col", value_name="value")
plot_dat['day'] = plot_dat['original_col'].str.split('_').str[1]
plot_dat['rep'] = plot_dat['original_col'].str.split('_').str[2]
plot_dat['mark'] = plot_dat['original_col'].str.split('_').str[3]


plot_dat['connected_gene'] = pd.Categorical(plot_dat['connected_gene'],\
                        categories=plot_dat.loc[plot_dat.mark == 'H3K36me2'].groupby('connected_gene').\
                            mean('value').sort_values(by='value').index)


## plot

g = sns.catplot(data=plot_dat, kind="box",palette="pastel",
            x="connected_gene", y="value", row="mark", col='day', aspect=1.4, 
            showfliers=False, linewidth=1.2)

#https://stackoverflow.com/questions/67309730/how-to-overlay-a-scatterplot-on-top-of-boxplot-with-sns-catplot
g.map_dataframe(sns.stripplot, x="connected_gene", y="value", 
                palette=["#404040"],
                alpha=0.6,dodge=True,
                size=3)
# g.set_xticklabels(rotation=90,rotation_mode='anchor',ha='left', va='top')
for ax in g.axes.flat:
    ax.axhline(0, color='black', linestyle='--', linewidth=0.5)
    
    for label in ax.get_xticklabels():
        label.set_rotation(90)
        label.set_horizontalalignment('center')

g.axes[1, 0].set_ylabel('log2(NSD2i/Vehicle)') 

plt.savefig('enhancer_normed_leastk36me2_topgenes.pdf', bbox_inches='tight')
plt.close()


## add gene expression
## load ensemble id to symbol matching table
ensembl_syb = pd.read_csv('../../rnaseq/umap/human_ensembl_syb.tsv', header=0, index_col=None, sep='\t')
ensembl_syb = dict(zip(ensembl_syb['ensembl_gene_id'], ensembl_syb['hgnc_symbol']))

## load tpm table
tpm = pd.read_csv('../../rnaseq/tpm/tximport-tpm.csv', header=0, index_col=0)

## change index to gene symbol
tpm = tpm.set_index(tpm.index.map(ensembl_syb))
tpm = tpm[tpm.index.notnull()]

## log2tpm
tpm_normed = np.log2(np.divide(*(tpm.iloc[:,np.where(tpm.columns.str.contains('_N'))[0]]+1).\
                    align((tpm.iloc[:,np.where(tpm.columns.str.contains('_C'))[0]]+1), axis=0)))

tpm_normed.columns = tpm_normed.columns.str.replace('N', 'Rep').str.replace('Day', 'D') + '_log2tpm'
tpm_long = tpm_normed.stack().reset_index()
tpm_long['day'] = tpm_long['level_1'].str.split('_').str[0]
tpm_long['rep'] = tpm_long['level_1'].str.split('_').str[1]

tpm_long.columns = ['connected_gene', 'original_col', 'value', 'day', 'rep']


tpm_long.connected_gene in plot_dat.connected_gene

plot_dat_expr = pd.merge(plot_dat, tpm_long.loc[tpm_long['connected_gene'].isin(list(topg_n_enhancer.keys())[0:20])],
         on=['connected_gene', 'day', 'rep'], how='right')


## plot
g = sns.catplot(data=plot_dat_expr, kind="box",palette="pastel",
            x="connected_gene", y="value_x", row="mark", col='day', aspect=1.4, 
            showfliers=False, linewidth=1.2)

#https://stackoverflow.com/questions/67309730/how-to-overlay-a-scatterplot-on-top-of-boxplot-with-sns-catplot
g.map_dataframe(sns.stripplot, x="connected_gene", y="value_x", 
                palette=["#404040"],
                alpha=0.6,dodge=True,
                size=3)

g.map_dataframe(sns.stripplot, x="connected_gene", y="value_y", 
                palette=["#CD5555"],
                alpha=0.5,dodge=True,
                s=5, marker="D",
                jitter=False)

for ax in g.axes.flat:
    ax.axhline(0, color='black', linestyle='--', linewidth=0.5)
    
    for label in ax.get_xticklabels():
        label.set_rotation(90)
        label.set_horizontalalignment('center')

g.axes[1, 0].set_ylabel('log2(NSD2i/Vehicle)') 

plt.savefig('mostenhancer_leastk36me2_topgenes_wtpm.pdf', bbox_inches='tight')
plt.close()




## plot genes lost most k36me2
toplostgenes = ac_normed_norna.loc[ac_normed_norna['connected_gene'].isin(topgenes[:30])].\
                            melt(id_vars=["genehancer_id", "connected_gene"], var_name="original_col", value_name="value")

toplostgenes['day'] = toplostgenes['original_col'].str.split('_').str[1]
toplostgenes['rep'] = toplostgenes['original_col'].str.split('_').str[2]
toplostgenes['mark'] = toplostgenes['original_col'].str.split('_').str[3]


toploset_expr = pd.merge(toplostgenes,\
                         tpm_long.loc[tpm_long['connected_gene'].isin(topgenes[:30])],\
                         on=['connected_gene', 'day', 'rep'], how='right')
toploset_expr['connected_gene'] = pd.Categorical(toploset_expr['connected_gene'], categories=topgenes[:30])


# plot
g = sns.catplot(data=toploset_expr, kind="box",palette="pastel",
            x="connected_gene", y="value_x", row="mark", col='day', aspect=1.4, 
            showfliers=False, linewidth=1.2)

#https://stackoverflow.com/questions/67309730/how-to-overlay-a-scatterplot-on-top-of-boxplot-with-sns-catplot
g.map_dataframe(sns.stripplot, x="connected_gene", y="value_x", 
                palette=["#404040"],
                alpha=0.6,dodge=True,
                size=3)

g.map_dataframe(sns.stripplot, x="connected_gene", y="value_y", 
                palette=["#CD5555"],
                alpha=0.5,dodge=True,
                s=5, marker="D",
                jitter=False)

for ax in g.axes.flat:
    ax.axhline(0, color='black', linestyle='--', linewidth=0.5)
    
    for label in ax.get_xticklabels():
        label.set_rotation(90)
        label.set_horizontalalignment('center')

g.axes[1, 0].set_ylabel('log2(NSD2i/Vehicle)')

plt.savefig('leastk36me2_topgenes_wtpm.pdf', bbox_inches='tight')
plt.close()
