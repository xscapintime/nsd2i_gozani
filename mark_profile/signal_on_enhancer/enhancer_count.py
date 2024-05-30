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