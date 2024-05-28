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
ac_normed_norna = ac_normed.loc[least_k36].loc[~ac_normed['connected_gene'].str.contains('lnc|HSALNG|piR-|LOC|ENSG|RF00|LINC0|CM0')]

ac_normed_norna.to_csv('ct_enhancer_normed_leastk36me2_norna.txt', header=True, index=False, sep='\t')


