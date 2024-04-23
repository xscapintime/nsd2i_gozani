import os
import pandas as pd
import numpy as np
import seaborn as sns
from matplotlib import pyplot as plt
from itertools import product



## peak number
pknum = pd.read_csv('cutnrun.seacr_noigg.stats.txt', sep='\t', header=None)
pknum = pknum[:-1]

## peak average length
pklen = pd.read_csv('cutnrun.seacr_noigg.peaklen.txt', sep='\t', header=0)

## FRiP
frip = pd.read_csv('../frip/nsd2i_frip.txt', sep='\t', header=None)


## merge
merged_dat = pknum.merge(pklen, left_on=1, right_on='sample', how='inner').merge(frip, left_on=1, right_on=0, how='inner')
merged_dat = merged_dat[[1, '0_x', 'ave_peaklen', '1_y', 2]]
merged_dat.columns = ['sample', 'pknum', 'ave_peaklen', 'all', 'inpeak']
merged_dat['frip'] = merged_dat['inpeak'] / merged_dat['all']

merged_dat.set_index('sample', inplace=True)

df1 = merged_dat.iloc[np.where(merged_dat.index.str.contains('NSD2i'))]
df2 = merged_dat.iloc[np.where(merged_dat.index.str.contains('Vehicle'))]

normed = np.divide(*df1.align(df2, axis=1))


## wide to long
normed_long = normed.reset_index().melt(id_vars='sample', var_name='metric', value_name='NSD2i/Vehicle')

normed_long['mark'] = normed_long['sample'].str.split('_').str[3]
normed_long['day'] = normed_long['sample'].str.split('_').str[1]
normed_long['trt'] = normed_long['sample'].str.split('_').str[0]
normed_long['rep'] = normed_long['sample'].str.split('_').str[2]


## plot
# remove read number
normed_long = normed_long.loc[(normed_long['metric']!='all') & (normed_long['metric']!='inpeak'), ]
normed_long['metric'] = normed_long['metric'].map({'pknum': 'Peak num', 'ave_peaklen': 'Ave. Peak len', 'frip': 'FRiP'})

for m in normed_long['mark'].unique():
    dat = normed_long.loc[normed_long['mark']==m, ] 

    # pointplot
    plt.figure(figsize=(5,4))
    g = sns.pointplot(data=dat, x="day", y="NSD2i/Vehicle", hue="metric", dodge=True,
                        markers=['o', 's', 'x'])
    g.get_legend().set_title("")
    g.set_title(m)
    g.set(xlabel="", ylabel="NSD2i/Vehicle")
    g.axhline(1, ls='--', c='grey')
    plt.tight_layout()
    plt.savefig(f'{m}_seacr_allstats.pdf')
    plt.close()