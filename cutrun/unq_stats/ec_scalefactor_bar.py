import os
import pandas as pd
import seaborn as sns
from matplotlib import pyplot as plt

## load
stats = pd.read_csv('miapaca2_cr_bt2.clean.stats.txt', sep='\t', header=0)

marksinuse = ['H3K36me2', 'H3K36me3', 'H3K27me3', 'H3K27Ac']

stats = stats[stats['samples'].str.contains('|'.join(marksinuse))]
stats['mark'] = stats['samples'].str.split('_').str[3].str.replace('Ac', 'ac')
stats['dayrep'] = stats['samples'].str.split('_H').str[0]
stats['condition'] = stats['samples'].str.split('_').str[0]

# mark order
stats['mark'] = pd.Categorical(stats['mark'], categories=['H3K36me2', 'H3K36me3', 'H3K27me3', 'H3K27ac'], ordered=True)

## constant/#ec
stats['cons_ec'] = 1e6/stats['uniq_ec']


## plot
plt.figure(figsize=(12,8))
g = sns.catplot(data=stats, x='dayrep', y='cons_ec', col='mark', kind='bar',\
                hue='condition', legend=True, height=5.5, aspect=.7, dodge=False,\
                sharey=False, palette={'NSD2i':'#5773CC', 'Vehicle':'#8B8989'})

for axes in g.axes.flat:
    _ = axes.set_xticklabels(axes.get_xticklabels(), rotation=90)

g.set_ylabels('Scaling factors: 1e6/E. coli reads')
g.set_xlabels('')
g.legend.set_title('')

# plt.tight_layout()
plt.savefig('x1e6ec_scalefactor_bar.pdf', bbox_inches='tight')
plt.close()

# export cons_ec
df = stats.iloc[:, [0, 7]]
df['samples'] = df['samples'].str.replace('Ac', 'ac')
df.to_csv('scalingfactors.x1e6ec.csv', sep=',', index=False)
