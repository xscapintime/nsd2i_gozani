import os
import pandas as pd
import seaborn as sns
from matplotlib import pyplot as plt


## load
stats = pd.read_csv('peaknum.stats', sep='\t', header=None)
stats[1] = stats[1].str.split('_peaks').str[0]
stats.columns = ['num', 'sam']
stats['day'] = stats['sam'].str.split('_').str[1:3].apply(lambda x:'_'.join(x))


g = sns.barplot(data=stats, y='num', x='sam', hue='day',dodge=False)
g.get_legend().remove()
g.set(xlabel="", ylabel="")
g.set_xticklabels(g.get_xticklabels(), 
                          rotation=330, 
                          rotation_mode='anchor',ha='left', va='top')

plt.tight_layout()
plt.savefig('niceseq_peaknum.pdf')
plt.close()