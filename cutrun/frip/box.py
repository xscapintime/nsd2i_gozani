import os
import pandas as pd
import seaborn as sns
from matplotlib import pyplot as plt


## load
stats = pd.read_csv('nsd2i_frip.txt', sep='\t', header=None)
stats['mark'] = stats[0].str.split('_').str[3]
stats.columns = ['sam', 'all', 'inpeak', 'mark']
stats['day'] = stats['sam'].str.split('_').str[1]
stats['trt'] = stats['sam'].str.split('_').str[0]
stats['drp'] = stats['sam'].str.split('_').str[1:3].str.join('_')

stats['prop'] = stats['inpeak'] / stats['all']

for m in stats['mark'].unique():
    dat = stats.loc[stats['mark']==m, ] 

    g = sns.catplot(data=dat, x='drp', y='prop', hue='day', col='trt', kind='bar',\
                    legend=False, height=5, aspect=.6, dodge=False)
    # g.set(title=m)
    g.fig.suptitle(m)
    g.set_axis_labels("", "FRiP")
    for axes in g.axes.flat:
        _ = axes.set_xticklabels(axes.get_xticklabels(), rotation=90)

    plt.tight_layout()
    plt.savefig(f'{m}_seacr_frip.pdf')
    plt.close()




