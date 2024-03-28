import os
import pandas as pd
import seaborn as sns
from matplotlib import pyplot as plt


## load
stats = pd.read_csv('cutnrun.seacr_noigg.stats.txt', sep='\t', header=None)
stats = stats[:-1]
stats['mark'] = stats[1].str.split('_').str[3]
stats.columns = ['num', 'sam', 'mark']
stats['day'] = stats['sam'].str.split('_').str[1]
stats['trt'] = stats['sam'].str.split('_').str[0]
stats['drp'] = stats['sam'].str.split('_').str[1:3].str.join('_')

for m in stats['mark'].unique():
    dat = stats.loc[stats['mark']==m, ] 

    g = sns.catplot(data=dat, x='drp', y='num', hue='day', col='trt', kind='bar',\
                    legend=False, height=5, aspect=.6)
    # g.set(title=m)
    g.fig.suptitle(m)
    g.set_axis_labels("", "Peak number")
    for axes in g.axes.flat:
        _ = axes.set_xticklabels(axes.get_xticklabels(), rotation=90)

    plt.tight_layout()
    plt.savefig(f'{m}_seacr_peaknum.pdf')
    plt.close()




