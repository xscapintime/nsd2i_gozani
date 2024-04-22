import os
import pandas as pd
import seaborn as sns
from matplotlib import pyplot as plt


## peak number
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
    g.set(title=m)
    g.fig.suptitle(m)
    g.set_axis_labels("", "Peak number")
    for axes in g.axes.flat:
        _ = axes.set_xticklabels(axes.get_xticklabels(), rotation=90)

    plt.tight_layout()
    plt.savefig(f'{m}_seacr_peaknum.pdf')
    plt.close()

    # pointplot
    plt.figure(figsize=(4,5))
    g = sns.pointplot(data=dat, x="day", y="num", hue="trt", dodge=True,
                  palette=["#5773CC","#8B8989" ], linestyles=["-", "--"])
    g.get_legend().set_title("")
    g.set_title(m)
    g.set(xlabel="", ylabel="Peak number")
    plt.tight_layout()
    plt.savefig(f'{m}_seacr_pknum_point.pdf')
    plt.close()


## peak average length
stats = pd.read_csv('cutnrun.seacr_noigg.peaklen.txt', sep='\t', header=0)
stats['mark'] = stats['sample'].str.split('_').str[3]
stats['day'] = stats['sample'].str.split('_').str[1]
stats['trt'] = stats['sample'].str.split('_').str[0]


for m in stats['mark'].unique():
    dat = stats.loc[stats['mark']==m, ] 

    # pointplot
    plt.figure(figsize=(4,5))
    g = sns.pointplot(data=dat, x="day", y="ave_peaklen", hue="trt", dodge=True,
                  palette=["#5773CC","#8B8989" ], linestyles=["-", "--"])
    g.get_legend().set_title("")
    g.set_title(m)
    g.set(xlabel="", ylabel="Ave. peak length")
    plt.tight_layout()
    plt.savefig(f'{m}_seacr_pklen_point.pdf')
    plt.close()

