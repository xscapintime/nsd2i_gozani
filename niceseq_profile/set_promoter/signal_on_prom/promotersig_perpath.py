import os,glob
import pandas as pd
import numpy as np
import seaborn as sns
import matplotlib
from matplotlib import pyplot as plt
import matplotlib.font_manager
from statannotations.Annotator import Annotator


plt.style.use('seaborn-poster')
matplotlib.rcParams['pdf.fonttype'] = 42
matplotlib.rcParams['ps.fonttype'] = 42

paths = ['E2F_TARGETS',
'EPITHELIAL_MESENCHYMAL_TRANSITION',
'HEDGEHOG',
'INTERFERON_GAMMA_RESPONSE',
'KRAS_SIGNALING_DN',
'KRAS_SIGNALING_UP',
'MYC_V1',
'MYC_V2']


for pn in paths:
    
    ## load bed files and aggregate
    bedfiles = sorted(glob.glob(f'*{pn}*.bed'))
    conditions = np.unique([x.split('.sig')[0].replace('MiaPaCa2.','') for x in bedfiles])

    lists_dict = {'.'.join(condition.split('.')[:2]): [] for condition in conditions}

    for c in conditions:
        for b in bedfiles:
            if c in b:
                df = pd.read_csv(b,sep='\t',header=None)
                df.columns = ['chr','start','end','name','score','strand', 'nice']
                df['gene'] = df['name'].str.split('\.').str[0]
                df['trscp'] = df['name'].str.split('\.').str[1]
                # df = df.groupby(['chr','start','end','gene','strand']).mean().reset_index()
                df['condition'] = c
                df['set'] = b.split('.')[5]
                lists_dict['.'.join(c.split('.')[:2])].append(df)

    nsd2i = []
    ck = []
    for l in lists_dict.items():
        df = pd.concat(l[1])
        df['trt'] = df['condition'].str.split('.').str[0]
        df['day'] = df['condition'].str.split('.').str[1]
        df = df.groupby(['chr','start','end','gene','strand','trt','day','set']).mean().reset_index()

        if 'vehicl' in l[0]:
            ck.append(df)
        else:
            nsd2i.append(df)

    nsd2_df = pd.concat(nsd2i)
    ck_df = pd.concat(ck)

    nomred_df = pd.merge(nsd2_df,ck_df,on=['chr','start','end','gene','strand','day','set'],suffixes=('_nsd2','_ck'))
    nomred_df['log2fc'] = np.log2((nomred_df['nice_nsd2']+1)/(nomred_df['nice_ck']+1))


    ## merge 5 controls
    plot_df = pd.concat([nomred_df[~nomred_df.set.str.contains('ctrl')][['chr','start','end','gene','strand','day','log2fc','set']],
                            nomred_df[nomred_df.set.str.contains('ctrl')].groupby(['chr','start','end','gene','strand','day']).mean().reset_index().assign(set='ctrls')[['chr','start','end','gene','strand','day','log2fc','set']]        
                            ])

    ## plot
    # plot
    plt.figure(figsize=(10, 5))
    g = sns.catplot(data=plot_df, kind="violin",palette=["#9F79EE", "#CDC9C9"],
                x="set", y='log2fc', col="day",
                saturation=0.7, linewidth=1, inner='box')
    #https://stackoverflow.com/questions/67309730/how-to-overlay-a-scatterplot-on-top-of-boxplot-with-sns-catplot
    g.map_dataframe(sns.stripplot, x="set", y='log2fc', 
                palette=["#404040"], jitter=True,
                alpha=0.3, dodge=True, size=2)

    # stats
    pairs = [tuple(set(plot_df.set))]
    ant = Annotator(None, pairs)
    kwargs = {
        'plot_params': { # this takes what normally goes into sns.barplot etc.
            'x': 'set',
            'y': 'log2fc',
            'palette':["#9F79EE", "#CDC9C9"]
        },
        'annotation_func': 'apply_test', # has three options
        'configuration': {'test': 't-test_welch', 'text_format' :'full'}, # this takes what normally goes into ant.configure
        'plot': 'violinplot'
    }
    g.map_dataframe(ant.plot_and_annotate_facets, **kwargs)
    for ax in g.axes.flat:
        ax.set_ylabel('log2(NSD2i/Vehicle) NiCE-seq')
        ax.set_xlabel('')
    g.fig.suptitle(f'{pn}', y=1.02, fontsize=12)
    plt.savefig(f'{pn}_nice_onpromoter_violin_mergedctrl.pdf', bbox_inches='tight')
    plt.close()


