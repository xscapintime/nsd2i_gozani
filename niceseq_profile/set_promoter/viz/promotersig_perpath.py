import os,glob
import pandas as pd
import numpy as np
import seaborn as sns
import matplotlib
from matplotlib import pyplot as plt
import matplotlib.font_manager
from statannotations.Annotator import Annotator
from itertools import combinations


plt.style.use('seaborn-poster')
matplotlib.rcParams['pdf.fonttype'] = 42
matplotlib.rcParams['ps.fonttype'] = 42


# load gene sets and ctrls table
gene_sets = sorted(glob.glob(os.path.join('../../../path_vsctrl_enhancernum/ctrlsets','*_ctrl_genes.txt')))

paths = [os.path.basename(x).replace('_ctrl_genes.txt', '') for x in gene_sets]


for pn in paths:
    ## load gene sets and ctrls table
    set_df = pd.read_csv(os.path.join('../../../path_vsctrl_enhancernum/ctrlsets/', f'{pn}_ctrl_genes.txt'), sep='\t')

    ## load bed files and aggregate
    bedfiles = sorted(glob.glob(f'../signal_on_prom/*{pn}*.bed'))
    conditions = np.unique([os.path.basename(x).split('.sig')[0].replace('MiaPaCa2.','') for x in bedfiles])

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
                df['set'] = os.path.basename(b).split('.')[5]
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


    ## merge 5 controls but kept each promoter
    plot_df = pd.concat([nomred_df[~nomred_df.set.str.contains('ctrl')][['chr','start','end','gene','strand','day','log2fc','set','nice_nsd2','nice_ck']],
                            nomred_df[nomred_df.set.str.contains('ctrl')].groupby(['chr','start','end','gene','strand','day']).mean().reset_index().assign(set='ctrls')[['chr','start','end','gene','strand','day','log2fc','set','nice_nsd2','nice_ck']]        
                            ])
    plot_df['set'] = plot_df['set'].str.split('_').str[0]
    
    ## plot
    g = sns.catplot(data=plot_df, kind="violin",palette=["#9F79EE", "#CDC9C9"],
                x="set", y='log2fc', col="day",
                saturation=0.7, linewidth=1, inner='box', aspect=0.7)
    #https://stackoverflow.com/questions/67309730/how-to-overlay-a-scatterplot-on-top-of-boxplot-with-sns-catplot
    g.map_dataframe(sns.stripplot, x="set", y='log2fc', 
                palette=["#404040"], jitter=True,
                alpha=0.2, dodge=True, size=2)

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

    # add mean
    # Iterate over each subplot
    for ax in g.axes.flat:
        # Get the data for this subplot
        data = plot_df[
            (plot_df['day'] == ax.get_title().split('=')[1].strip()) #&
            # (plot_df['trt'] == ax.get_title().split('|')[0].split('=')[1].strip())
        ]

        # Iterate over each category in the x-axis
        for category in data['set'].unique():
            subset = data[data['set'] == category]
            mean_val = subset['log2fc'].mean()

            # Get the position of the category on the x-axis
            x_position = data['set'].unique().tolist().index(category)

            # Draw a short horizontal line at the mean value
            ax.plot([x_position - 0.2, x_position + 0.2], [mean_val, mean_val], color='#FFD700', linestyle=':', linewidth=1.5)

            # Annotate the mean value
            ax.text(
                x=x_position,
                y=mean_val,
                s=f'{mean_val:.2f}',
                color='black',
                ha='center',
                va='bottom',
                fontsize='small'
            )
    
    g.set_axis_labels("", "log2(NSD2i/Vehicle) NiCE-seq")
    g.fig.suptitle(f'{pn}', y=1.02, fontsize=16)
    # plt.tight_layout()
    plt.savefig(f'{pn}_nice_onpromoter_violin_mergedctrl.pdf', bbox_inches='tight')
    # plt.savefig(f'{pn}_nice_onpromoter_violin_mergedctrl.png', bbox_inches='tight', dpi=600)
    plt.close()



    ### merge 5 controls and promoters based ref gene
    mapping_dict = {ctrl: dict(zip(set_df[ctrl], set_df[pn])) for ctrl in set_df.columns[:-1]}
    nomred_df['refgene'] = nomred_df.apply(lambda row: next((mapping_dict[ctrl].get(row['gene'], row['gene']) for ctrl in set_df.columns[:-1] if ctrl in row['set']), row['gene']), axis=1)


    ## merge 5 controls
    # merge promoter of each gene
    # merge promoter of ctrl genes based ref gene
    plot_df_mer = pd.concat([nomred_df[~nomred_df.set.str.contains('ctrl')].groupby(['refgene','day'])[['log2fc','nice_nsd2', 'nice_ck']].mean().reset_index().assign(set=f'{pn}'),
                        nomred_df[nomred_df.set.str.contains('ctrl')].groupby(['refgene','day'])[['log2fc','nice_nsd2', 'nice_ck']].mean().reset_index().assign(set='ctrls')        
                        ])
    plot_df_mer['set'] = plot_df_mer['set'].str.split('_').str[0]
    
        ## plot
    g = sns.catplot(data=plot_df_mer, kind="violin",palette=["#9F79EE", "#CDC9C9"],
                x="set", y='log2fc', col="day",
                saturation=0.7, linewidth=1, inner='box', aspect=0.7)
    #https://stackoverflow.com/questions/67309730/how-to-overlay-a-scatterplot-on-top-of-boxplot-with-sns-catplot
    g.map_dataframe(sns.stripplot, x="set", y='log2fc', 
                palette=["#404040"], jitter=True,
                alpha=0.3, dodge=True, size=2)

    # stats
    pairs = [tuple(set(plot_df_mer.set))]
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


    # add mean
    # Iterate over each subplot
    for ax in g.axes.flat:
        # Get the data for this subplot
        data = plot_df_mer[
            (plot_df_mer['day'] == ax.get_title().split('=')[1].strip()) #&
            # (plot_df['trt'] == ax.get_title().split('|')[0].split('=')[1].strip())
        ]

        # Iterate over each category in the x-axis
        for category in data['set'].unique():
            subset = data[data['set'] == category]
            mean_val = subset['log2fc'].mean()

            # Get the position of the category on the x-axis
            x_position = data['set'].unique().tolist().index(category)

            # Draw a short horizontal line at the mean value
            ax.plot([x_position - 0.2, x_position + 0.2], [mean_val, mean_val], color='#FFD700', linestyle=':', linewidth=1.5)

            # Annotate the mean value
            ax.text(
                x=x_position,
                y=mean_val,
                s=f'{mean_val:.2f}',
                color='black',
                ha='center',
                va='bottom',
                fontsize='small'
            )

    g.set_axis_labels("", "log2(NSD2i/Vehicle) NiCE-seq")


    g.fig.suptitle(f'{pn}', y=1.02, fontsize=16)
    # plt.tight_layout()
    plt.savefig(f'{pn}_nice_onpromoter_violin_mergedctrl_meredpromo.pdf', bbox_inches='tight')
    # plt.savefig(f'{pn}_nice_onpromoter_violin_mergedctrl_meredpromo.png', bbox_inches='tight', dpi=600)
    plt.close()


    # ## plot nsd2i and ctrl separately
    # sep_df = plot_df_mer.melt(id_vars=['refgene', 'day', 'log2fc', 'set'], 
    #                      value_vars=['nice_nsd2', 'nice_ck'], 
    #                      var_name='source', 
    #                      value_name='nice_value')


    # g = sns.catplot(data=sep_df, kind="violin", hue='source', palette=["#9F79EE", "#CDC9C9"],
    #             x="set", y='nice_value', col="day",
    #             saturation=0.7, linewidth=1, inner='box', aspect=0.7,
    #             split=True)
    # #https://stackoverflow.com/questions/67309730/how-to-overlay-a-scatterplot-on-top-of-boxplot-with-sns-catplot
    # g.map_dataframe(sns.stripplot, x="set", y='nice_value', hue='source',
    #             palette=["#404040"], jitter=True,
    #             alpha=0.3, dodge=True, size=2)


    # # stats
    # pairs = list(combinations(list(product(sep_df.set.unique(), sep_df.source.unique())),2))
    # ant = Annotator(None, pairs)
    # kwargs = {
    #     'plot_params': { # this takes what normally goes into sns.barplot etc.
    #         'x': 'set',
    #         'y': 'nice_value',
    #         'palette':["#9F79EE", "#CDC9C9"],
    #         'hue': 'source',
    #         'split': True
    #     },
    #     'annotation_func': 'apply_test', # has three options
    #     'configuration': {'test': 't-test_welch', 'text_format' :'full'}, # this takes what normally goes into ant.configure
    #     'plot': 'violinplot'
    # }
    g.map_dataframe(ant.plot_and_annotate_facets, **kwargs)





