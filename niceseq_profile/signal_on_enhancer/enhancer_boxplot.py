import os,glob
import pandas as pd
import numpy as np
import seaborn as sns
import matplotlib
from matplotlib import pyplot as plt
import matplotlib.font_manager
from statannotations.Annotator import Annotator
from itertools import product


plt.style.use('seaborn-poster')
matplotlib.rcParams['pdf.fonttype'] = 42
matplotlib.rcParams['ps.fonttype'] = 42


## load nice signal on enhancers
nice_enh = pd.read_csv('niceseq_on_enhancer.txt', sep='\t', header=0)

## norm
nsd2i =  nice_enh[sorted(nice_enh.columns[nice_enh.columns.str.contains('NSD2i')])]
nsd2i = nsd2i.groupby(nsd2i.columns.str.extract(r'(.+?)\.rp', expand=False), axis=1).mean()

ck =  nice_enh[sorted(nice_enh.columns[nice_enh.columns.str.contains('vehicle')])]
ck = ck.groupby(ck.columns.str.extract(r'(.+?)\.rp', expand=False), axis=1).mean()

nice_normed = np.log2((nsd2i+1)/(ck+1).values)
# nice_normed = nice_normed[np.sum(nice_normed,axis=1) > 0]

# add gene symbol
nice_normed = pd.merge(nice_enh[['connected_gene']], nice_normed, left_index=True, right_index=True)

# gene symbol as index
nice_normed.set_index('connected_gene', inplace=True)



## load gene sets
gene_sets = sorted(glob.glob(os.path.join('../../path_vsctrl_enhancernum/ctrlsets','*.txt')))

for s in gene_sets:
    gene_set = pd.read_csv(s, sep='\t', header=0)
    pn = os.path.basename(s).replace('_ctrl_genes.txt','')

    path_set = nice_normed[nice_normed.index.isin(gene_set[pn])].stack().to_frame().assign(geneset=f'{pn}'.split('_')[0]).reset_index()
    path_set['refgene'] = path_set['connected_gene']

    ctrl_set1 = nice_normed[nice_normed.index.isin(gene_set.iloc[:,0])].stack().to_frame().assign(geneset='ctrl1').reset_index()
    ctrl_set1['refgene'] = ctrl_set1['connected_gene'].map(dict(zip(gene_set['ctrl1'],gene_set[f'{pn}'])))

    ctrl_set2 = nice_normed[nice_normed.index.isin(gene_set.iloc[:,1])].stack().to_frame().assign(geneset='ctrl2').reset_index()
    ctrl_set2['refgene'] = ctrl_set2['connected_gene'].map(dict(zip(gene_set['ctrl2'],gene_set[f'{pn}'])))

    ctrl_set3 = nice_normed[nice_normed.index.isin(gene_set.iloc[:,2])].stack().to_frame().assign(geneset='ctrl3').reset_index()
    ctrl_set3['refgene'] = ctrl_set3['connected_gene'].map(dict(zip(gene_set['ctrl3'],gene_set[f'{pn}'])))

    ctrl_set4 = nice_normed[nice_normed.index.isin(gene_set.iloc[:,3])].stack().to_frame().assign(geneset='ctrl4').reset_index()
    ctrl_set4['refgene'] = ctrl_set4['connected_gene'].map(dict(zip(gene_set['ctrl4'],gene_set[f'{pn}'])))

    ctrl_set5 = nice_normed[nice_normed.index.isin(gene_set.iloc[:,4])].stack().to_frame().assign(geneset='ctrl5').reset_index()
    ctrl_set5['refgene'] = ctrl_set5['connected_gene'].map(dict(zip(gene_set['ctrl5'],gene_set[f'{pn}'])))


    ## concat all sets
    all_sets = pd.concat([path_set,ctrl_set1,ctrl_set2,ctrl_set3,ctrl_set4,ctrl_set5]).reset_index()
    all_sets['day'] = all_sets['level_1'].str.split('.').str[2]


    # merge 5 control sets
    # `index` represents different enhancers of the same gene
    plot_long_df = pd.concat([all_sets[~all_sets.geneset.str.contains('ctrl')].groupby(['refgene', 'day','index']).mean().reset_index().assign(geneset=f'{pn}'.split('_')[0]),
        all_sets[all_sets.geneset.str.contains('ctrl')].groupby(['connected_gene', 'day','index']).mean().reset_index().assign(geneset='ctrls')])


    # plot
    # plt.figure(figsize=(10, 5))

    g = sns.catplot(data=plot_long_df, kind="violin",palette=["#7CCD7C", "#CDC9C9"],
                x="geneset", y=0, col="day",
                saturation=0.7, linewidth=1, inner='box',aspect=0.7)

    #https://stackoverflow.com/questions/67309730/how-to-overlay-a-scatterplot-on-top-of-boxplot-with-sns-catplot
    g.map_dataframe(sns.stripplot, x="geneset", y=0, 
                    palette=["#404040"], jitter=True,
                   alpha=0.3, dodge=True, size=2)


    # stats
    pairs = [tuple(set(plot_long_df.geneset))]

    ant = Annotator(None, pairs)
    kwargs = {
        'plot_params': { # this takes what normally goes into sns.barplot etc.
            'x': 'geneset',
            'y': 0,
            'palette':["#7CCD7C", "#CDC9C9"]
        },
        'annotation_func': 'apply_test', # has three options
        'configuration': {'test': 't-test_welch', 'text_format' :'full'}, # this takes what normally goes into ant.configure
        'plot': 'violinplot'
    }

    g.map_dataframe(ant.plot_and_annotate_facets, **kwargs)


    for ax in g.axes.flat:
        ax.set_ylabel('log2(NSD2i/Vehicle) NiCE-seq')
        ax.set_xlabel('')

    g.fig.suptitle(f'{pn}', y=1.02, fontsize=16)

    plt.savefig(f'{pn}_nice_onenhancer_violin_mergedctrl.pdf', bbox_inches='tight')
    plt.savefig(f'{pn}_nice_onenhancer_violin_mergedctrl.png', bbox_inches='tight',  dpi=600)
    plt.close()



    ## merge controls by refgene
    plot_long_df_mer = pd.concat([all_sets[~all_sets.geneset.str.contains('ctrl')].groupby(['refgene', 'day']).mean().reset_index().assign(geneset=f'{pn}'.split('_')[0]),
        all_sets[all_sets.geneset.str.contains('ctrl')].groupby(['refgene', 'day']).mean().reset_index().assign(geneset='ctrls')])


    # plot
    # plt.figure(figsize=(10, 5))

    g = sns.catplot(data=plot_long_df_mer, kind="violin",palette=["#7CCD7C", "#CDC9C9"],
                x="geneset", y=0, col="day",
                saturation=0.7, linewidth=1, inner='box', aspect=0.7)

    #https://stackoverflow.com/questions/67309730/how-to-overlay-a-scatterplot-on-top-of-boxplot-with-sns-catplot
    g.map_dataframe(sns.stripplot, x="geneset", y=0, 
                    palette=["#404040"], jitter=True,
                   alpha=0.3, dodge=True, size=2)


    # stats
    pairs = [tuple(set(plot_long_df_mer.geneset))]

    ant = Annotator(None, pairs)
    kwargs = {
        'plot_params': { # this takes what normally goes into sns.barplot etc.
            'x': 'geneset',
            'y': 0,
            'palette':["#7CCD7C", "#CDC9C9"]
        },
        'annotation_func': 'apply_test', # has three options
        'configuration': {'test': 't-test_welch', 'text_format' :'full'}, # this takes what normally goes into ant.configure
        'plot': 'violinplot'
    }

    g.map_dataframe(ant.plot_and_annotate_facets, **kwargs)


    for ax in g.axes.flat:
        ax.set_ylabel('log2(NSD2i/Vehicle) NiCE-seq')
        ax.set_xlabel('')

    g.fig.suptitle(f'{pn}', y=1.02, fontsize=16)

    plt.savefig(f'{pn}_nice_onenhancer_violin_mergedctrl_mergeden.pdf', bbox_inches='tight')
    plt.savefig(f'{pn}_nice_onenhancer_violin_mergedctrl_mergeden.png', bbox_inches='tight',  dpi=600)
    plt.close()
