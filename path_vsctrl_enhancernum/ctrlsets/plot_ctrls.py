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


## load tpm
## load ensemble id to symbol matching table
ensembl_syb = pd.read_csv('../../rnaseq/umap/human_ensembl_syb.tsv', header=0, index_col=None, sep='\t')
ensembl_syb = dict(zip(ensembl_syb['ensembl_gene_id'], ensembl_syb['hgnc_symbol']))

## load tpm table
tpm = pd.read_csv('../../rnaseq/tpm/tximport-tpm.csv', header=0, index_col=0)

## change index to gene symbol
tpm = tpm.set_index(tpm.index.map(ensembl_syb))
tpm = tpm[tpm.index.notnull()]

# merge transcripts, only ~100 genes with more than 1 transcript
tpm = tpm.groupby(tpm.index).mean()

## log2tpm
tpm_normed = np.log2(np.divide(*(tpm.iloc[:,np.where(tpm.columns.str.contains('_N'))[0]]+1).\
                    align((tpm.iloc[:,np.where(tpm.columns.str.contains('_C'))[0]]+1), axis=0)))

tpm_normed.columns = tpm_normed.columns.str.replace('N', 'Rep').str.replace('Day', 'D') + '_log2tpm'



## load ctrl gene tables

ctrl_files = glob.glob('*_ctrl_genes.txt')

pathname = [ f.replace('_ctrl_genes.txt', '') for f in ctrl_files ]

for pn in pathname:
# ['KRAS_SIGNALING_UP', 'KRAS_SIGNALING_DN', 'EPITHELIAL_MESENCHYMAL_TRANSITION', 'INTERFERON_GAMMA_RESPONSE', 'E2F_TARGETS',
        #    'HEDGEHOG', 'MYC_V1', 'MYC_V2', 'PRC2_targ']:

    # path = eval(pn)
    ctrl_genes = pd.read_csv(f'{pn}_ctrl_genes.txt', sep='\t', header=0)
    path = list(ctrl_genes[f'{pn}'])

    ## check raw expression
    path_set_raw = tpm.loc[path].stack().to_frame().assign(geneset=f'{pn}'.split('_')[0])
    ctrl_set1_raw = tpm.loc[ctrl_genes.iloc[:,0]].stack().to_frame().assign(geneset='ctrl1')
    ctrl_set2_raw = tpm.loc[ctrl_genes.iloc[:,1]].stack().to_frame().assign(geneset='ctrl2')
    ctrl_set3_raw = tpm.loc[ctrl_genes.iloc[:,2]].stack().to_frame().assign(geneset='ctrl3')
    ctrl_set4_raw = tpm.loc[ctrl_genes.iloc[:,3]].stack().to_frame().assign(geneset='ctrl4')
    ctrl_set5_raw = tpm.loc[ctrl_genes.iloc[:,4]].stack().to_frame().assign(geneset='ctrl5')    

    ## concat all sets
    all_sets_raw = pd.concat([path_set_raw,ctrl_set1_raw,ctrl_set2_raw,ctrl_set3_raw,ctrl_set4_raw,ctrl_set5_raw]).reset_index()
    all_sets_raw['day'] = all_sets_raw['level_1'].str.split('_').str[0]
    all_sets_raw['rep'] = all_sets_raw['level_1'].str.split('_').str[1].str[1]
    all_sets_raw['trt'] = all_sets_raw['level_1'].str.split('_').str[1].str[0]
    all_sets_raw['log2'] = np.log2(all_sets_raw[0]+1)


    # box log2 tpm +1 
    g = sns.catplot(data=all_sets_raw, kind="violin",
                palette=["#4876FF", "#CDC9C9", "#CDC9C9", "#CDC9C9","#CDC9C9","#CDC9C9"],
                x="geneset", y='log2', col="day",row= 'trt',aspect=1, 
                linewidth=1)
    g.set_axis_labels("", "log2(TPM+1)")
    g.fig.suptitle(f'{pn}', y=1.02, fontsize=16)


    # stats
    pairs = list(product([pn.split('_')[0]], [s for s in np.unique(all_sets_raw.geneset) if 'ctrl' in s]))
    ant = Annotator(None, pairs)
    kwargs = {
        'plot_params': { # this takes what normally goes into sns.barplot etc.
            'x': 'geneset',
            'y': 'log2',
            'palette':["#4876FF", "#CDC9C9", "#CDC9C9", "#CDC9C9","#CDC9C9","#CDC9C9"]
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
        data = all_sets_raw[
            (all_sets_raw['day'] == ax.get_title().split('|')[1].split('=')[1].strip()) &
            (all_sets_raw['trt'] == ax.get_title().split('|')[0].split('=')[1].strip())
        ]

        # Iterate over each category in the x-axis
        for category in data['geneset'].unique():
            subset = data[data['geneset'] == category]
            mean_val = subset['log2'].mean()

            # Get the position of the category on the x-axis
            x_position = data['geneset'].unique().tolist().index(category)

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

    g.set_axis_labels("", "log2(TPM+1)")
    plt.savefig(f'{pn}_log2tpm_from_D1.pdf', bbox_inches='tight')
    plt.close()




    # no log2 transformation
    g = sns.catplot(data=all_sets_raw, kind="violin",
                palette=["#4876FF", "#CDC9C9", "#CDC9C9", "#CDC9C9","#CDC9C9","#CDC9C9"],
                x="geneset", y=0, col="day",row= 'trt',aspect=0.8, 
                linewidth=1.2)
    g.set_axis_labels("", "TPM")
    g.fig.suptitle(f'{pn}', y=1.02, fontsize=16)

    # stats
    pairs = list(product([pn.split('_')[0]], [s for s in np.unique(all_sets_raw.geneset) if 'ctrl' in s]))
    ant = Annotator(None, pairs)
    kwargs = {
        'plot_params': { # this takes what normally goes into sns.barplot etc.
            'x': 'geneset',
            'y': all_sets_raw[0],
            'palette':["#4876FF", "#CDC9C9", "#CDC9C9", "#CDC9C9","#CDC9C9","#CDC9C9"]
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
        data = all_sets_raw[
            (all_sets_raw['day'] == ax.get_title().split('|')[1].split('=')[1].strip()) &
            (all_sets_raw['trt'] == ax.get_title().split('|')[0].split('=')[1].strip())
        ]

        # Iterate over each category in the x-axis
        for category in data['geneset'].unique():
            subset = data[data['geneset'] == category]
            mean_val = subset[0].mean()

            # Get the position of the category on the x-axis
            x_position = data['geneset'].unique().tolist().index(category)

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

    g.set_axis_labels("", "TPM")
    plt.savefig(f'{pn}_tpm_from_D1.pdf', bbox_inches='tight')
    plt.close()



    ## box plot of kras dn vs 5 control sets using nomred logfc tpm
    path_set = tpm_normed.loc[path].stack().to_frame().assign(geneset=f'{pn}'.split('_')[0])
    ctrl_set1 = tpm_normed.loc[ctrl_genes.iloc[:,0]].stack().to_frame().assign(geneset='ctrl1')
    ctrl_set2 = tpm_normed.loc[ctrl_genes.iloc[:,1]].stack().to_frame().assign(geneset='ctrl2')
    ctrl_set3 = tpm_normed.loc[ctrl_genes.iloc[:,2]].stack().to_frame().assign(geneset='ctrl3')
    ctrl_set4 = tpm_normed.loc[ctrl_genes.iloc[:,3]].stack().to_frame().assign(geneset='ctrl4')
    ctrl_set5 = tpm_normed.loc[ctrl_genes.iloc[:,4]].stack().to_frame().assign(geneset='ctrl5')


    ## concat all sets
    all_sets = pd.concat([path_set,ctrl_set1,ctrl_set2,ctrl_set3,ctrl_set4,ctrl_set5]).reset_index()
    all_sets['day'] = all_sets['level_1'].str.split('_').str[0]
    all_sets['rep'] = all_sets['level_1'].str.split('_').str[1]


    ## plot the violin
    g = sns.catplot(data=all_sets, kind="violin",
                palette=["#4876FF", "#CDC9C9", "#CDC9C9", "#CDC9C9","#CDC9C9","#CDC9C9"],
                x="geneset", y=0, col="day", aspect=0.8, 
                linewidth=1.2)
    g.set_axis_labels("", "log2(NSD2i/Vehicle) TPM")
    g.fig.suptitle(f'{pn}', y=1.02, fontsize=16)

    plt.savefig(f'{pn}_ctrl_from_D1.pdf', bbox_inches='tight')
    plt.close()


    ## count enhancers
    count_df = pd.DataFrame()

    for column in ctrl_genes.columns:
        count_df[column] = ctrl_genes[column].map(enhancer_count).fillna(0)

    count_df.index = ctrl_genes[f'{pn}']
    count_df = count_df[[f'{pn}','ctrl1', 'ctrl2', 'ctrl3', 'ctrl4', 'ctrl5']]


    ## combine 5 controls
    plot_df_mergedcrtl = pd.concat([count_df[[f'{pn}']],
                                    count_df.loc[:,count_df.columns.str.contains('ctrl')].\
                                        mean(axis=1)],axis=1).stack().reset_index()

    plot_df_mergedcrtl.columns = ['gene', 'set', 'enhancer_n']
    plot_df_mergedcrtl.set = plot_df_mergedcrtl.set.replace({0: 'ctrl'})
    plot_df_mergedcrtl.set = plot_df_mergedcrtl.set.str.split('_').str[0]



    ## plot the violin plot
    plt.figure(figsize=(5, 5))
    g = sns.violinplot(data=plot_df_mergedcrtl, x='set', y='enhancer_n',
                   palette = ["#4876FF", "#CDC9C9", "#CDC9C9", "#CDC9C9","#CDC9C9","#CDC9C9"],
                   saturation=0.7, linewidth=1, inner='box')

    g = sns.stripplot(data=plot_df_mergedcrtl, x="set", y='enhancer_n', jitter=True,
                  color='black', alpha=0.3, dodge=True, size=2)

    df_mean = plot_df_mergedcrtl.groupby('set', sort=False)['enhancer_n'].mean()
    _ = [g.hlines(y, i-.12, i+.12, zorder=2, colors='white', linewidth=2) for i, y in df_mean.reset_index()['enhancer_n'].items()]

    nobs = ["mean: " + str(f"{i:.2f}") for i in df_mean.to_list()]

    pos = range(len(nobs))
    for tick, label in zip(pos, g.get_xticklabels()):
       g.text(pos[tick], df_mean[tick] + 1, nobs[tick],
                horizontalalignment='center',
                size='small',
                color='black')


    plt.ylabel('Enhancer number per gene')
    plt.xlabel('')
    plt.title(f'{pn}')

    # stats
    comb = [tuple(set(plot_df_mergedcrtl.set))]

    annotator = Annotator(g, comb, data=plot_df_mergedcrtl, x='set', y='enhancer_n')
    annotator.configure(test="t-test_paired", text_format="full",loc='inside')
    annotator.apply_and_annotate()

    plt.savefig(f'{pn}_enhancer_count_violin_mergedctrl.pdf', bbox_inches='tight')



