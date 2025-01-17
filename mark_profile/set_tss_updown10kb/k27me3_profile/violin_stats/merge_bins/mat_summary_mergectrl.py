import os,glob, gzip
import math
import pandas as pd
import numpy as np
import seaborn as sns
import matplotlib
from matplotlib import pyplot as plt
import matplotlib.font_manager
from itertools import product
from statannotations.Annotator import Annotator
from matplotlib import ticker as mticker


# plt.style.use('seaborn-poster')
matplotlib.rcParams['pdf.fonttype'] = 42
matplotlib.rcParams['ps.fonttype'] = 42




def split_and_assign(df, row_split, col_split, row_group, col_group):
    split_df = df.iloc[row_split, col_split].copy()
    split_df['row_group'] = row_group
    split_df['col_group'] = col_group
    return split_df


for m in glob.glob('../../*.tss.10k.gz'):
    pn = os.path.basename(m).replace('.tss.10k.gz','').split('.')[0]
    day = os.path.basename(m).replace('.tss.10k.gz','').split('.')[1]

    ## load gene sets and ctrls table
    set_df = pd.read_csv(os.path.join('../../../../../path_vsctrl_filtered_noncrna/ctrlsets_include0/', f'{pn}_ctrl_genes.txt'), sep='\t')


    ### merge 5 controls and promoters based ref gene
    # mapping_dict = {ctrl: dict(zip(set_df[ctrl], set_df[pn])) for ctrl in set_df.columns[:-1]}
    # mapping_dict = {k: v for sub_dict in mapping_dict.values() for k, v in sub_dict.items()}
    control_to_ref = {}

    # Loop through dataframe to populate the dictionary
    for index, row in set_df.iterrows():
        ref_gene = row[pn]
        for control_col in set_df.columns[:-1]:
            control_gene = row[control_col]
            if control_gene not in control_to_ref:
                control_to_ref[control_gene] = set()
            control_to_ref[control_gene].add(ref_gene)

    # Convert sets to lists for easier mapping later
    control_to_ref = {key: list(value) for key, value in control_to_ref.items()}


    ## load deeptools matrix
    df = pd.read_csv(m, header=None, skiprows=1,sep='\t')
    df[3] = df[3].str.split('.').str[0]

    ## meta data
    with gzip.open(m, 'rb') as f:
        meta = f.readline()

    group_boundaries = int(meta.decode("utf-8").split('":')[20].split(',"')[0].split(',')[1])
    sample_boundaries = int(meta.decode("utf-8").split('":')[22].split(',')[1])

    # Define splits
    row_splits = [slice(None, group_boundaries), slice(group_boundaries, None)]
    col_splits = [slice(None, sample_boundaries), slice(sample_boundaries, None)]
    row_groups = [f'{pn}', 'Controls']
    col_groups = [f'NSD2i.day_{day}', f'vehicle.day_{day}']

    split_dfs = [
        split_and_assign(df.iloc[:, 6:], row_split, col_split, row_group, col_group)
        for row_split, row_group in zip(row_splits, row_groups)
        for col_split, col_group in zip(col_splits, col_groups)
        ]
    
    for i, split_df in enumerate(split_dfs):
        split_df['gene'] = df.iloc[split_df.index, 3]

        if split_df['row_group'].unique()[0] == 'Controls':
            split_df['refgene'] = split_df.gene.map(control_to_ref)

            # Explode the DataFrame and update the corresponding entry in split_dfs
            split_dfs[i] = split_df.explode('refgene', ignore_index=True)
        else:
            split_df['refgene'] = split_df['gene']

        print(split_dfs[i]['refgene'].nunique())

    # Concatenate all parts
    df_combined = pd.DataFrame(np.vstack([pd.concat([split_dfs[0], split_dfs[2]], axis=0),
                            pd.concat([split_dfs[1], split_dfs[3]], axis=0).values]))
    
    # average the center 2.5kb，[:, 38:63] 25bins
    # this will take the mean of all the genes 
    # df_center = pd.concat([df_combined.iloc[:, 38:62].mean(axis=1), df_combined.iloc[:,100:]], axis=1)

    ## average by genes in the geneset (refgene)
    df_exploded = df_combined.drop(columns=[102])
    
    df_center = pd.concat([df_exploded.iloc[:, 38:63], df_exploded.iloc[:,100:]], axis=1).groupby([100,101,103]).mean().reset_index()

    # mean of all the bins
    df_center = df_center.groupby([100, 101, 103]).apply(lambda group: group.iloc[:, 3:].mean().mean()).reset_index(name='mean') 

    # deeptools or extracting promoter might missed some genes
    min_gene_count = df_center.groupby([100, 101])[103].nunique().min()

    balanced_df = (
        df_center.groupby([100, 101])
        .apply(lambda group: group[group[103].isin(group[103].unique()[:min_gene_count])])
        .reset_index(drop=True)
    )


    balanced_df = balanced_df.rename(columns={100: 'set', 101: 'trt',
                                            103: 'refgene', 'mean': 'value'})

    balanced_df['set'] = balanced_df['set'].apply(lambda x: x.split('_')[0])
    balanced_df['set'] = pd.Categorical(balanced_df['set'], categories=['Controls', balanced_df['set'].unique()[balanced_df['set'].unique() != 'Controls'][0]], ordered=True)

    # skip zero
    # long_df = long_df[long_df['value'] != 0]
    # long_df['value1p'] = long_df['value'] + 1
    balanced_df['logvalue'] = np.log10(balanced_df['value'] + 0.001)


    # ## to save the tables
    # # remove the marker for the duplicated control genes
    # df_center[103] = df_center[103].str.split('.').str[0]
    # df_center.to_csv(f'{pn}.{day}.center500bp.refgene.txt', sep='\t', index=False)



    ## log transfrom

    for stat_me in ['Mann-Whitney', 'Wilcoxon']:
        g = sns.catplot(data=balanced_df, kind="violin",palette=["#edb9aa", "#7c84f4"],
                    x='set', y='logvalue', col='trt',
                    saturation=0.7, linewidth=.1, inner='box',
                    aspect=.8)
        # stats
        pairs = [tuple(set(balanced_df.set))]

        ant = Annotator(None, pairs)
        kwargs = {
            'plot_params': { # this takes what normally goes into sns.barplot etc.
                'x': 'set',
                'y': 'logvalue',
                'palette':["#edb9aa", "#7c84f4"]
            },
            'annotation_func': 'apply_test', # has three options
            'configuration': {'test': stat_me, 'text_format' :'full'}, # this takes what normally goes into ant.configure
            'plot': 'violinplot'
        }

        g.map_dataframe(ant.plot_and_annotate_facets, **kwargs)

        g.set(ylim=(-3.55,2.99))

        for ax in g.axes.flat:
            ax.yaxis.set_major_formatter(mticker.StrMethodFormatter("$10^{{{x:.0f}}}$"))
            ymin, ymax = ax.get_ylim()
            # tick_range = np.arange(np.floor(ymin), ymax)
            tick_range = np.arange(ymin, ymax)
            # ax.yaxis.set_ticks(tick_range)
            ax.yaxis.set_ticks([np.log10(x) for p in tick_range for x in np.linspace(10 ** p, 10 ** (p + 1), 10)], minor=True)

        # add mean
        # Iterate over each subplot
        for ax in g.axes.flat:
            # Get the data for this subplot
            data = balanced_df[
                (balanced_df['trt'] == ax.get_title().split('=')[1].strip()) #&
                # (plot_df['trt'] == ax.get_title().split('|')[0].split('=')[1].strip())
            ]

            # Iterate over each category in the x-axis
            for category in data['set'].unique():
                subset = data[data['set'] == category]
                mean_val = subset['value'].mean()
                logmean_val = subset['logvalue'].mean()

                # Get the position of the category on the x-axis
                x_position = data['set'].unique().tolist().index(category)

                # Draw a short horizontal line at the mean value
                ax.plot([x_position - 0.2, x_position + 0.2], [logmean_val, logmean_val],
                    color='black', linestyle=':', linewidth=1.5)

                # Annotate the mean value
                ax.text(
                    x=x_position,
                    y=logmean_val,
                    s=f'mean={mean_val:.2f}',
                    color='black',
                    ha='center',
                    va='bottom',
                    fontsize='small'
                )

        g.set_ylabels('H3K27me3 (log scale)')
        g.set_xlabels('')

        g.fig.suptitle(f'{pn}'.split('.')[0], y=1.02, fontsize=12)

        plt.savefig(f'{pn}_{day}_tsscentral2.4kbp_k27me3_log_{stat_me}.pdf', bbox_inches='tight')
        plt.close()


        # no log transfrom
        g = sns.catplot(data=balanced_df, kind="violin",palette=["#edb9aa", "#7c84f4"],
                    x='set', y='value', col='trt',
                    saturation=0.7, linewidth=.1, inner='box',
                    aspect=.8)
        # stats
        pairs = [tuple(set(balanced_df.set))]

        ant = Annotator(None, pairs)
        kwargs = {
            'plot_params': { # this takes what normally goes into sns.barplot etc.
                'x': 'set',
                'y': 'value',
                'palette':["#edb9aa", "#7c84f4"]
            },
            'annotation_func': 'apply_test', # has three options
            'configuration': {'test': stat_me, 'text_format' :'full'}, # this takes what normally goes into ant.configure
            'plot': 'violinplot'
        }

        g.map_dataframe(ant.plot_and_annotate_facets, **kwargs)

        # g.set(ylim=(-3.55,2.99))

        # add mean
        # Iterate over each subplot
        for ax in g.axes.flat:
            # Get the data for this subplot
            data = balanced_df[
                (balanced_df['trt'] == ax.get_title().split('=')[1].strip()) #&
                # (plot_df['trt'] == ax.get_title().split('|')[0].split('=')[1].strip())
            ]

            # Iterate over each category in the x-axis
            for category in data['set'].unique():
                subset = data[data['set'] == category]
                mean_val = subset['value'].mean()
                # logmean_val = subset['logvalue'].mean()

                # Get the position of the category on the x-axis
                x_position = data['set'].unique().tolist().index(category)

                # Draw a short horizontal line at the mean value
                ax.plot([x_position - 0.2, x_position + 0.2], [mean_val, mean_val],
                    color='black', linestyle=':', linewidth=1.5)

                # Annotate the mean value
                ax.text(
                    x=x_position,
                    y=mean_val,
                    s=f'mean={mean_val:.2f}',
                    color='black',
                    ha='center',
                    va='bottom',
                    fontsize='small'
                )

        g.set_ylabels('H3K27me3')
        g.set_xlabels('')

        g.fig.suptitle(f'{pn}'.split('.')[0], y=1.02, fontsize=12)

        plt.savefig(f'{pn}_{day}tsscentral2.4kbp_k27me3_{stat_me}.pdf', bbox_inches='tight')
        plt.close()
