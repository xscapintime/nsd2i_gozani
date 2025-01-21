import os,glob, gzip
import pandas as pd
import numpy as np
import seaborn as sns
from matplotlib import pyplot as plt
import matplotlib.font_manager
from itertools import product
from statannotations.Annotator import Annotator
from matplotlib import ticker as mticker


plt.style.use('seaborn-v0_8-poster')
plt.rcParams['pdf.fonttype'] = 42
plt.rcParams['ps.fonttype'] = 42


def split_and_assign(df, row_split, col_split, row_group, col_group):
    split_df = df.iloc[row_split, col_split].copy()
    split_df['row_group'] = row_group
    split_df['col_group'] = col_group
    return split_df

paths = list(set(sorted([ os.path.basename(m).split('.')[0] for m in glob.glob('../../*.tss.10k.gz')])))

for p in paths:
    # mats = [ f for f in glob.glob(f'../{p}.d*.tss.10k.gz') if 'd9' not in f ]
    mats = glob.glob(f'../../{p}.d*.tss.10k.gz')

    ## load gene sets and ctrls table
    set_df = pd.read_csv(os.path.join('../../../../../path_vsctrl_filtered_noncrna/ctrlsets_include0/', f'{p}_ctrl_genes.txt'), sep='\t')


    ### merge 5 controls and promoters based ref gene
    # mapping_dict = {ctrl: dict(zip(set_df[ctrl], set_df[pn])) for ctrl in set_df.columns[:-1]}
    # mapping_dict = {k: v for sub_dict in mapping_dict.values() for k, v in sub_dict.items()}
    control_to_ref = {}

    # Loop through dataframe to populate the dictionary
    for index, row in set_df.iterrows():
        ref_gene = row[p]
        for control_col in set_df.columns[:-1]:
            control_gene = row[control_col]
            if control_gene not in control_to_ref:
                control_to_ref[control_gene] = set()
            control_to_ref[control_gene].add(ref_gene)

    # Convert sets to lists for easier mapping later
    control_to_ref = {key: list(value) for key, value in control_to_ref.items()}

    # deeptools matrix
    mats_list = []
    for m in mats:
        pn = os.path.basename(m).replace('.tss.10k.gz','')
        day = pn.split('.')[1]

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
        row_groups = [f'{p}', 'Controls']
        col_groups = [f'NSD2i.{day}', f'vehicle.{day}']

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


        mats_list.append(balanced_df)
    
    mat_df_long = pd.concat(mats_list)
    mat_df_long['value'] = mat_df_long.value.astype(float)
    mat_df_long['set'] = pd.Categorical(
                    mat_df_long['set'], 
                    categories=['Controls', mat_df_long['set'].unique()[mat_df_long['set'].unique() != 'Controls'][0]], 
                    ordered=True)
    mat_df_long['day'] = mat_df_long['trt'].str.split('.').str[1]
    mat_df_long['day'] = pd.Categorical(mat_df_long['day'], categories=['d1', 'd5', 'd9'], ordered=True)
    mat_df_long['trt'] = mat_df_long['trt'].str.split('.').str[0]
    mat_df_long['trt'] = pd.Categorical(mat_df_long['trt'], categories=[f'vehicle', f'NSD2i'], ordered=True)



    # longt to wide
    df_pivot = mat_df_long.pivot_table(index=['set', 'refgene','day'], columns='trt', values='value').reset_index()
    df_pivot['diff'] = df_pivot[f'NSD2i'] - df_pivot[f'vehicle']

    df_pivot['set'] = df_pivot['set'].apply(lambda x: x.split('_')[0])
    df_pivot['set'] = pd.Categorical(
                    df_pivot['set'], 
                    categories=['Controls', df_pivot['set'].unique()[df_pivot['set'].unique() != 'Controls'][0]], 
                    ordered=True)

    #     # plot
    # plt.figure(figsize=(4, 5))

    # g = sns.violinplot(data=long_df,x='variable',y='value', inner='box',
    #     palette=["#4682B4", "#CDC9C9"],  saturation=0.7, linewidth=1)
    # g = sns.stripplot(data=long_df, x="variable", y='value', jitter=True,
    #                   color='black', alpha=0.3, dodge=True, size=2)

    # df_mean = long_df.groupby('variable', sort=False)['value'].mean()
    # _ = [g.hlines(y, i-.12, i+.12, zorder=2, colors='#FFD700', linestyle=':', linewidth=1.5) for i, y in df_mean.reset_index()['value'].items()]

    # nobs = ["mean=" + str(f"{i:.2f}") for i in df_mean.to_list()]
    # pos = range(len(nobs))
    # for tick, label in zip(pos, g.get_xticklabels()):
    #    g.text(pos[tick], df_mean[tick] + .05, nobs[tick],
    #             horizontalalignment='center',
    #             size='small',
    #             color='black')
    # plt.ylabel('H3K27me3 (Gene set - Control)')
    # plt.xlabel('')
    # plt.title(f'{pn}', fontdict={'fontsize':7})
    
    # # stats
    # comb = [tuple(set(long_df.variable))]
    # annotator = Annotator(g, comb, data=long_df, x='variable', y='value')
    # annotator.configure(test="Wilcoxon", text_format="full",loc='inside')
    # annotator.apply_and_annotate()
    # plt.savefig(f'{pn}_cprdiff_tss_central4.8kbp_k27me3.pdf', bbox_inches='tight')


    ## plot
    # g = sns.catplot(data=mat_df_long, kind="violin",palette={'geneset':"#4682B4","ctrlset":"#CDC9C9"},
    #             x='variable', y='value', col='day',
    #             saturation=0.7, linewidth=.3, inner='box',
    #             aspect=.8)
    # # g.map_dataframe(sns.stripplot, x="variable", y="value", 
    # #             palette=["black"],
    # #             alpha=0.3,dodge=True,
    # #             size=3)


    # # stats
    # pairs = [tuple(set(mat_df_long.variable))]

    # ant = Annotator(None, pairs)
    # kwargs = {
    #     'plot_params': { # this takes what normally goes into sns.barplot etc.
    #         'x': 'variable',
    #         'y': 'value',
    #         'palette':{'geneset':"#4682B4","ctrlset":"#CDC9C9"}
    #     },
    #     'annotation_func': 'apply_test', # has three options
    #     'configuration': {'test': 'Wilcoxon', 'text_format' :'full'}, # this takes what normally goes into ant.configure
    #     'plot': 'violinplot'
    # }

    # g.map_dataframe(ant.plot_and_annotate_facets, **kwargs)


    # # add mean
    # # Iterate over each subplot
    # for ax in g.axes.flat:
    #     # Get the data for this subplot
    #     data = mat_df_long[
    #         (mat_df_long['day'] == ax.get_title().split('=')[1].strip()) #&
    #         # (plot_df['trt'] == ax.get_title().split('|')[0].split('=')[1].strip())
    #     ]

    #     # Iterate over each category in the x-axis
    #     for category in  data['variable'].cat.categories:
    #         subset = data[data['variable'] == category]
    #         mean_val = subset['value'].mean()

    #         # Get the position of the category on the x-axis
    #         x_position = data['variable'].cat.categories.tolist().index(category)

    #         # Draw a short horizontal line at the mean value
    #         ax.plot([x_position - 0.2, x_position + 0.2], [mean_val, mean_val], color='#FFD700', linestyle=':', linewidth=1.5)

    #         # Annotate the mean value
    #         ax.text(
    #             x=x_position,
    #             y=mean_val,
    #             s=f'mean={mean_val:.2f}',
    #             color='black',
    #             ha='center',
    #             va='bottom',
    #             fontsize='small'
    #         )

    # g.set_ylabels('H3K27me3 (NSD2i - Control)')
    # g.set_xlabels('')

    # g.fig.suptitle(f'{p}', y=1.02, fontsize=12)

    # plt.savefig(f'{p}_nsd2subveh_tss_central4.8kbp_k27me3.pdf', bbox_inches='tight')
    # plt.close()



    ### do a tomato and blue palette plot
    g = sns.catplot(data=df_pivot, kind="violin",palette=["#edb9aa", "#7c84f4"],
                x='set', y='diff', col='day',
                saturation=0.7, linewidth=.3, inner='box',
                aspect=.8)
    # g.map_dataframe(sns.stripplot, x="variable", y="value", 
    #             palette=["black"],
    #             alpha=0.3,dodge=True,
    #             size=3)

    # stats
    pairs = [tuple(set(df_pivot.set))]
    

    ant = Annotator(None, pairs)
    kwargs = {
        'plot_params': { # this takes what normally goes into sns.barplot etc.
            'x': 'set',
            'y': 'diff',
            'palette':["#edb9aa", "#7c84f4"]
        },
        'annotation_func': 'apply_test', # has three options
        'configuration': {'test': 'Wilcoxon', 'text_format' :'full'}, # this takes what normally goes into ant.configure
        'plot': 'violinplot'
    }

    g.map_dataframe(ant.plot_and_annotate_facets, **kwargs)


    # add mean
    # Iterate over each subplot
    for ax in g.axes.flat:
        # Get the data for this subplot
        data = df_pivot[
            (df_pivot['day'] == ax.get_title().split('=')[1].strip()) #&
            # (plot_df['trt'] == ax.get_title().split('|')[0].split('=')[1].strip())
        ]

        # Iterate over each category in the x-axis
        for category in  data['set'].cat.categories:
            subset = data[data['set'] == category]
            mean_val = subset['diff'].mean()

            # Get the position of the category on the x-axis
            x_position = data['set'].cat.categories.tolist().index(category)

            # Draw a short horizontal line at the mean value
            ax.plot([x_position - 0.2, x_position + 0.2], [mean_val, mean_val], color='black', linestyle=':', linewidth=1.5)

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

    g.set_ylabels('H3K27me3 (NSD2i - Control)')
    g.set_xlabels('')
    g.set(ylim=(-10, 73))


    g.fig.suptitle(f'{p}', y=1.02, fontsize=12)

    plt.savefig(f'{p}_delta_tsscentral4.8kbp_k27me3.pdf', bbox_inches='tight')
    plt.close()
