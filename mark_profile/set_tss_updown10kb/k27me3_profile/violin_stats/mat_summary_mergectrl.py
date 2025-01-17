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



# def append_suffix_to_duplicates(column):
#     # Create a dictionary to track counts for each duplicate value
#     value_counts = {}
    
#     # Iterate over each value in the column
#     for idx, value in column.items():
#         if value not in value_counts:
#             value_counts[value] = 0  # Initialize the count
        
#         value_counts[value] += 1  # Increment the count
        
#         # If it's a duplicate (count > 1), append the suffix
#         if value_counts[value] > 1:
#             column.at[idx] = f"{value}.{value_counts[value] - 1}"
    
#     return column


def split_and_assign(df, row_split, col_split, row_group, col_group):
    split_df = df.iloc[row_split, col_split].copy()
    split_df['row_group'] = row_group
    split_df['col_group'] = col_group
    return split_df


for m in glob.glob('../*.tss.10k.gz'):
    pn = os.path.basename(m).replace('.tss.10k.gz','').split('.')[0]
    day = os.path.basename(m).replace('.tss.10k.gz','').split('.')[1]

    ## load gene sets and ctrls table
    set_df = pd.read_csv(os.path.join('../../../../path_vsctrl_filtered_noncrna/ctrlsets_include0/', f'{pn}_ctrl_genes.txt'), sep='\t')
    
    # for col in set_df.columns:
    #     # Identify all duplicates including the first occurrence
    #     duplicated_indices = set_df.index[set_df[col].duplicated(keep=False)]

    #     if not duplicated_indices.empty:
    #         # Dictionary to keep track of counts for each duplicate value
    #         duplicate_counts = {}

    #         for idx in duplicated_indices:
    #             value = set_df.loc[idx, col]

    #             # Initialize the count for this value if not already tracked
    #             if value not in duplicate_counts:
    #                 duplicate_counts[value] = 0

    #             # Increment the count for this value
    #             duplicate_counts[value] += 1

    #             # Append the count to the original string value
    #             set_df.loc[idx, col] = f"{value}.{duplicate_counts[value]}"



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

    for split_df in split_dfs:
        split_df['gene'] = df.iloc[split_df.index, 3]

        if split_df['row_group'].unique()[0] == 'Controls':
            split_df['refgene'] = split_df.gene.map(control_to_ref)
            # split_df['refgene'] = split_df['refgene'].apply(lambda refs: ', '.join(refs))

        else:
            split_df['refgene'] = split_df['gene']
        

    # Concatenate all parts
    df_combined = pd.DataFrame(np.vstack([pd.concat([split_dfs[0], split_dfs[2]], axis=0),
                            pd.concat([split_dfs[1], split_dfs[3]], axis=0).values]))
    
    # average the center 2.5kb，[:, 38:63] 25bins
    # this will take the mean of all the genes 
    # df_center = pd.concat([df_combined.iloc[:, 38:62].mean(axis=1), df_combined.iloc[:,100:]], axis=1)

    ## average by genes in the geneset (refgene)
    # explode the dataframe to create a row for each reference gene
    df_exploded = df_combined.explode(103, ignore_index=True)
    df_exploded = df_exploded.drop(columns=[102])
    
    df_center = pd.concat([df_exploded.iloc[:, 38:63], df_exploded.iloc[:,100:]], axis=1).groupby([100,101,103]).mean().reset_index()

    # mean of all the bins
    df_center = df_center.groupby([100, 101, 103]).apply(lambda group: group.iloc[:, 3:].mean().mean()).reset_index(name='mean') 


    long_df = pd.melt(df_center, id_vars=[100, 101, 103], var_name='column', value_name='value')
    long_df = long_df.rename(columns={100: 'set', 101: 'trt'})
    long_df['value'] = long_df.value.astype(float)

    long_df['set'] = long_df['set'].apply(lambda x: x.split('_')[0])
    long_df['set'] = pd.Categorical(long_df['set'], categories=['Controls', long_df['set'].unique()[long_df['set'].unique() != 'Controls'][0]], ordered=True)

    # skip zero
    # long_df = long_df[long_df['value'] != 0]
    long_df['logvalue'] = np.log10(long_df['value'] + 0.001)


    # ## to save the tables
    # # remove the marker for the duplicated control genes
    # df_center[103] = df_center[103].str.split('.').str[0]
    # df_center.to_csv(f'{pn}.{day}.center500bp.refgene.txt', sep='\t', index=False)



    ## log transfrom
    g = sns.catplot(data=long_df, kind="violin",palette=["#edb9aa", "#7c84f4"],
                x='set', y='logvalue', col='trt',
                saturation=0.7, linewidth=.3, inner='box',
                aspect=.8)
    # stats
    pairs = [tuple(set(long_df.set))]

    ant = Annotator(None, pairs)
    kwargs = {
        'plot_params': { # this takes what normally goes into sns.barplot etc.
            'x': 'set',
            'y': 'logvalue',
            'palette':["#edb9aa", "#7c84f4"]
        },
        'annotation_func': 'apply_test', # has three options
        'configuration': {'test': 'Mann-Whitney', 'text_format' :'full'}, # this takes what normally goes into ant.configure
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
        data = long_df[
            (long_df['trt'] == ax.get_title().split('=')[1].strip()) #&
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
            ax.plot([x_position - 0.2, x_position + 0.2], [logmean_val, logmean_val], color='#FFD700', linestyle=':', linewidth=1.5)

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

    g.fig.suptitle(f'{pn}'.split('.')[0], y=1.02, fontsize=12)

    plt.savefig(f'{pn}_tsscentral2.4kbp_k27me3_log.pdf', bbox_inches='tight')
    plt.close()
     