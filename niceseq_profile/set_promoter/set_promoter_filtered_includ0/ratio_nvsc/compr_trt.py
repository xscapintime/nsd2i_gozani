import os,glob, gzip
import pandas as pd
import numpy as np
import seaborn as sns
import matplotlib
from matplotlib import pyplot as plt
import matplotlib.font_manager
from itertools import product
from statannotations.Annotator import Annotator
from matplotlib import ticker as mticker


plt.style.use('seaborn-poster')
matplotlib.rcParams['pdf.fonttype'] = 42
matplotlib.rcParams['ps.fonttype'] = 42


def split_and_assign(df, row_split, col_split, row_group, col_group):
    split_df = df.iloc[row_split, col_split].copy()
    split_df['row_group'] = row_group
    split_df['col_group'] = col_group
    return split_df

paths = list(set(sorted([ os.path.basename(m).split('.')[0] for m in glob.glob('../*.promoter.3k.gz')])))

for p in paths:
    mats = sorted(glob.glob(f'../{p}.d*.promoter.3k.gz'))

    mats_list = []
    for m in mats:
        pn = os.path.basename(m).replace('.promoter.3k.gz','')
        day = pn.split('.')[1]

        ## load deeptools matrix
        df = pd.read_csv(m, header=None, skiprows=1,sep='\t')
        # df['gene'] = df[3].str.split('.').str[0]
        # df.groupby('gene').mean()

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

        # Concatenate all parts
        df_combined = pd.DataFrame(np.vstack([pd.concat([split_dfs[0], split_dfs[2]], axis=0),
                                pd.concat([split_dfs[1], split_dfs[3]], axis=0).values]))

        # average all the genes in the same set/trt at center
        df_center = pd.concat([df_combined.iloc[:, 55:65].mean(axis=1), df_combined.iloc[:,120:]], axis=1)#.groupby([120,121])#.mean()

        df_cetersub = pd.DataFrame()

        df_cetersub[f"{p}".split('_')[0]] = df_center[(df_center[120] == f'{pn}') & (df_center[121] == f'NSD2i.day_{day}')].reset_index()[0] \
                    / df_center[(df_center[120] == f'{pn}') & (df_center[121] == f'vehicle.day_{day}')][0].to_list()

        df_cetersub['Controls'] = df_center[(df_center[120] == 'Controls') & (df_center[121] == f'NSD2i.day_{day}')].reset_index()[0] \
                    / df_center[(df_center[120] == 'Controls') & (df_center[121] == f'vehicle.day_{day}')][0].to_list()        
        

        long_df = pd.melt(df_cetersub)
        long_df['day'] = day

        mats_list.append(long_df)
    
    mat_df_long = pd.concat(mats_list)

    #remove NA
    mat_df_long = mat_df_long.dropna()
    # remove inf
    mat_df_long = mat_df_long[~mat_df_long.value.isin([np.inf, -np.inf])]
    # remove 0
    mat_df_long = mat_df_long[mat_df_long.value != 0]


    # log
    mat_df_long['logvalue'] = np.log10(mat_df_long['value'])




    # ## plot
    # g = sns.catplot(data=mat_df_long, kind="violin",palette=["#7CCD7C", "#CDC9C9"],
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
    #         'palette':["#7CCD7C", "#CDC9C9"]
    #     },
    #     'annotation_func': 'apply_test', # has three options
    #     'configuration': {'test': 'Mann-Whitney', 'text_format' :'full'}, # this takes what normally goes into ant.configure
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
    #     for category in data['variable'].unique():
    #         subset = data[data['variable'] == category]
    #         mean_val = subset['value'].mean()

    #         # Get the position of the category on the x-axis
    #         x_position = data['variable'].unique().tolist().index(category)

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

    # g.set_ylabels('NicE-seq (NSD2i / Vehicle)')

    # g.fig.suptitle(f'{p}', y=1.02, fontsize=12)

    # plt.savefig(f'{p}_ratio_tss_central500bp_nice.pdf', bbox_inches='tight')
    # plt.close()



    ## log transfrom
    g = sns.catplot(data=mat_df_long, kind="violin",palette=["#7CCD7C", "#CDC9C9"],
                x='variable', y='logvalue', col='day',
                saturation=0.7, linewidth=.3, inner='box',
                aspect=.8)
    # stats
    pairs = [tuple(set(mat_df_long.variable))]

    ant = Annotator(None, pairs)
    kwargs = {
        'plot_params': { # this takes what normally goes into sns.barplot etc.
            'x': 'variable',
            'y': 'logvalue',
            'palette':["#7CCD7C", "#CDC9C9"]
        },
        'annotation_func': 'apply_test', # has three options
        'configuration': {'test': 'Mann-Whitney', 'text_format' :'full'}, # this takes what normally goes into ant.configure
        'plot': 'violinplot'
    }

    g.map_dataframe(ant.plot_and_annotate_facets, **kwargs)

    for ax in g.axes.flat:
        ax.yaxis.set_major_formatter(mticker.StrMethodFormatter("$10^{{{x:.0f}}}$"))
        ymin, ymax = ax.get_ylim()
        tick_range = np.arange(np.floor(ymin), ymax)
        tick_range = np.arange(np.floor(ymin), ymax)
        ax.yaxis.set_ticks(tick_range)
        ax.yaxis.set_ticks([np.log10(x) for p in tick_range for x in np.linspace(10 ** p, 10 ** (p + 1), 10)], minor=True)



    # add mean
    # Iterate over each subplot
    for ax in g.axes.flat:
        # Get the data for this subplot
        data = mat_df_long[
            (mat_df_long['day'] == ax.get_title().split('=')[1].strip()) #&
            # (plot_df['trt'] == ax.get_title().split('|')[0].split('=')[1].strip())
        ]

        # Iterate over each category in the x-axis
        for category in data['variable'].unique():
            subset = data[data['variable'] == category]
            mean_val = subset['value'].mean()
            logmean_val = subset['logvalue'].mean()

            # Get the position of the category on the x-axis
            x_position = data['variable'].unique().tolist().index(category)

            # Draw a short horizontal line at the mean value
            ax.plot([x_position - 0.2, x_position + 0.2], [logmean_val, logmean_val], color='#FFD700', linestyle=':', linewidth=1.5)

            # Annotate the mean value
            ax.text(
                x=x_position,
                y=logmean_val,
                s=f'mean={mean_val:.3f}',
                color='black',
                ha='center',
                va='bottom',
                fontsize='small'
            )


    g.set_ylabels('NicE-seq (NSD2i / Vehicle)')

    g.fig.suptitle(f'{p}', y=1.02, fontsize=12)

    plt.savefig(f'{p}_ratio_tss_central500bp_nice_log.pdf', bbox_inches='tight')
    plt.close()
    

