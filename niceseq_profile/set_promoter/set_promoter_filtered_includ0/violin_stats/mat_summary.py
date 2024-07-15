import os,glob, gzip
import json
import pandas as pd
import numpy as np
import seaborn as sns
import matplotlib
from matplotlib import pyplot as plt
import matplotlib.font_manager
from itertools import product


plt.style.use('seaborn-poster')
matplotlib.rcParams['pdf.fonttype'] = 42
matplotlib.rcParams['ps.fonttype'] = 42


def split_and_assign(df, row_split, col_split, row_group, col_group):
    split_df = df.iloc[row_split, col_split].copy()
    split_df['row_group'] = row_group
    split_df['col_group'] = col_group
    return split_df


for m in glob.glob('../*.promoter.3k.gz'):
    pn = os.basename(m).replace('.promoter.3k.gz','')

    ## load deeptools matrix
    df = pd.read_csv(m, header=None, skiprows=1,sep='\t')

    ## meta data
    with gzip.open(m, 'rb') as f:
        meta = f.readline()

    group_boundaries = int(meta.decode("utf-8").split('":')[20].split(',"')[0].split(',')[1])
    sample_boundaries = int(meta.decode("utf-8").split('":')[22].split(',')[1])


    # Define splits
    row_splits = [slice(None, group_boundaries), slice(group_boundaries, None)]
    col_splits = [slice(None, sample_boundaries), slice(sample_boundaries, None)]
    row_groups = [f'{pn}', 'Controls']
    col_groups = ['NSD2i.day_1', 'vehicle.day_1']

    split_dfs = [
        split_and_assign(df, row_split, col_split, row_group, col_group)
        for row_split, row_group in zip(row_splits, row_groups)
        for col_split, col_group in zip(col_splits, col_groups)
    ]

    # Concatenate all parts
    df_combined = pd.concat(split_dfs, axis=0)

    # Add row identifier before melting
    df_combined['row'] = df_combined.index

    # Transform to long format
    long_df = pd.melt(df_combined, id_vars=['row', 'row_group', 'col_group'], var_name='column', value_name='value')
