import os 
import pandas as pd
from matplotlib import pyplot as plt
import matplotlib
import seaborn as sns

matplotlib.rc('pdf', fonttype=42)

## read barcode stats
# Read the kmet
kmet = pd.read_csv('kmet_nuccount.txt', sep='\t', index_col=0)

# Read the kac file
kac = pd.read_csv('kac_nuccount.txt', sep='\t', index_col=0)

# Merge 
nuc_stats = pd.concat([kmet,kac]).dropna(axis=1)

# total barcode number
tot_barcode = pd.DataFrame(nuc_stats.sum())
tot_barcode['samples'] = tot_barcode.index.str.split('.').str[0]
tot_barcode = tot_barcode.groupby('samples').sum()

## unique reads number
uniq_stats = pd.read_csv('miapaca2_cr_bt2.stats.txt',sep='\t',index_col=0)


## scaling factor via nucloesome barcodes
# normalize by unique reads
# barcodes scaling factor
nuc_norm_hg = tot_barcode.apply(lambda x : x / uniq_stats['uniq_hg'][x.index])
nuc_norm_all = tot_barcode.apply(lambda x : x / uniq_stats['uniq_all'][x.index])

# export
nuc_norm_hg.to_csv('nuc_norm_hg.txt', sep='\t', header=False)
nuc_norm_all.to_csv('nuc_norm_all.txt', sep='\t', header=False)


## on-targt normalization
sum_ab = pd.DataFrame()

for i in range(0, len(nuc_stats.columns), 2):
    # Check if there's a next column to avoid IndexError
    if i+1 < len(nuc_stats.columns):
        nuc_stats.columns = nuc_stats.columns.str.replace('Ac','ac')
        nuc_stats.index = nuc_stats.index.str.replace('-','_')
        
        col_name =  nuc_stats.columns[i].split('.')[0]
        sum_ab[col_name] = nuc_stats.iloc[:, i] + nuc_stats.iloc[:, i+1]        
        
# Sum up every two rows
sum_ab = sum_ab.groupby(sum_ab.index.str.replace('_A$','').str.replace('_B$','')).sum()


def on_target_ratio(df):
    new_df = df.copy()  # Create a copy of the input dataframe
    for column in new_df.columns:
        parts = column.split('_')
        if len(parts) >= 4 and parts[3] == 'IgG':
            column_sum = new_df[column].sum()
            new_df[column] = new_df[column] / column_sum
        else:
            targetsum = new_df[new_df.index.str.split('_').str[0] == parts[3]][column].sum()
            new_df[column] = new_df[column] / targetsum
    return new_df

ontar_norm = on_target_ratio(sum_ab)


# export
ontar_norm.to_csv('ontar_norm.txt', sep='\t', header=True)