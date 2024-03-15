######### nomalizing to the most ##########


import os, sys, glob
import pandas as pd
from itertools import product


## read ms table
msdat = pd.read_csv('miapaca2_ms.csv', sep=',', header=0, index_col=0)



## heatmap of intereseted marks
# marks have C&R data
marks = ['H3K27ac', 'H3K27me3', 'H3K36me1', 'H3K36me2', 'H3K36me3', 'H3K4me1',
         'H3K4me3', 'H3K9me3']



## calculate scaling facotr for cut & run
sf_df = msdat.loc[[row in marks for row in msdat.index],~msdat.columns.str.contains('Rep3')]
sf_df =  sf_df.apply(lambda x: x/x.max(), axis=1)
sf_df_long = sf_df.reset_index().melt(id_vars='index', var_name='Sample', value_name='Scaling Factor')
sf_df_long['sample'] = sf_df_long['Sample'].str.replace('DMSO','Vehicle') + \
    '_' + sf_df_long['index'].str.replace('ac','Ac')


# sf_df_long.to_csv('ms_scaling_factor.txt', index=False, sep='\t', header=False)


## calculate subsampling factor for cut & run
## unique reads number
uniq_stats = pd.read_csv('miapaca2_cr_bt2.clean.stats.txt',sep='\t',index_col=0)

subfactor_df = pd.merge(sf_df_long, uniq_stats[['uniq_hg']], left_on='sample', right_index=True, how='inner')

subfactor_df['expcted_reads'] = subfactor_df.apply(lambda row: row['uniq_hg'] if row['Scaling Factor'] == 1 \
                                                   else row['Scaling Factor'] * subfactor_df.loc[(subfactor_df['index'] == row['index']) & (subfactor_df['Scaling Factor'] == 1), 'uniq_hg'].values[0], axis=1)

subfactor_df['subsample_fraction'] = subfactor_df['expcted_reads'] / subfactor_df['uniq_hg']

subfactor_df.to_csv('test.txt', index=False, sep='\t', header=False)