######### nomalizing to the least ##########
### merged rep1-rep3 of each mass spec data ####

import os, sys, glob
import pandas as pd



## read ms table
msdat = pd.read_csv('../miapaca2_ms.csv', sep=',', header=0, index_col=0)

## wide to long
ms_long = msdat.reset_index().melt(id_vars='index', var_name='Sample', value_name='ms')

## marks have C&R data
marks = ['H3K27ac', 'H3K27me3', 'H3K36me1', 'H3K36me2', 'H3K36me3', 'H3K4me1',
         'H3K4me3', 'H3K9me3']


ms_long = ms_long[ms_long['index'].isin(marks)]


## merged rep1-rep3
ms_long['condition'] = ms_long['Sample'].str.split('_').str[0]
ms_long['day'] = ms_long['Sample'].str.split('_').str[1]
ms_long['rep'] = ms_long['Sample'].str.split('_').str[2]

ms_long = ms_long.groupby(['index', 'condition', 'day']).agg({'ms':'mean'}).reset_index()


## smallest as 1
ms_long['scale'] = ms_long.groupby('index')['ms'].transform(lambda x: x / x.min())

ms_long['condition'] = ms_long['condition'].str.replace('DMSO','Vehicle')
ms_long['index'] = ms_long['index'].str.replace('ac','Ac')


ms_long.to_csv('ms_merged_scaling_factor.txt', index=False, sep='\t', header=False)

