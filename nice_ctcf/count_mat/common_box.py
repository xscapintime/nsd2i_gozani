import os
import glob
import numpy as np
import pandas as pd

## list files
# ctcf read counts in commmon peaks
ctcf_files = glob.glob('count_ctcf/*.111111.txt')

# nice-seq read counts in common peaks
nice_files = glob.glob('count_nice/*.111111.txt')

# agg ctcf mats
ctcf_mats = []
for f in ctcf_files:
    day = f.split('/')[-1].split('.')[5]
    ctcfcon = f.split('/')[-1].split('.')[3]
    mat = pd.read_csv(f, sep=' ', header=0, index_col=0)
    mat = mat.loc[:,mat.columns.str.contains(f'{ctcfcon}_{day}')]

    ctcf_mats.append(mat)

# agg nice mats
nice_mats = []
for f in nice_files:
    day = f.split('/')[-1].split('.')[5]
    ctcfcon = f.split('/')[-1].split('.')[3]
    mat = pd.read_csv(f, sep=' ', header=0, index_col=0)
    mat = mat.loc[:,mat.columns.str.replace('Day','D').str.replace('VEH','Vehicle').str.contains(f'{day}_{ctcfcon}')]

    nice_mats.append(mat)


# agg nice and ctcf
all_dat = []

for n, c in zip(nice_mats, ctcf_mats):
    n = n.stack().reset_index()
    n.level_1 = n.level_1.str.replace('Day','D').str.replace('-full','').str.replace('_rep','_Rep').str.replace('VEH','Vehicle')
    n.level_1 = n.level_1.str.split('_').apply(lambda x: '_'.join([x[1], x[0], x[2], 'NiCE']))
    n['attr'] = n.level_1.str.replace('_NiCE','')

    c = c.stack().reset_index()
    c['attr'] = c.level_1.str.replace('_CTCF','')

    # merge
    merged = pd.merge(n, c, on=['Geneid', 'attr'])
    merged['day'] = merged.level_1_y.str.split('_').apply(lambda x: x[1])
    merged['trt'] = merged.level_1_y.str.split('_').apply(lambda x: x[0])
    merged['rep'] = merged.level_1_y.str.split('_').apply(lambda x: x[2])
    merged.drop(columns=['level_1_x', 'level_1_y', 'attr'], inplace=True)
    merged.columns = ['Geneid', 'nice', 'ctcf', 'day', 'trt', 'rep']

    all_dat.append(merged)

all_dat = pd.concat(all_dat)