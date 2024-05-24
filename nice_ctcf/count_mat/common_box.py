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
    mat = pd.read_csv(f, sep=' ', header=0, index_col=0)
    mat = mat.loc[:,mat.columns.str.contains(day)]

    ctcf_mats.append(mat)

# agg nice mats
nice_mats = []
for f in nice_files:
    day = f.split('/')[-1].split('.')[5]
    mat = pd.read_csv(f, sep=' ', header=0, index_col=0)
    mat = mat.loc[:,mat.columns.str.replace('Day','D').str.contains(day)]

    nice_mats.append(mat)
