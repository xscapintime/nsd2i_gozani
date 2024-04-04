import os,glob
import pandas as pd
import numpy as np

files = glob.glob('*superset.bed')

samplelist = ['NSD2i_D1', 'NSD2i_D5', 'NSD2i_D9',
              'Vehicle_D1', 'Vehicle_D5', 'Vehicle_D9']



for f in files:
    sam = f.split('.')[0]

    set = pd.read_csv(f, header=None, sep='\t')
    set[3] = set[3].str.split(',').apply(lambda x: sorted(np.unique(x)))
    set['binary_code'] = set[3].apply(lambda x : ''.join(['1' if item in x else '0' for item in samplelist]))
    set[3] = set[3].apply(lambda x : ','.join(x))
    set.to_csv(f'{sam}.binarycode.bed', sep='\t', header=False, index=False)