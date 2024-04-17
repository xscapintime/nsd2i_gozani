import os,glob
import pandas as pd
import numpy as np
import seaborn as sns
import matplotlib.pyplot as plt

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

    code_order = ['111111',
                  '100000', '110000', '111000', '110100', '111110', '111100',
                  '000111', '100111', '110111', '000110', '000100', '100110',
                  '010000', '011000', '001000',
                  '000010', '000011', '000001',
                  '100100', '110110', '101101', '001001', '011011', '010010']
    plot_dat = set[set.binary_code.isin(code_order)]
    plot_dat['group'] = np.where(plot_dat.binary_code.isin(['100000', '110000', '111000', '110100', '111110', '111100']), 'gain',
                                 np.where(plot_dat.binary_code.isin(['000111', '100111', '110111', '000110', '000100', '100110']), 'loss',
                                 np.where(plot_dat.binary_code.isin(['010000', '011000', '001000', '000010', '000011', '000001']), 'unique',
                                 np.where(plot_dat.binary_code.isin(['100100', '110110', '101101', '001001', '011011', '010010']), 'no change', 'shared'))))

    plt.figure(figsize=(5,6))
    sns.countplot(data=plot_dat, y='binary_code', order=code_order,
                    hue=plot_dat['group'], dodge=False,
                    hue_order=['shared', 'gain', 'loss', 'unique', 'no change'])
    plt.legend(bbox_to_anchor=(1.02, 1), loc='upper left', borderaxespad=0)
    plt.tight_layout()
    plt.savefig(f'{sam}.binarycode.pdf')

    plot_dat['peaklen'] = plot_dat[2] - plot_dat[1]
    sns.catplot(data=plot_dat,  y='peaklen', kind='violin', col='group',
    col_order=['shared', 'gain', 'loss', 'unique', 'no change'], col_wrap=3)
    plt.tight_layout()
    plt.savefig(f'{sam}.peaklength.pdf')

