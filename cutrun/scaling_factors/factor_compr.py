import os, glob
import pandas as pd

## mass spec scaling factors
ms = pd.read_csv('../../masspec/ms_scaling_factor.txt', sep='\t', header=None)

## all-nucleosome scaling factors
nuc = pd.read_csv('../nucleosome_norm.clean/nuc_norm_hg.txt', sep='\t', header=None)

## target-nucleosome scaling factors
tarnuc = pd.read_csv('../nucleosome_norm.clean/tarnuc_norm_hg.txt', sep='\t', header=None)
tarnuc[0] = tarnuc[0].str.replace('ac', 'Ac')


## e.coli scaling factors
uniq_reads = pd.read_csv('../nucleosome_norm.clean/miapaca2_cr_bt2.clean.stats.txt', sep='\t', header=0)
uniq_reads['ec_norm'] = uniq_reads['uniq_hg'] / uniq_reads['uniq_ec']


## merge all scaling factors
all_factors = nuc.merge(tarnuc, left_on=0, right_on=0, how='left').\
    merge(uniq_reads[['samples', 'ec_norm']], left_on=0, right_on='samples', how='left').\
    merge(ms[[3,2]], left_on=0, right_on=3, how='left')


all_factors = all_factors[['samples',2,'1_x','1_y','ec_norm']]

all_factors.columns = ['samples', 'ms_norm', 'nuc_norm', 'tarnuc_norm', 'ec_norm']

## export
all_factors.to_csv('all_scaling_factors.txt', sep='\t', index=False, header=True)