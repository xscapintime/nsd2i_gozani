import os,glob
import pandas as pd
import numpy as np
import seaborn as sns
import matplotlib
from matplotlib import pyplot as plt
import matplotlib.font_manager
from statannotations.Annotator import Annotator
from itertools import product


plt.style.use('seaborn-poster')
matplotlib.rcParams['pdf.fonttype'] = 42
matplotlib.rcParams['ps.fonttype'] = 42


## load enhancer annotaion to get the genes that eventually will be using
#ct_enhancer_normed = pd.read_csv('../../mark_profile/signal_on_enhancer/ct_enhancer_normed.txt', header=0, index_col=0, sep='\t')
#ct_enhancer_normed.columns = ct_enhancer_normed.columns.str.replace('NSD2i_', '')

gennehancer = pd.read_csv('GeneHancer_v5.bestgene.bed', header=None, sep='\t')
gennehancer.set_index(gennehancer[4], inplace=True)


## load ctrl gene tables
ctrl_files = glob.glob('../ctrlsets/*_ctrl_genes.txt')

pathname = [ os.path.basename(f).replace('_ctrl_genes.txt', '') for f in ctrl_files ]


for pn in pathname:

    ctrl_genes = pd.read_csv(f'../ctrlsets/{pn}_ctrl_genes.txt', sep='\t', header=0)
    path = list(ctrl_genes[f'{pn}'])

    pathset_df = gennehancer.loc[gennehancer.index.isin(path)]

    all_ctrls = []
    for i in range(0, 5):
        ctrl_set = gennehancer.loc[gennehancer.index.isin(ctrl_genes.iloc[:,i])]
        all_ctrls.append(ctrl_set)

    all_ctrls_df = pd.concat(all_ctrls)

    # export
    pathset_df.to_csv(f'{pn}_enhancer.bed', sep='\t', header=False, index=False)
    all_ctrls_df.to_csv(f'{pn}_allctrls_enhancer.bed', sep='\t', header=False, index=False)    

