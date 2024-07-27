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


## load tpm
# filtered TPM
# tpm = pd.read_csv('../tpm_symbol_filtered_noncrna.txt', header=0, index_col=0, sep='\t')
tpm = pd.read_csv('../tpm_symbol.txt', header=0, index_col=0, sep='\t')

# only keep day 1 vihicle
tpm = tpm.loc[:,tpm.columns.str.contains('Day1_C')]


## gene list and controls
glfiles =[ file for file in  glob.glob('*.txt') if 'controls' not in file ]

for f in glfiles:
    pn = f.split('.t')[0]

    path = pd.read_csv(f, header=None, sep='\t')[0].to_list()
    ctrl = pd.read_csv(f'{pn}.controls.txt', header=None, sep='\t')[0].to_list()

    # path_set_raw = tpm.loc[tpm.index[tpm.index.isin(path)]].stack().to_frame().assign(geneset=f'{pn}'.split('_')[0])
    # ctrl_raw = tpm.loc[tpm.index[tpm.index.isin(ctrl)]].stack().to_frame().assign(geneset='ctrl')

    path_set_raw = tpm.loc[path].stack().to_frame().assign(geneset=f'{pn}'.split('_')[0])
    ctrl_raw = tpm.loc[ctrl].stack().to_frame().assign(geneset='ctrl')

    # all_sets = pd.concat([path_set_raw, ctrl_raw]).reset_index().rename(columns={'level_1':'sample', 0:'tpm'})
    # all_sets['day'] = all_sets['sample'].str.split('_').str[0]
    # all_sets['trt'] = all_sets['sample'].str.split('_').str[1].str[0]

    # match the control gene to the ref geen, avoid duplicated control genes get merged
    all_sets = pd.concat([path_set_raw.reset_index(), ctrl_raw.reset_index()], axis=1).reset_index()
    all_sets = all_sets.drop(all_sets.columns[0], axis=1)
    all_sets.columns = ['refgene', 'sample', 'tpm_ref', 'geneset', 'symbol_ctrl', 'sample_ctrl', 'tpm_ctrl', 'geneset_ctrl']
    
    # merge by refgene
    plotdat = all_sets.groupby('refgene').mean().reset_index()

    # wide to long
    plotdat = plotdat.melt(id_vars='refgene', value_vars=['tpm_ref', 'tpm_ctrl'])
    plotdat.columns = ['refgene','geneset', 'tpm']
    plotdat.geneset = plotdat.geneset.str.replace('tpm_ref', f'{pn}').str.replace('tpm_ctrl', 'ctrl')
    plotdat.geneset = pd.Categorical(plotdat.geneset, categories=[f'{pn}'.split('_')[0], 'ctrl'], ordered=True)


    ## plot
    plt.figure(figsize=(4, 5))
    g = sns.violinplot(data=plotdat, x='geneset', y='tpm',
                   palette = ["#CD9B9B", "#CDC9C9" ],
                   saturation=0.9, linewidth=1, inner='box')

    g = sns.stripplot(data=plotdat, x="geneset", y='tpm', jitter=True,
                  color='black', alpha=0.3, dodge=True, size=2)

    df_mean = plotdat.groupby('geneset', sort=False)['tpm'].mean()
    _ = [g.hlines(y, i-.12, i+.12, zorder=2, colors='#FFD700', linestyle=':', linewidth=1.5) for i, y in df_mean.reset_index()['tpm'].items()]

    nobs = ["mean=" + str(f"{i:.2f}") for i in df_mean.to_list()]

    pos = range(len(nobs))
    for tick, label in zip(pos, g.get_xticklabels()):
       g.text(pos[tick], df_mean[tick] + 1, nobs[tick],
                horizontalalignment='center',
                size='small',
                color='black')

    # plt.ylim([-4.99, None])
    plt.ylabel('TPM')
    plt.xlabel('')
    plt.title(f'{pn}', fontdict={'fontsize':7})

    # stats

    if plotdat.geneset.value_counts()[0] == plotdat.geneset.value_counts()[1]:
        statstest = 'Wilcoxon'
    else:
        statstest = 'Mann-Whitney'

    comb = [tuple(set(plotdat.geneset))]
    annotator = Annotator(g, comb, data=plotdat, x='geneset', y='tpm')
    annotator.configure(test=statstest, text_format="full",loc='inside')
    annotator.apply_and_annotate()

    plt.savefig(f'{pn}_TPM_vsctrl_d1veh_unfilteredinput.pdf', bbox_inches='tight')