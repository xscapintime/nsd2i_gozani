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
## load ensemble id to symbol matching table
ensembl_syb = pd.read_csv('../../rnaseq/umap/human_ensembl_syb.tsv', header=0, index_col=None, sep='\t')
ensembl_syb = dict(zip(ensembl_syb['ensembl_gene_id'], ensembl_syb['hgnc_symbol']))

## load tpm table
tpm = pd.read_csv('../../rnaseq/tpm/tximport-tpm.csv', header=0, index_col=0)

## change index to gene symbol
tpm = tpm.set_index(tpm.index.map(ensembl_syb))
tpm = tpm[tpm.index.notnull()]

# merge transcripts, only ~100 genes with more than 1 transcript
tpm = tpm.groupby(tpm.index).mean()

## log2tpm
tpm_normed = np.log2(np.divide(*(tpm.iloc[:,np.where(tpm.columns.str.contains('_N'))[0]]+1).\
                    align((tpm.iloc[:,np.where(tpm.columns.str.contains('_C'))[0]]+1), axis=0)))

tpm_normed.columns = tpm_normed.columns.str.replace('N', 'Rep').str.replace('Day', 'D') + '_log2tpm'


## load enhancer annotaion to get the genes that eventually will be using
#ct_enhancer_normed = pd.read_csv('../../mark_profile/signal_on_enhancer/ct_enhancer_normed.txt', header=0, index_col=0, sep='\t')
#ct_enhancer_normed.columns = ct_enhancer_normed.columns.str.replace('NSD2i_', '')

gennehancer = pd.read_csv('GeneHancer_v5.bestgene.bed', header=None, sep='\t')



## load ctrl gene tables
ctrl_files = glob.glob('../ctrlsets/*_ctrl_genes.txt')

pathname = [ os.path.basename(f).replace('_ctrl_genes.txt', '') for f in ctrl_files ]


for pn in pathname:

    ctrl_genes = pd.read_csv(f'../ctrlsets/{pn}_ctrl_genes.txt', sep='\t', header=0)
    path = list(ctrl_genes[f'{pn}'])

