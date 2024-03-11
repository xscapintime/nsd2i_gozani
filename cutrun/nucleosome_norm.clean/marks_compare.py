import os 
import pandas as pd
from matplotlib import pyplot as plt
import matplotlib
import seaborn as sns

matplotlib.rc('pdf', fonttype=42)


## load
barcode_stats = pd.read_csv('barcode_stats_df.txt', sep='\t', header=0, index_col=0)

## intertesed marks
marks = ['H3K27ac', 'H3K27me3', 'H3K36me2', 'H3K36me3']

## subset
dats = barcode_stats.loc[barcode_stats.index.str.contains('H3K27ac|H3K27me3|H3K36me2|H3K36me3')]
dats = dats.loc[:,dats.columns.str.contains('uniq_hg|tot_barcode|nuc_norm_hg|target_barcode|target_norm_hg|H3K27ac|H3K27me3|H3K36me2|H3K36me3')]

dats = dats[dats.columns.drop('H3K27ac+S28phos')]


## dot plot
for mark in marks:
    dat = dats[dats.index.str.contains(mark)]

    dat['Trt'] = dat.index.str.split('_').str[0]
    dat['Day'] = dat.index.str.split('_').str[1]

    x=dat.target_barcode/dat.tot_barcode
    y=dat.nuc_norm_hg

    plt.figure(figsize=(5,4.6))
    ax = sns.scatterplot(data=dat,\
        x=x,
        y=y,
        hue=dat['Day'],
        style=dat['Trt'])

    # for i, txt in enumerate(dat.index.to_list()):
    #     ax.text(x[i]*0.8, y[i]-2, txt)
    plt.title(mark)
    plt.xlabel('Target/all barcodes')
    plt.ylabel('Human reads/all barcodes')
    
    plt.tight_layout()
    plt.savefig(f'{mark}_allbarcodes_dotplot.pdf')
    plt.close()



## dot plot for target barcoes norm
for mark in marks:
    dat = dats[dats.index.str.contains(mark)]

    dat['Trt'] = dat.index.str.split('_').str[0]
    dat['Day'] = dat.index.str.split('_').str[1]

    x=dat.target_barcode/dat.tot_barcode
    y=dat.target_norm_hg

    plt.figure(figsize=(5,4.6))
    ax = sns.scatterplot(data=dat,\
        x=x,
        y=y,
        hue=dat['Day'],
        style=dat['Trt'])

    # for i, txt in enumerate(dat.index.to_list()):
    #     ax.text(x[i]*0.8, y[i]-2, txt)
    plt.title(mark)
    plt.xlabel('Target/all barcodes')
    plt.ylabel('Human reads/target barcode')

    plt.tight_layout()
    plt.savefig(f'{mark}_tarbarcodes_dotplot.pdf')
    plt.close()
