import os, glob
import pandas as pd
import numpy as np

## load bed files
files = glob.glob('*genehancer.bed')

gene_list = []
for f in files:
    bed = pd.read_csv(f, sep='\t', header=None)
    bed = pd.concat([bed[[0,1,2,5,8]], bed[11].str.split(';', expand=True)[[0,1,2]].apply(lambda x: x.str.split('=', expand=True)[1])], axis=1)
    bed.columns = ['chr', 'start', 'end', 'feature', 'score', 'genehancer_id', 'connected_gene', 'gene_score']
    bed_maxscore = bed.loc[bed.groupby(['chr', 'start', 'end'])['score'].idxmax()]
    bed_maxscore['gene_score'] = bed_maxscore['gene_score'].astype(float)

    gene_list.append(bed_maxscore[['connected_gene', 'gene_score']])

gene_score = pd.concat(gene_list)

gene_score = gene_score.groupby('connected_gene')['gene_score'].mean().sort_values(ascending=False)

gene_score.to_csv('genehancer_score.csv', header=False, sep='\t')