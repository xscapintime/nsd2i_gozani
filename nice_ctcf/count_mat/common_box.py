import os
import glob
import numpy as np
import pandas as pd
import seaborn as sns
from matplotlib import pyplot as plt


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
    mat = mat.loc[:,mat.columns.str.replace('Vehicle','vehicle').str.contains(f'{ctcfcon}_{day}')]

    ctcf_mats.append(mat)

# agg nice mats
nice_mats = []
for f in nice_files:
    day = f.split('/')[-1].split('.')[5]
    ctcfcon = f.split('/')[-1].split('.')[3]
    mat = pd.read_csv(f, sep=' ', header=0, index_col=0)
    mat = mat.loc[:,mat.columns.str.replace('Day','D').str.replace('VEH','vehicle').str.contains(f'{day}_{ctcfcon}')]

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


### plot
melt_dat = pd.melt(all_dat, id_vars=['Geneid', 'day', 'trt', 'rep'], value_vars=['nice', 'ctcf'], var_name='variable', value_name='value')
melt_dat['log2count'] = melt_dat['value'].apply(lambda x: np.log2(x+1))
melt_dat['day'] = pd.Categorical(melt_dat['day'], categories=['D1', 'D5', 'D9'])

## NSD2i and vehicle sepearate
plt.figure(figsize=(4,5))
g = sns.catplot(
    data=melt_dat, kind="violin", hue="variable",
    x="day", y='log2count', col="trt",
    height=4, aspect=.5,
)
g.set_axis_labels("", "log2(read count+1)")
# plt.tight_layout()
plt.savefig('counts_commonpks.pdf')
plt.close()



## NSD2i/vehicle
nsd2i_data = all_dat[all_dat['trt'] == 'NSD2i']
vehicle_data = all_dat[all_dat['trt'] == 'Vehicle']

# Perform the division
merged_data = pd.merge(nsd2i_data, vehicle_data, on=['Geneid', 'day', 'rep'], suffixes=('_NSD2i', '_Vehicle'))
merged_data['nice_ratio'] = np.log2((merged_data['nice_NSD2i'] +1)/(merged_data['nice_Vehicle']+1))
merged_data['ctcf_ratio'] = np.log2((merged_data['ctcf_NSD2i'] +1)/(merged_data['ctcf_Vehicle']+1))
merged_data.drop(columns=['trt_NSD2i', 'trt_Vehicle'], inplace=True)


# Melt the data
melted_melt_nice_ratiodata = pd.melt(merged_data, id_vars=['Geneid', 'day', 'rep'], value_vars=['nice_ratio', 'ctcf_ratio'], var_name='variable', value_name='value')
melted_melt_nice_ratiodata['day'] = pd.Categorical(melted_melt_nice_ratiodata['day'], categories=['D1', 'D5', 'D9'])

# export the table
melted_melt_nice_ratiodata['type'] = 'common'
melted_melt_nice_ratiodata.to_csv('common_ratio.txt', index=False, sep='\t')


# Plot the ratio
plt.figure(figsize=(4, 5))
g = sns.violinplot(data=melted_melt_nice_ratiodata, x="day", y="value",
                   hue="variable", split=True, inner="quart")
g.set_xlabel("")
g.set_ylabel("log2(NSD2i/Vehicle)")
g.legend_.set_title(None)
plt.savefig('ratio_commonpks.pdf')
plt.close()



