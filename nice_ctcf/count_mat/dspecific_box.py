import os
import glob
import numpy as np
import pandas as pd
import seaborn as sns
from matplotlib import pyplot as plt
from scipy import stats


## list files
# ctcf read counts in specific peaks
ctcf_files = [f for f in glob.glob('count_ctcf/*.txt') if '111111' not in f]

# nice-seq read counts in specific peaks
nice_files = [f for f in glob.glob('count_nice/*.txt') if '111111' not in f]

# agg ctcf mats
ctcf_mats = []
for f in ctcf_files:
    day = f.split('/')[-1].split('.')[6]
    ctcfcon = f.split('/')[-1].split('.')[3]
    type = f.split('/')[-1].split('.')[5].replace('cloes','close')
    mat = pd.read_csv(f, sep=' ', header=0, index_col=0)
    mat = mat.loc[:,mat.columns.str.replace('Vehicle','vehicle').str.contains(f'{ctcfcon}_{day}')]
    mat.columns = [col + f'_{type}' for col in mat.columns]

    ctcf_mats.append(mat)

# agg nice mats
nice_mats = []
for f in nice_files:
    day = f.split('/')[-1].split('.')[6]
    ctcfcon = f.split('/')[-1].split('.')[3]
    type = f.split('/')[-1].split('.')[5].replace('cloes','close')
    mat = pd.read_csv(f, sep=' ', header=0, index_col=0)
    mat = mat.loc[:,mat.columns.str.replace('Day','D').str.replace('VEH','vehicle').str.contains(f'{day}_{ctcfcon}')]
    mat.columns = [col + f'_{type}' for col in mat.columns]

    nice_mats.append(mat)


# agg nice and ctcf
all_dat = []

for n, c in zip(nice_mats, ctcf_mats):
    n = n.stack().reset_index()
    n.level_1 = n.level_1.str.replace('Day','D').str.replace('-full','').str.replace('_rep','_Rep').str.replace('VEH','Vehicle')
    n.level_1 = n.level_1.str.split('_').apply(lambda x: '_'.join([x[1], x[0], x[2], 'NiCE', x[3]]))
    n['attr'] = n.level_1.str.replace('_NiCE','')

    c = c.stack().reset_index()
    c['attr'] = c.level_1.str.replace('_CTCF','')

    # merge
    merged = pd.merge(n, c, on=['Geneid', 'attr'])
    merged['day'] = merged.level_1_y.str.split('_').apply(lambda x: x[1])
    merged['trt'] = merged.level_1_y.str.split('_').apply(lambda x: x[0])
    merged['rep'] = merged.level_1_y.str.split('_').apply(lambda x: x[2])
    merged['type'] = merged.level_1_y.str.split('_').apply(lambda x: x[4])
    merged.drop(columns=['level_1_x', 'level_1_y', 'attr'], inplace=True)
    merged.columns = ['Geneid', 'nice', 'ctcf', 'day', 'trt', 'rep', 'type']

    all_dat.append(merged)

all_dat = pd.concat(all_dat)


### plot
## NSD2i/vehicle
nsd2i_data = all_dat[all_dat['trt'] == 'NSD2i']
vehicle_data = all_dat[all_dat['trt'] == 'Vehicle']

# Perform the division
merged_data = pd.merge(nsd2i_data, vehicle_data, on=['Geneid', 'day', 'rep', 'type'], suffixes=('_NSD2i', '_Vehicle'))
merged_data['nice_ratio'] = np.log2((merged_data['nice_NSD2i'] +1)/(merged_data['nice_Vehicle']+1))
merged_data['ctcf_ratio'] = np.log2((merged_data['ctcf_NSD2i'] +1)/(merged_data['ctcf_Vehicle']+1))
merged_data.drop(columns=['trt_NSD2i', 'trt_Vehicle'], inplace=True)


# Melt the data
melted_melt_nice_ratiodata = pd.melt(merged_data, id_vars=['Geneid', 'day', 'rep', 'type'], value_vars=['nice_ratio', 'ctcf_ratio'], var_name='variable', value_name='value')
melted_melt_nice_ratiodata['day'] = pd.Categorical(melted_melt_nice_ratiodata['day'], categories=['D1', 'D5', 'D9'])

# read common peaks data
common_data = pd.read_csv('common_ratio.txt', sep='\t', header=0)
plot_dat = pd.concat([melted_melt_nice_ratiodata, common_data])
plot_dat['day'] = pd.Categorical(plot_dat['day'], categories=['D1', 'D5', 'D9'])
plot_dat['type'] = pd.Categorical(plot_dat['type'], categories=['common', 'open', 'close'])



## Plot the ratio
# violin
plt.figure(figsize=(4, 5))
g = sns.catplot(data=plot_dat, kind="violin", hue="variable",
            x="day", y='value', col="type", split=True,height=4, aspect=.5)
g.set_axis_labels("","log2(NSD2i/Vehicle)")
# g.refline(y=0, color='black', linestyle='--')
for ax in g.axes.flat:
    ax.axhline(0, color='black', linestyle='--', linewidth=0.5)
# g.set_ylabel()
# g.legend_.set_title(None)
plt.savefig('ratio_dayspecific_common.pdf')
plt.close()


# scatter
pivot_for_scatter = plot_dat.pivot_table(index=['Geneid', 'day', 'rep', 'type'], columns='variable', values='value').reset_index()

plt.figure(figsize=(5, 5))
g = sns.FacetGrid(pivot_for_scatter, col="type")
g.map(sns.scatterplot, "nice_ratio", "ctcf_ratio")


def annotate(data, **kws):
    x, y = data['nice_ratio'], data['ctcf_ratio']
    r, _ = stats.pearsonr(x, y)  # Calculate the Pearson correlation coefficient
    ax = plt.gca()  # Get current axis
    ax.text(0.05, 0.95, f'r={r:.2f}, p={_:.2f}', transform=ax.transAxes, verticalalignment='top', fontweight='bold')


g = sns.FacetGrid(pivot_for_scatter, col="type", row='day')
g.map_dataframe(sns.regplot, "nice_ratio", "ctcf_ratio", scatter_kws={'alpha': 0.5, 's': 2}, line_kws={'color': None})
g.map_dataframe(annotate)
g.add_legend()
# g.set_axis_labels("NiCE-seq log2(NSD2i/Vehicle)", "CTCF log2(NSD2i/Vehicle)")  # Set labels for x and y axes
g.axes[-1, 1].set_xlabel('NiCE-seq log2(NSD2i/Vehicle)')  # Only the bottom left subplot gets the x-axis label
g.axes[1, 0].set_ylabel('CTCF log2(NSD2i/Vehicle)')  # On

plt.tight_layout()
plt.savefig('ratio_scatter.pdf')
plt.close()
