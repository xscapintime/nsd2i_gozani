import os,glob
import pandas as pd
import numpy as np
from itertools import product
import seaborn as sns
from matplotlib import pyplot as plt



## CTCF IDR peak number
ctcf_idr = pd.read_csv('ctcf_idr/ctcf_idr.stats',sep='\t',header=None)
ctcf_idr['trtday'] = ctcf_idr[1].apply(lambda x: x.replace('.CTCF','').split('.idr')[0])


## NiCE-seq CTCF intersection
nice_ctcf = pd.read_csv('nice_ctcf/nice_ctcf.stats',sep='\t',header=None)

def assign_value(row, combinations):
    for combo in combinations:
        if all(item in row for item in combo):
            return f"{combo[0]}.{combo[1]}"
    return 'No match'

combs = list(product(['NSD2i', 'vehicle'], ['D1','D5','D9']))

nice_ctcf['trtday'] = nice_ctcf[1].apply(lambda x: assign_value(x, combs))


## merge
nice_ctcf = nice_ctcf.merge(ctcf_idr, on='trtday', how='outer')
nice_ctcf['prop'] = nice_ctcf['0_x']/nice_ctcf['0_y']*100
nice_ctcf['day'] = nice_ctcf['trtday'].apply(lambda x: x.split('.')[1])
nice_ctcf['trt'] = nice_ctcf['trtday'].apply(lambda x: x.split('.')[0])

nice_ctcf['peak_type'] = np.where(nice_ctcf['1_x'].str.contains('111111'), 'common',
                                  np.where(nice_ctcf['1_x'].str.contains('open'), 'open', 'closed'))


## plot
plt.figure(figsize=(4,5))
g = sns.catplot(
    data=nice_ctcf, kind="bar", hue="peak_type",
    x="day", y="prop", col="trt",
    height=4, aspect=.5,
)
g.set_axis_labels("", "%CTCF peaks")
# plt.tight_layout()
plt.savefig('nice_in_ctcf.pdf')
plt.close()