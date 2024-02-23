import os
import pandas as pd


union = pd.read_csv('union.bed', header=None, sep='\t')

samplelist = ['NSD2i_Day_1', 'NSD2i_Day_5', 'NSD2i_Day_9',
              'Vehicle_Day_1', 'Vehicle_Day_5', 'Vehicle_Day_9']

## name sorting and cleaning
union[3] = union[3].apply(lambda x : ','.join(list(sorted(set(x.split(','))))))

intersect_all = union[union[3] == ','.join(samplelist)]
intersect_nsd2i = union[union[3] == ','.join(samplelist[:3])]
intersect_vehicle = union[union[3] == ','.join(samplelist[3:6])]

## intersection
for df in [ intersect_all, intersect_nsd2i, intersect_vehicle ]:
    fn = [name for name, obj in globals().items() if obj is df][0]
    df.to_csv(f'{fn}.bed', sep='\t', header=False, index=False)


## treatment/time unique
for sam in samplelist:
    uniq = union[union[3] == sam]
    
    uniq.to_csv(f'{sam}.unique.bed', sep='\t', header=False, index=False)