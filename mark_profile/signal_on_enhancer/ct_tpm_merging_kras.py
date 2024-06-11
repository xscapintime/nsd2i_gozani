import os,glob
import pandas as pd
import numpy as np
import seaborn as sns
import matplotlib
from matplotlib import pyplot as plt
import matplotlib.font_manager

plt.style.use('seaborn-poster')
matplotlib.rcParams['pdf.fonttype'] = 42
matplotlib.rcParams['ps.fonttype'] = 42


## load log2fc cut&run
ct_enhancer_normed = pd.read_csv('ct_enhancer_normed.txt', header=0, index_col=0, sep='\t')
ct_enhancer_normed.columns = ct_enhancer_normed.columns.str.replace('NSD2i_', '')


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


# merge
# CT value on enhancers are multiple, so TPM would be dupicated
tpm_enhancersig = tpm_normed.merge(ct_enhancer_normed, left_index=True, right_index=True, how='inner')


## rank genes by k36me2 loss and low expression

tpm_enhancersig['avetpm'] = tpm_enhancersig.loc[:, tpm_enhancersig.columns.str.contains('log2tpm')].\
            apply(lambda x: np.mean(x), axis=1)
tpm_enhancersig['avek36me2'] = tpm_enhancersig.loc[:, tpm_enhancersig.columns.str.contains('K36me2')].\
            apply(lambda x: np.mean(x), axis=1)

tpm_enhancersig.sort_values(by=['avetpm', 'avek36me2',], ascending=[True, True], inplace=True)

tpm_enhancersig.drop(['avetpm', 'avek36me2'], axis=1, inplace=True)



## wide to long for plotting
# top 100 genes
plot_dat = tpm_enhancersig[:300].stack().reset_index()
plot_dat.columns = ['connected_gene', 'mark', 'value']
plot_dat['day'] = plot_dat['mark'].str.split('_').str[0]
plot_dat['rep'] = plot_dat['mark'].str.split('_').str[1]
plot_dat['mark'] = plot_dat['mark'].str.split('_').str[2]


# split tpm as another variable
tpmcol  = plot_dat[plot_dat['mark'] == 'log2tpm']
plot_dat = plot_dat[plot_dat['mark'] != 'log2tpm']
plot_dat = pd.merge(plot_dat, tpmcol, on=['connected_gene', 'day', 'rep'], how='left')


## plot top 20 genes, box plot of enhancer signal and tpm
# unique genes while keeping the order
seen = set()
top20genes = [x for x in tpm_enhancersig.index if x not in seen and not seen.add(x)][:20]

g = sns.catplot(data=plot_dat[plot_dat['connected_gene'].isin(top20genes)], kind="box",palette="pastel",
            x="connected_gene", y="value_x", row="mark_x", col='day', aspect=1.4, 
            showfliers=False, linewidth=1.2)

#https://stackoverflow.com/questions/67309730/how-to-overlay-a-scatterplot-on-top-of-boxplot-with-sns-catplot
g.map_dataframe(sns.stripplot, x="connected_gene", y="value_x", 
                palette=["#404040"],
                alpha=0.6,dodge=True,
                size=3)

g.map_dataframe(sns.stripplot, x="connected_gene", y="value_y", 
                palette=["#CD5555"],
                alpha=0.5,dodge=True,
                s=5, marker="D",
                jitter=False)

for ax in g.axes.flat:
    ax.axhline(0, color='black', linestyle='--', linewidth=0.5)
    
    for label in ax.get_xticklabels():
        label.set_rotation(90)
        label.set_horizontalalignment('center')

g.axes[1, 0].set_ylabel('log2(NSD2i/Vehicle)')

plt.savefig('lowtpm_lowk36me2.pdf', bbox_inches='tight')
plt.close()



## plot tpm vs ct signal of kras pathway genes
## kras gene
kras_up = ['ABCB1','ACE','ADAM17','ADAM8','ADAMDEC1','ADGRA2','ADGRL4','AKAP12','AKT2','ALDH1A2','ALDH1A3','AMMECR1','ANGPTL4',
 'ANKH','ANO1','ANXA10','APOD','ARG1','ATG10','AVL9','BIRC3','BMP2','BPGM','BTBD3','BTC','C3AR1','CA2','CAB39L','CBL',
 'CBR4','CBX8','CCL20','CCND2','CCSER2','CD37','CDADC1','CFB','CFH','CFHR2','CIDEA','CLEC4A','CMKLR1','CPE','CROT','CSF2',
 'CSF2RA','CTSS','CXCL10','CXCR4','DCBLD2','DNMBP','DOCK2','DUSP6','EMP1','ENG','EPB41L3','EPHB2','EREG','ERO1A','ETS1',
 'ETV1','ETV4','ETV5','EVI5','F13A1','F2RL1','FBXO4','FCER1G','FGF9','FLT4','FUCA1','G0S2','GABRA3','GADD45G','GALNT3',
 'GFPT2','GLRX','GNG11','GPNMB','GPRC5B','GUCY1A1','GYPC','H2BC3','HBEGF','HDAC9','HKDC1','HOXD11','HSD11B1','ID2','IGF2',
 'IGFBP3','IKZF1','IL10RA','IL1B','IL1RL2','IL2RG','IL33','IL7R','INHBA','IRF8','ITGA2','ITGB2','ITGBL1','JUP','KCNN4',
 'KIF5C','KLF4','LAPTM5','LAT2','LCP1','LIF','LY96','MAFB','MALL','MAP3K1','MAP4K1','MAP7','MMD','MMP10','MMP11','MMP9',
 'MPZL2','MTMR10','MYCN','NAP1L2','NGF','NIN','NR0B2','NR1H4','NRP1','PCP4','PCSK1N','PDCD1LG2','PECAM1','PEG3','PIGR',
 'PLAT','PLAU','PLAUR','PLEK2','PLVAP','PPBP','PPP1R15A','PRDM1','PRELID3B','PRKG2','PRRX1','PSMB8','PTBP2','PTCD2','PTGS2',
 'PTPRR','RABGAP1L','RBM4','RBP4','RELN','RETN','RGS16','SATB1','SCG3','SCG5','SCN1B','SDCCAG8','SEMA3B','SERPINA3','SLPI',
 'SNAP25','SNAP91','SOX9','SPARCL1','SPON1','SPP1','SPRY2','ST6GAL1','STRN','TFPI','TLR8','TMEM100','TMEM158','TMEM176A',
 'TMEM176B','TNFAIP3','TNFRSF1B','TNNT2','TOR1AIP2','TPH1','TRAF1','TRIB1','TRIB2','TSPAN1','TSPAN13','TSPAN7','USH1C',
 'USP12','VWA5A','WDR33','WNT7A','YRDC','ZNF277','ZNF639']

# kras down
kras_dn = ['ABCB11','ABCG4','ACTC1','ADRA2C','AKR1B10','ALOX12B','AMBN','ARHGDIG','ARPP21','ASB7','ATP4A','ATP6V1B1','BARD1',
            'BMPR1B','BRDT','BTG2','C5','CACNA1F','CACNG1','CALCB','CALML5','CAMK1D','CAPN9','CCDC106','CCNA1','CCR8','CD207',
            'CD40LG','CD80','CDH16','CDKAL1','CELSR2','CHRNG','CHST2','CKM','CLDN16','CLDN8','CLPS','CLSTN3','CNTFR','COL2A1',
            'COPZ2','COQ8A','CPA2','CPB1','CPEB3','CYP11B2','CYP39A1','DCC','DLK2','DTNB','EDAR','EDN1','EDN2','EFHD1','EGF',
            'ENTPD7','EPHA5','FGF16','FGF22','FGFR3','FGGY','FSHB','GAMT','GDNF','GP1BA','GP2','GPR19','GPR3','GPRC5C','GRID2',
            'GTF3C5','HNF1A','HSD11B2','HTR1B','HTR1D','IDUA','IFI44L','IFNG','IGFBP2','IL12B','IL5','INSL5','IRS4','ITGB1BP2',
            'ITIH3','KCND1','KCNE2','KCNMB1','KCNN1','KCNQ2','KLHDC8A','KLK7','KLK8','KMT2D','KRT1','KRT13','KRT15','KRT4','KRT5',
            'LFNG','LGALS7','LYPD3','MACROH2A2','MAGIX','MAST3','MEFV','MFSD6','MSH5','MTHFR','MX1','MYH7','MYO15A','MYOT','NGB',
            'NOS1','NPHS1','NPY4R','NR4A2','NR6A1','NRIP2','NTF3','NUDT11','OXT','P2RX6','P2RY4','PAX3','PAX4','PCDHB1','PDCD1',
            'PDE6B','PDK2','PKP1','PLAG1','PNMT','PRKN','PRODH','PROP1','PTGFR','PTPRJ','RGS11','RIBC2','RSAD2','RYR1','RYR2',
            'SCGB1A1','SCN10A','SELENOP','SERPINA10','SERPINB2','SGK1','SHOX2','SIDT1','SKIL','SLC12A3','SLC16A7','SLC25A23',
            'SLC29A3','SLC30A3','SLC38A3','SLC5A5','SLC6A14','SLC6A3','SMPX','SNCB','SNN','SOX10','SPHK2','SPRR3','SPTBN2',
            'SSTR4','STAG3','SYNPO','TAS2R4','TCF7L1','TCL1A','TENM2','TENT5C','TEX15','TFAP2B','TFCP2L1','TFF2','TG','TGFB2',
            'TGM1','THNSL2','THRB','TLX1','TNNI3','TSHB','UGT2B17','UPK3B','VPREB1','VPS50','WNT16','YBX2','YPEL1','ZBTB16',
            'ZC2HC1C','ZNF112']




## box of kras dn
krasdn_ordered = tpm_normed.loc[tpm_normed.index.isin(kras_dn)].apply(lambda x: np.mean(x), axis=1).sort_values().index

kras_plot_dat = tpm_enhancersig[tpm_enhancersig.index.isin(kras_dn)].stack().reset_index()
kras_plot_dat.columns = ['connected_gene', 'mark', 'value']
kras_plot_dat['day'] = kras_plot_dat['mark'].str.split('_').str[0]
kras_plot_dat['rep'] = kras_plot_dat['mark'].str.split('_').str[1]
kras_plot_dat['mark'] = kras_plot_dat['mark'].str.split('_').str[2]

tpmcol  = kras_plot_dat[kras_plot_dat['mark'] == 'log2tpm']
kras_plot_dat = kras_plot_dat[kras_plot_dat['mark'] != 'log2tpm']
kras_plot_dat = pd.merge(kras_plot_dat, tpmcol, on=['connected_gene', 'day', 'rep'], how='left')


kras_plot_dat['connected_gene'] = pd.Categorical(kras_plot_dat['connected_gene'],\
                                                 categories=krasdn_ordered[krasdn_ordered.isin(kras_plot_dat['connected_gene'])],
                                                 ordered=True)
krasdn_20 =  krasdn_ordered[krasdn_ordered.isin(kras_plot_dat['connected_gene'])][:20]


# plot
pplot_dat = kras_plot_dat.loc[kras_plot_dat['connected_gene'].isin(krasdn_20)]
pplot_dat['connected_gene'] = pd.Categorical(pplot_dat['connected_gene'], categories=krasdn_20, ordered=True)

g = sns.catplot(data=pplot_dat, kind="box",palette="pastel",
            x="connected_gene", y="value_x", row="mark_x", col='day', aspect=1.4, 
            showfliers=False, linewidth=1.2)

#https://stackoverflow.com/questions/67309730/how-to-overlay-a-scatterplot-on-top-of-boxplot-with-sns-catplot
g.map_dataframe(sns.stripplot, x="connected_gene", y="value_x", 
                palette=["#404040"],
                alpha=0.6,dodge=True,
                size=3)

g.map_dataframe(sns.stripplot, x="connected_gene", y="value_y", 
                palette=["#CD5555"],
                alpha=0.5,dodge=True,
                s=5, marker="D",
                jitter=False)

for ax in g.axes.flat:
    ax.axhline(0, color='black', linestyle='--', linewidth=0.5)
    
    for label in ax.get_xticklabels():
        label.set_rotation(90)
        label.set_horizontalalignment('center')

g.axes[1, 0].set_ylabel('log2(NSD2i/Vehicle)')

plt.savefig('topdn20_krasdn.pdf', bbox_inches='tight')
plt.close()



### load control gene sets
kras_ctrl_sets = pd.read_csv('../../kras_vsctrl/ctrlsets_from_allnormedtpm.csv')


kras_plot_dat = kras_plot_dat.assign(geneset='krasdn')



ctrl_sets = {}

# Loop to create each DataFrame
for i in range(1, 6):
    key = f'ctrl{i}'
    
    # Use .isin to filter and stack, reset index for the desired operation
    df = tpm_enhancersig.loc[kras_ctrl_sets[kras_ctrl_sets.iloc[:, i-1].isin(tpm_enhancersig.index)].iloc[:, i-1]]\
        .stack().reset_index()
    df.columns = ['connected_gene', 'mark', 'value']
    df['day'] = df['mark'].str.split('_').str[0]
    df['rep'] = df['mark'].str.split('_').str[1]
    df['mark'] = df['mark'].str.split('_').str[2]

    tpmcol = df[df['mark'] == 'log2tpm']
    df = df[df['mark'] != 'log2tpm']
    df = pd.merge(df, tpmcol, on=['connected_gene', 'day', 'rep'], how='left')

    df = df.assign(set=key)

    # Assign the DataFrame to the dictionary
    ctrl_sets[key] = df


ctrl1 = ctrl_sets['ctrl1']
ctrl2 = ctrl_sets['ctrl2']
ctrl3 = ctrl_sets['ctrl3']
ctrl4 = ctrl_sets['ctrl4']
ctrl5 = ctrl_sets['ctrl5']


## concat all
kras_ctrl_ct_tmp = pd.concat([kras_plot_dat, ctrl1, ctrl2, ctrl3, ctrl4, ctrl5], axis=0)


## plot grid distribution



