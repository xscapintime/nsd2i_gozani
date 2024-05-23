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

## log2tpm
tpm_normed = np.log2(np.divide(*(tpm.iloc[:,np.where(tpm.columns.str.contains('_N'))[0]]+1).\
                    align((tpm.iloc[:,np.where(tpm.columns.str.contains('_C'))[0]]+1), axis=0)))

tpm_normed.columns = tpm_normed.columns.str.replace('N', 'Rep').str.replace('Day', 'D') + '_log2tpm'


# merge
tpm_enhancersig = tpm_normed.merge(ct_enhancer_normed, left_index=True, right_index=True)
# tpm_enhancersig_long = tpm_enhancersig.stack().reset_index()
# tpm_enhancersig_long['group'] = tpm_enhancersig_long['level_1'].str.split('_').str[-1]


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


#### plot #####
## mass spec normed signal

for path in [kras_up,kras_dn]:

    pathname = 'UP' if path == kras_up else 'DN'

    for mark in ['H3K27Ac','H3K36me2', 'H3K27me3']:

        plot_dat = tpm_enhancersig_long[(tpm_enhancersig_long['level_0'].isin(path)) & \
                                    (tpm_enhancersig_long['group'].str.contains(mark))]
    
    
        plot_dat['group'] = pd.Categorical(plot_dat['group'],
                                           categories=np.unique(plot_dat.group),
                                           ordered=True)
    
    
        plt.figure(figsize=(10,8), dpi=300)
        g = sns.boxplot(x='group', y='ctsignal',
                    data=plot_dat,
                    showfliers=False, width=0.7, linewidth=2, saturation=0.8)
        plt.xlabel('')
        plt.ylabel(f'Mass Spec normalized {mark} on enahncers')
        plt.xticks(rotation=330, ha='left', va='top', rotation_mode='anchor')
        plt.title(f'HALLMARK_KRAS_SIGNALING_{pathname}')
    
        plt.tight_layout()
        plt.savefig(f'{mark}_enhancer_{pathname}.pdf')
        plt.close()




plot_dat = tpm_enhancersig.iloc[tpm_enhancersig.index.isin(path),tpm_enhancersig.columns.str.contains('log2tpm|H3K27Ac')]

sns.scatterplot(x=plot_dat['D1_Rep1_log2tpm'], y=plot_dat['D1_Rep1_H3K27Ac'])


## REST
rest_info = tpm_enhancersig.loc['REST'].to_frame()
rest_info['group'] = pd.Categorical(rest_info.index.str.split('_').str[-1])
rest_info['day'] = pd.Categorical(rest_info.index.str.split('_').str[0])


sns.catplot(x='group', y='REST', data=rest_info, hue='group', col='day',
            kind='box')
plt.xticks(rotation=330, ha='left', va='top', rotation_mode='anchor')




mark_df =  tpm_enhancersig_long[(tpm_enhancersig_long['level_0'].isin(path)) & \
                                    (tpm_enhancersig_long['group'].str.contains(mark))]
mark_df['day'] = mark_df['level_1'].str.split('_').str[0] 
mark_df = mark_df.groupby(['day','level_0']).mean().reset_index()

tpm_df = tpm_enhancersig_long[(tpm_enhancersig_long['level_0'].isin(path)) & \
                                    (tpm_enhancersig_long['group'].str.contains('log2tpm'))]
tpm_df['day'] = tpm_df['level_1'].str.split('_').str[0] 

tpm_df = tpm_df.groupby(['day','level_0']).mean().reset_index()






mean_tpm_ensig = tpm_enhancersig.groupby(tpm_enhancersig.columns.str.replace('_Rep[0-9]',''), axis=1).mean()


mean_tpm_ensig.loc[:,mean_tpm_ensig.columns.str.contains(mark)]

mean_tpm_ensig.loc[:,mean_tpm_ensig.columns.str.contains('tpm')]
