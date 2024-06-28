import os,glob
import pandas as pd
import numpy as np
# import seaborn as sns
# import matplotlib
# from matplotlib import pyplot as plt
# import matplotlib.font_manager
# from statannotations.Annotator import Annotator
# from itertools import product


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



## merged gene set
# gene_uni = set(ct_enhancer_normed.index) & set(tpm_normed.index)
gene_uni = set(tpm_normed.index)


# pathways
KRAS_SIGNALING_UP = ['ABCB1','ACE','ADAM17','ADAM8','ADAMDEC1','ADGRA2','ADGRL4','AKAP12','AKT2','ALDH1A2','ALDH1A3','AMMECR1','ANGPTL4',
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



KRAS_SIGNALING_DN = ['ABCB11','ABCG4','ACTC1','ADRA2C','AKR1B10','ALOX12B','AMBN','ARHGDIG','ARPP21','ASB7','ATP4A','ATP6V1B1','BARD1',
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



EPITHELIAL_MESENCHYMAL_TRANSITION = ['ABI3BP','ACTA2','ADAM12','ANPEP','APLP1','AREG','BASP1','BDNF','BGN','BMP1','CADM1',
                                     'CALD1','CALU','CAP2','CAPG','CCN1','CCN2','CD44','CD59','CDH11','CDH2','CDH6','COL11A1',
                                     'COL12A1','COL16A1','COL1A1','COL1A2','COL3A1','COL4A1','COL4A2','COL5A1','COL5A2','COL5A3',
                                     'COL6A2','COL6A3','COL7A1','COL8A2','COLGALT1','COMP','COPA','CRLF1','CTHRC1','CXCL1','CXCL12',
                                     'CXCL6','CXCL8','DAB2','DCN','DKK1','DPYSL3','DST','ECM1','ECM2','EDIL3','EFEMP2','ELN','EMP3',
                                     'ENO2','FAP','FAS','FBLN1','FBLN2','FBLN5','FBN1','FBN2','FERMT2','FGF2','FLNA','FMOD','FN1','FOXC2',
                                     'FSTL1','FSTL3','FUCA1','FZD8','GADD45A','GADD45B','GAS1','GEM','GJA1','GLIPR1','GPC1','GPX7','GREM1',
                                     'HTRA1','ID2','IGFBP2','IGFBP3','IGFBP4','IL15','IL32','IL6','INHBA','ITGA2','ITGA5','ITGAV','ITGB1',
                                     'ITGB3','ITGB5','JUN','LAMA1','LAMA2','LAMA3','LAMC1','LAMC2','LGALS1','LOX','LOXL1','LOXL2','LRP1',
                                     'LRRC15','LUM','MAGEE1','MATN2','MATN3','MCM7','MEST','MFAP5','MGP','MMP1','MMP14','MMP2','MMP3','MSX1',
                                     'MXRA5','MYL9','MYLK','NID2','NNMT','NOTCH2','NT5E','NTM','OXTR','P3H1','PCOLCE','PCOLCE2','PDGFRB',
                                     'PDLIM4','PFN2','PLAUR','PLOD1','PLOD2','PLOD3','PMEPA1','PMP22','POSTN','PPIB','PRRX1','PRSS2','PTHLH',
                                     'PTX3','PVR','QSOX1','RGS4','RHOB','SAT1','SCG2','SDC1','SDC4','SERPINE1','SERPINE2','SERPINH1',
                                     'SFRP1','SFRP4','SGCB','SGCD','SGCG','SLC6A8','SLIT2','SLIT3','SNAI2','SNTB1','SPARC','SPOCK1',
                                     'SPP1','TAGLN','TFPI2','TGFB1','TGFBI','TGFBR3','TGM2','THBS1','THBS2','THY1','TIMP1','TIMP3','TNC',
                                     'TNFAIP3','TNFRSF11B','TNFRSF12A','TPM1','TPM2','TPM4','VCAM1','VCAN','VEGFA','VEGFC','VIM','WIPF1','WNT5A']



INTERFERON_GAMMA_RESPONSE = ['ADAR','APOL6','ARID5B','ARL4A','AUTS2','B2M','BANK1','BATF2','BPGM','BST2','BTG1','C1R','C1S','CASP1','CASP3',
                             'CASP4','CASP7','CASP8','CCL2','CCL5','CCL7','CD274','CD38','CD40','CD69','CD74','CD86','CDKN1A','CFB','CFH',
                             'CIITA','CMKLR1','CMPK2','CMTR1','CSF2RB','CXCL10','CXCL11','CXCL9','DDX60','DHX58','EIF2AK2','EIF4E3','EPSTI1',
                             'FAS','FCGR1A','FGL2','FPR1','GBP4','GBP6','GCH1','GPR18','GZMA','HELZ2','HERC6','HIF1A','HLA-A','HLA-B','HLA-DMA',
                             'HLA-DQA1','HLA-DRB1','HLA-G','ICAM1','IDO1','IFI27','IFI30','IFI35','IFI44','IFI44L','IFIH1','IFIT1','IFIT2','IFIT3',
                             'IFITM2','IFITM3','IFNAR2','IL10RA','IL15','IL15RA','IL18BP','IL2RB','IL4R','IL6','IL7','IRF1','IRF2','IRF4','IRF5','IRF7',
                             'IRF8','IRF9','ISG15','ISG20','ISOC1','ITGB7','JAK2','KLRK1','LAP3','LATS2','LCP2','LGALS3BP','LY6E','LYSMD2','MARCHF1',
                             'METTL7B','MT2A','MTHFD2','MVP','MX1','MX2','MYD88','NAMPT','NCOA3','NFKB1','NFKBIA','NLRC5','NMI','NOD1','NUP93','OAS2',
                             'OAS3','OASL','OGFR','P2RY14','PARP12','PARP14','PDE4B','PELI1','PFKP','PIM1','PLA2G4A','PLSCR1','PML','PNP','PNPT1','PSMA2',
                             'PSMA3','PSMB10','PSMB2','PSMB8','PSMB9','PSME1','PSME2','PTGS2','PTPN1','PTPN2','PTPN6','RAPGEF6','RBCK1','RIGI','RIPK1',
                             'RIPK2','RNF213','RNF31','RSAD2','RTP4','SAMD9L','SAMHD1','SECTM1','SELP','SERPING1','SLAMF7','SLC25A28','SOCS1','SOCS3','SOD2',
                             'SP110','SPPL2A','SRI','SSPN','ST3GAL5','ST8SIA4','STAT1','STAT2','STAT3','STAT4','TAP1','TAPBP','TDRD7','TNFAIP2','TNFAIP3',
                             'TNFAIP6','TNFSF10','TOR1B','TRAFD1','TRIM14','TRIM21','TRIM25','TRIM26','TXNIP','UBE2L6','UPP1','USP18','VAMP5','VAMP8','VCAM1',
                             'WARS1','XAF1','XCL1','ZBP1','ZNFX1']



E2F_TARGETS = ['AK2','ANP32E','ASF1A','ASF1B','ATAD2','AURKA','AURKB','BARD1','BIRC5','BRCA1','BRCA2','BRMS1L','BUB1B','CBX5','CCNB2','CCNE1','CCP110',
               'CDC20','CDC25A','CDC25B','CDCA3','CDCA8','CDK1','CDK4','CDKN1A','CDKN1B','CDKN2A','CDKN2C','CDKN3','CENPE','CENPM','CHEK1','CHEK2','CIT',
               'CKS1B','CKS2','CNOT9','CSE1L','CTCF','CTPS1','DCK','DCLRE1B','DCTPP1','DDX39A','DEK','DEPDC1','DIAPH3','DLGAP5','DNMT1','DONSON','DSCC1',
               'DUT','E2F8','EED','EIF2S1','ESPL1','EXOSC8','EZH2','GINS1','GINS3','GINS4','GSPT1','H2AX','H2AZ1','HELLS','HMGA1','HMGB2','HMGB3','HMMR',
               'HNRNPD','HUS1','ILF3','ING3','IPO7','JPT1','KIF18B','KIF22','KIF2C','KIF4A','KPNA2','LBR','LIG1','LMNB1','LUC7L3','LYAR','MAD2L1','MCM2',
               'MCM3','MCM4','MCM5','MCM6','MCM7','MELK','MKI67','MLH1','MMS22L','MRE11','MSH2','MTHFD2','MXD3','MYBL2','MYC','NAA38','NAP1L1','NASP','NBN',
               'NCAPD2','NME1','NOLC1','NOP56','NUDT21','NUP107','NUP153','NUP205','ORC2','ORC6','PA2G4','PAICS','PAN2','PCNA','PDS5B','PHF5A','PLK1','PLK4',
               'PMS2','PNN','POLA2','POLD1','POLD2','POLD3','POLE','POLE4','POP7','PPM1D','PPP1R8','PRDX4','PRIM2','PRKDC','PRPS1','PSIP1','PSMC3IP','PTTG1',
               'RACGAP1','RAD1','RAD21','RAD50','RAD51AP1','RAD51C','RAN','RANBP1','RBBP7','RFC1','RFC2','RFC3','RNASEH2A','RPA1','RPA2','RPA3','RRM2','SHMT1',
               'SLBP','SMC1A','SMC3','SMC4','SMC6','SNRPB','SPAG5','SPC24','SPC25','SRSF1','SRSF2','SSRP1','STAG1','STMN1','SUV39H1','SYNCRIP','TACC3','TBRG4',
               'TCF19','TFRC','TIMELESS','TIPIN','TK1','TMPO','TOP2A','TP53','TRA2B','TRIP13','TUBB','TUBG1','UBE2S','UBE2T','UBR7','UNG','USP1','WDR90','WEE1',
               'XPO1','XRCC6','ZW10']


HEDGEHOG = ['CHE','ADGRG1','AMOT','CDK5R1','CDK6','CELSR1','CNTFR','CRMP1','DPYSL2','ETS2','GLI1','HEY1','HEY2','L1CAM','LDB1','MYH9',
            'NF1','NKX6-1','NRCAM','NRP1','NRP2','OPHN1','PLG','PML','PTCH1','RASA1','RTN1','SCG2','SHH','SLIT1','THY1','TLE1','TLE3',
            'UNC5C','VEGFA','VLDLR']


MYC_V1 = ['ABCE1','ACP1','AIMP2','AP3S1','APEX1','BUB3','C1QBP','CAD','CANX','CBX3','CCNA2','CCT2','CCT3','CCT4','CCT5',
          'CCT7','CDC20','CDC45','CDK2','CDK4','CLNS1A','CNBP','COPS5','COX5A','CSTF2','CTPS1','CUL1','CYC1','DDX18',
          'DDX21','DEK','DHX15','DUT','EEF1B2','EIF1AX','EIF2S1','EIF2S2','EIF3B','EIF3D','EIF3J','EIF4A1','EIF4E',
          'EIF4G2','EIF4H','EPRS1','ERH','ETF1','EXOSC7','FAM120A','FBL','G3BP1','GLO1','GNL3','GOT2','GSPT1','H2AZ1',
          'HDAC2','HDDC2','HDGF','HNRNPA1','HNRNPA2B1','HNRNPA3','HNRNPC','HNRNPD','HNRNPR','HNRNPU','HPRT1','HSP90AB1',
          'HSPD1','HSPE1','IARS1','IFRD1','ILF2','IMPDH2','KARS1','KPNA2','KPNB1','LDHA','LSM2','LSM7','MAD2L1','MCM2',
          'MCM4','MCM5','MCM6','MCM7','MRPL23','MRPL9','MRPS18B','MYC','NAP1L1','NCBP1','NCBP2','NDUFAB1','NHP2','NME1',
          'NOLC1','NOP16','NOP56','NPM1','ODC1','ORC2','PA2G4','PABPC1','PABPC4','PCBP1','PCNA','PGK1','PHB1','PHB2',
          'POLD2','POLE3','PPIA','PPM1G','PRDX3','PRDX4','PRPF31','PRPS2','PSMA1','PSMA2','PSMA4','PSMA6','PSMA7',
          'PSMB2','PSMB3','PSMC4','PSMC6','PSMD1','PSMD14','PSMD3','PSMD7','PSMD8','PTGES3','PWP1','RACK1','RAD23B',
          'RAN','RANBP1','RFC4','RNPS1','RPL14','RPL18','RPL22','RPL34','RPL6','RPLP0','RPS10','RPS2','RPS3','RPS5',
          'RPS6','RRM1','RRP9','RSL1D1','RUVBL2','SERBP1','SET','SF3A1','SF3B3','SLC25A3','SMARCC1','SNRPA','SNRPA1',
          'SNRPB2','SNRPD1','SNRPD2','SNRPD3','SNRPG','SRM','SRPK1','SRSF1','SRSF2','SRSF3','SRSF7','SSB','SSBP1',
          'STARD7','SYNCRIP','TARDBP','TCP1','TFDP1','TOMM70','TRA2B','TRIM28','TUFM','TXNL4A','TYMS','U2AF1','UBA2',
          'UBE2E1','UBE2L3','USP1','VBP1','VDAC1','VDAC3','XPO1','XPOT','XRCC6','YWHAE','YWHAQ']

MYC_V2 = ['AIMP2','BYSL','CBX3','CDK4','DCTPP1','DDX18','DUSP2','EXOSC5','FARSA','GNL3','GRWD1','HK2','HSPD1','HSPE1',
          'IMP4','IPO4','LAS1L','MAP3K6','MCM4','MCM5','MPHOSPH10','MRTO4','MYBBP1A','MYC','NDUFAF4','NIP7','NOC4L',
          'NOLC1','NOP16','NOP2','NOP56','NPM1','PA2G4','PES1','PHB1','PLK1','PLK4','PPAN','PPRC1','PRMT3','PUS1',
          'RABEPK','RCL1','RRP12','RRP9','SLC19A1','SLC29A2','SORD','SRM','SUPV3L1','TBRG4','TCOF1','TFB2M','TMEM97',
          'UNG','UTP20','WDR43','WDR74']


PRC2_targ = pd.read_csv('PRC2_targets.txt', header=None, index_col=None, sep='\t')[0].to_list()


## slecting random geens from control d1
def euclidean_distance(row1, row2):
    return np.sqrt(np.sum((row1 - row2) ** 2))

# base expression
all_d1ctrl = tpm.loc[gene_uni,tpm.columns.str.contains('Day1_C')]

for pn in ['KRAS_SIGNALING_UP', 'KRAS_SIGNALING_DN', 'EPITHELIAL_MESENCHYMAL_TRANSITION', 'INTERFERON_GAMMA_RESPONSE', 'E2F_TARGETS',
           'HEDGEHOG', 'MYC_V1', 'MYC_V2', 'PRC2_targ']:
    
    path = eval(pn)

    path = [x for x in path if x in gene_uni]

    path_base_d1 = tpm.loc[path,tpm.columns.str.contains('Day1_C')]

    ctrl_sets = []
    for i in range(path_base_d1.shape[0]):
        distances = all_d1ctrl.apply(lambda row: euclidean_distance(row.values, path_base_d1.iloc[[i]].values), axis=1)
        sorted_indices = np.argsort(distances)

        # keep the 5 rows with the smallest distances
        closest_indices = [idx for idx in sorted_indices if all_d1ctrl.index[idx] in gene_uni][1:6]
        ctrl_sets.append(closest_indices)

    ## convert index to genes
    ctrl_set_genes = [all_d1ctrl.index[idx].to_list() for idx in ctrl_sets]
    ctrl_genes = pd.DataFrame(ctrl_set_genes)

    ctrl_genes.columns = ['ctrl1','ctrl2','ctrl3','ctrl4','ctrl5']
    ctrl_genes[f'{pn}'] = path

    # save to file
    ctrl_genes.to_csv(f'{pn}_ctrl_genes.txt', index=False, sep='\t')

