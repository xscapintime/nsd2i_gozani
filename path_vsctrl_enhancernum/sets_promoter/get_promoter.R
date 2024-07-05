rm(list = ls())
library(tidyverse)

library(TxDb.Hsapiens.UCSC.hg38.knownGene)
library(Organism.dplyr)
library(BSgenome.Hsapiens.UCSC.hg38)
library(org.Hs.eg.db)
library(stringr)

source("rush_promo.R")


src <- src_ucsc("Homo sapiens")

## list gene sets and ctrls
setfiles <- list.files("../ctrlsets/", pattern = "*_ctrl_genes.txt")

sets <- lapply(setfiles, function(file) read.table(paste0("../ctrlsets/", file), sep = "\t",  header = T))
names(sets) <- str_split(setfiles, "_ct", simplify = T)[,1]



## get promoter bed
for (i in 1:length(sets)) {
    setname = names(sets[i])
    
    genelist = sets[[i]][setname][,1]

    prom <- rush_promo(db = src, upstream = 1000, downstream = 0, genes = genelist, sequence = F)
    rtracklayer::export.bed(unlist(prom), paste0("prom_bed/", setname, "_prom_1k.bed"))

    for (j in 1:(ncol(sets[[i]])-1)) {
        ctrlset <- sets[[i]][,j]
    
        prom <- rush_promo(db = src, upstream = 1000, downstream = 0, genes = ctrlset, sequence = F)
        rtracklayer::export.bed(unlist(prom), paste0("prom_bed/", setname, "_", colnames(sets[[i]])[j], "_prom_1k.bed"))

    }
    

}
