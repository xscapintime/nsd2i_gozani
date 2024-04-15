# GSEA

rm(list = ls())
#options(scipen = 999)
library(tidyverse)
library(biomaRt)
library(fgsea)
library(stringr)

### gsea
## MSigDB mouse hallmark
args <- commandArgs(TRUE)
gmt <- args[1]
# gmt <- "gmt/c6.all.v2023.2.Hs.symbols.gmt"
pathways <- gmtPathways(gmt)

pathname <- (sub("gmt/","", gmt) %>% str_split(., pattern = ".v", simplify=T))[1]



## DEG table
files <- list.files(path = ".", pattern = "_deg_id.tsv$", recursive = T)

for (i in 1:length(files)) {

    res <- read.csv(files[i], sep = "\t")

    fn <- gsub("_deg_id.tsv", "", files[i])

    # rank plot
    ranks <- deframe(res["stat"])
    names(ranks) <- res$hgnc_symbol

    pdf(paste0(fn, "/", fn, "_stat-rank.pdf"), width = 5, height = 5.5)
    barplot(sort(ranks, decreasing = T))
    dev.off()

    ## gsea
    fgseaRes <- fgsea(pathways = pathways, ranks, minSize=15, maxSize = 500)
    fgseaResTidy <- fgseaRes %>%
        as_tibble() %>%
        arrange(desc(NES))

    for (j in 1:nrow(fgseaResTidy)) {
        fgseaResTidy$leadingEdgechr[j] <- paste(fgseaResTidy$leadingEdge[[j]], collapse=', ')
    }

    
    write.table(fgseaResTidy[,-grep("leadingEdge$",colnames(fgseaResTidy))],
        file = paste0(fn, "/", pathname,"_gsea.tsv"), quote = F, sep = "\t", row.names = F)

    # plot top significant
    pdf(paste0(fn, "/", pathname, "_pathways.pdf"), width = 5, height = 3.5)

    for (i in 1:nrow(fgseaResTidy[fgseaResTidy$pval < 0.005,])) {
        path <- as.character(fgseaResTidy[order(fgseaResTidy$pval),][i,"pathway"])    
        plt  <- plotEnrichment(pathways[[path]], ranks) + labs(subtitle=path)
        print(plt)
    }

    dev.off()

}


## gene list
gene_score <- read.csv("genehancer_score.txt", sep = "\t", header = F)


# convert gene id to symbol
# if SSL error
# https://github.com/grimbough/biomaRt/issues/31#issuecomment-718210434
httr::set_config(httr::config(ssl_verifypeer = FALSE))

## ensembl id and symble
mart <- biomaRt::useEnsembl(biomart = "genes", dataset = "hsapiens_gene_ensembl")




gene_id <- getBM(attributes = c("ensembl_gene_id", "hgnc_symbol"),
                     filters = "ensembl_gene_id",
                     values = res$ensembl_gene_id, mart = mart)




