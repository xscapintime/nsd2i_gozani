# GSEA

rm(list = ls())
#options(scipen = 999)
library(tidyverse)
library(biomaRt)
library(fgsea)
library(stringr)

## gene list enhancer corresponding
genehancer_score <- read.csv("genehancer_score.txt", sep = "\t", header = F)

# rank plot
ranks <- deframe(genehancer_score["V2"])
names(ranks) <- genehancer_score$V1

pdf("stat_rank.pdf", width = 5, height = 5.5)
barplot(sort(ranks, decreasing = T))
dev.off()

# # convert gene id to symbol
# # if SSL error
# # https://github.com/grimbough/biomaRt/issues/31#issuecomment-718210434
# httr::set_config(httr::config(ssl_verifypeer = FALSE))

# ## ensembl id and symble
# humanmart <- biomaRt::useEnsembl(biomart = "genes", dataset = "hsapiens_gene_ensembl")


# gene_id <- getBM(attributes = c("ensembl_gene_id", "hgnc_symbol"),
#                 filters = "ensembl_gene_id",
#                 values = genehancer_score$V1[grep("ENSG", genehancer_score$V1)],
#                 mart = humanmart)
## unconverted ENSG IDs are not in the biomaRt database or TEC




### gsea
## MSigDB
#args <- commandArgs(TRUE)
# gmt <- args[1]
gmtlist <- list.files(path = ".", pattern = ".gmt$", recursive = T)

for (g in 1:length(gmtlist)) {
    gmt <- gmtlist[g]
    pathways <- gmtPathways(gmt)
    pathname <- (sub("gmt/","", gmt) %>% str_split(., pattern = ".v", simplify=T))[1]

    ## gsea
    fgseaRes <- fgsea(pathways = pathways, ranks, minSize=15, maxSize = 500)

    fgseaResTidy <- fgseaRes %>%
        as_tibble() %>%
        arrange(desc(NES))

    for (j in 1:nrow(fgseaResTidy)) {
        fgseaResTidy$leadingEdgechr[j] <- paste(fgseaResTidy$leadingEdge[[j]], collapse=', ')
    }

    # write.table(fgseaResTidy[,-grep("leadingEdge$",colnames(fgseaResTidy))],
    #     file = paste0(pathname,"_gsea.tsv"), quote = F, sep = "\t", row.names = F)


    # # plot top significant
    # pdf(paste0(pathname, "_pathways.pdf"), width = 5, height = 3.5)

    # for (i in 1:nrow(fgseaResTidy[fgseaResTidy$pval < 0.005,])) {
    #     path <- as.character(fgseaResTidy[order(fgseaResTidy$pval),][i,"pathway"])    
    #     plt  <- plotEnrichment(pathways[[path]], ranks) + labs(subtitle=path)
    #     print(plt)
    # }

    # dev.off()


    # plot table
    topPathwaysUp <- fgseaRes[ES > 0][head(order(pval), n=10), pathway]
    topPathwaysDown <- fgseaRes[ES < 0][head(order(pval), n=10), pathway]
    topPathways <- c(topPathwaysUp, rev(topPathwaysDown))
    
    pdf(paste0(pathname, "_topPathways.pdf"), width = 7, height = 10)
    ttb <- plotGseaTable(pathways[topPathways], ranks, fgseaRes,
                  gseaParam=0.5)
    print(ttb)
    dev.off()

    # collapsedPathways <- collapsePathways(fgseaRes[order(pval)][padj < 0.01], 
    #                                       pathways, ranks)
    # mainPathways <- fgseaRes[pathway %in% collapsedPathways$mainPathways][
    #                          order(-NES), pathway]
    
    # pdf(paste0(pathname, "_collapsePathways.pdf"), width = 7, height = 10)
    # ctb <- plotGseaTable(pathways[mainPathways], ranks, fgseaRes, 
    #               gseaParam=0.5)
    # print(ctb)
    # dev.off()


}


