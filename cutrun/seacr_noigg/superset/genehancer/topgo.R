rm(list = ls())
library(topGO)
library(org.Mm.eg.db)
library(stringr)
library(ggplot2)
library(ggthemes)
library(tidyverse)


## gene list enhancer corresponding
genehancer_score <- read.csv("genehancer_score.txt", sep = "\t", header = F)


## read and GO BP enrichment
# seperate up and down
detype <- "Up"

upgores_list <- list()
for (i in 1:length(detbs)) {
    res <- read.csv(detbs[i], sep = "\t", header=T)


    gene_universe <- as.numeric(ifelse(res$de == detype,TRUE,FALSE)) %>%
        factor() %>% setNames(res$geneIDs)

    GOdata <- new("topGOdata",
                    ontology = "BP",
                    allGenes = gene_universe,
                    annot = annFUN.org,
                    mapping = "org.Mm.eg",
                    ID = "ensembl")

    res_fisher <- runTest(GOdata, algorithm = "classic", statistic = "fisher")
    allres_bp <- GenTable(GOdata, classicFisher = res_fisher,
                            orderBy = "classicFisher",
                            ranksOf = "classicFisher",
                            topNodes = length(usedGO(GOdata)),numChar=1000)
    upgores_list[[i]] <- allres_bp
    
    print(paste0("Exporting...", detype, "...", gsub("_transcript_de.txt", "", detbs[i])))
    write.table(allres_bp, file = gsub("_transcript_de.txt",paste0("_gobp_",detype,".txt"), detbs[i]), quote = F, sep = "\t")

}

names(upgores_list) <- gsub("_transcript_de.txt",paste0("_gobp_",detype), detbs) %>% gsub("/","_",.)


## down DEGs
detype <- "Down"

downgores_list <- list()
for (i in 1:length(detbs)) {
    res <- read.csv(detbs[i], sep = "\t", header=T)


    gene_universe <- as.numeric(ifelse(res$de == detype,TRUE,FALSE)) %>%
        factor() %>% setNames(res$geneIDs)

    GOdata <- new("topGOdata",
                    ontology = "BP",
                    allGenes = gene_universe,
                    annot = annFUN.org,
                    mapping = "org.Mm.eg",
                    ID = "ensembl")

    res_fisher <- runTest(GOdata, algorithm = "classic", statistic = "fisher")
    allres_bp <- GenTable(GOdata, classicFisher = res_fisher,
                            orderBy = "classicFisher",
                            ranksOf = "classicFisher",
                            topNodes = length(usedGO(GOdata)),numChar=1000)
    downgores_list[[i]] <- allres_bp
    print(paste0("Exporting...", detype, "...", gsub("_transcript_de.txt", "", detbs[i])))
    write.table(allres_bp, file = gsub("_transcript_de.txt",paste0("_gobp_",detype,".txt"), detbs[i]), quote = F, sep = "\t")

}

names(downgores_list) <- gsub("_transcript_de.txt",paste0("_gobp_",detype), detbs) %>% gsub("/","_",.)

 


## dotplot 
# https://www.biostars.org/p/471549/
# set theme
theme_set(theme_bw(base_size = 24) +
  theme(
    plot.title = element_text(angle = 0, size = 16, vjust = 1),
    plot.subtitle = element_text(angle = 0, size = 14, vjust = 1),
    plot.caption = element_text(angle = 0, size = 12, vjust = 1),

    axis.text.x = element_text(angle = 0, size = 14, hjust = 1.10, colour = "black"),
    axis.text.y = element_text(angle = 0, size = 14, vjust = 0.5, colour = "black"),
    axis.title = element_text(size = 14),
    axis.title.x = element_text(size = 14, colour = "black"),
    axis.title.y = element_text(size = 14, colour = "black"),
    # axis.line = element_line(colour = 'black'),

    #Legend
    legend.position = 'right',
    legend.background = element_rect(),
    legend.key = element_blank(), # removes the border
    legend.key.size = unit(1, "cm"), # Sets overall area/size of the legend
    legend.text = element_text(size = 14), # Text size
    title = element_text(size = 14))
)



### for ups
for (j in 1:length(upgores_list)) {

ntop <- 20
ggdata <- (upgores_list[[j]] %>% filter(Annotated != 1 & Significant !=1))[1:ntop,]
ggdata$Term <- factor(ggdata$Term, levels = unique(rev(ggdata$Term))) # fixes order
ggdata$classicFisher <- as.numeric(ggdata$classicFisher)


p <- ggplot(ggdata,
  aes(y = Term, x = -log10(classicFisher), size = Significant/Annotated, fill = -log10(classicFisher))) +

  expand_limits(x = c(1,3.5)) +
  geom_point(shape = 21) +
  scale_size(range = c(2.5,12.5)) +
  scale_fill_continuous(low = 'royalblue', high = 'red4') +

  ylab('') + xlab('Enrichment score') +
  labs(
    title = 'GO biological process',
    #subtitle = 'Top 50 terms ordered by Kolmogorov-Smirnov p-value',
    subtitle = 'Top 20 siginificant terms',
    caption = 'Cut-off lines drawn at equivalents of p=0.05, p=0.01, p=0.001',
    fill = "-log10(p-value)") +

  guides(fill = guide_colourbar(order = 1),
         size = guide_legend(order = 2)) +

  geom_vline(xintercept = c(-log10(0.05), -log10(0.01), -log10(0.001)),
    linetype = c("dashed", "longdash", "solid"),
    colour = c("black", "black", "black"),
    size = c(0.5, 1, 2)) +
    coord_fixed()

p_fixed <- egg::set_panel_size(p, width  = unit(4, "in"), height = unit(7, "in"))


ggsave(paste0(str_split(names(upgores_list[j]),"_",simplify = T)[1], "/", 
        names(upgores_list[j]), ".pdf"), p_fixed, 
                device = NULL,
                height = 10,
                width = 16)

}



### for downs
for (j in 1:length(downgores_list)) {

ntop <- 20
ggdata <- (downgores_list[[j]] %>% filter(Annotated != 1 & Significant !=1))[1:ntop,]
ggdata$Term <- factor(ggdata$Term, levels = unique(rev(ggdata$Term))) # fixes order
ggdata$classicFisher <- as.numeric(ggdata$classicFisher)


p <- ggplot(ggdata,
  aes(y = Term, x = -log10(classicFisher), size = Significant/Annotated, fill = -log10(classicFisher))) +

  expand_limits(x = c(1,3.5)) +
  geom_point(shape = 21) +
  scale_size(range = c(2.5,12.5)) +
  scale_fill_continuous(low = 'royalblue', high = 'red4') +

  ylab('') + xlab('Enrichment score') +
  labs(
    title = 'GO biological process',
    #subtitle = 'Top 50 terms ordered by Kolmogorov-Smirnov p-value',
    subtitle = 'Top 20 siginificant terms',
    caption = 'Cut-off lines drawn at equivalents of p=0.05, p=0.01, p=0.001',
    fill = "-log10(p-value)") +

  guides(fill = guide_colourbar(order = 1),
         size = guide_legend(order = 2)) +

  geom_vline(xintercept = c(-log10(0.05), -log10(0.01), -log10(0.001)),
    linetype = c("dashed", "longdash", "solid"),
    colour = c("black", "black", "black"),
    size = c(0.5, 1, 2)) +
    coord_fixed()

p_fixed <- egg::set_panel_size(p, width  = unit(4, "in"), height = unit(7, "in"))


ggsave(paste0(str_split(names(downgores_list[j]),"_",simplify = T)[1], "/", 
        names(downgores_list[j]), ".pdf"), p_fixed, 
                device = NULL,
                height = 10,
                width = 16)

}


