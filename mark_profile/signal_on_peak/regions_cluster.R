rm(list = ls())
library(tidyverse)
library(ComplexHeatmap)
library(stringr)
library(RColorBrewer)
library(circlize)

## load
mat <- read.csv("k36me2_changes.txt", header = T, row.names = 1, sep = "\t")

mat <- as.matrix(mat)

## density heatmap
pdf("k36me2_density.pdf")
densityHeatmap(mat, title = " ", ylab= "K36me2 log2(NSD2i/DMSO)")
dev.off()


# htcol <- colorRamp2(c(fivenum(mat)[1:4], 0, max(mat)/2, max(mat)), c(rev(brewer.pal(5, "Blues"))[-5],"#FFFFFF", brewer.pal(3, "Reds")[-1]))

## annotation
# peak id
instc_id <- read.csv("../../niceseq_profile/signal_on_peak/intersect_all.bed", header = F, sep = "\t")

# promoter and enhancer
promoters <- read.csv("../../niceseq/peaks/peak_breakdown/cage_intersect_all.bed", sep = "\t", header = F)[,1:3] 
enhancers <- read.csv("../../niceseq/peaks/peak_breakdown/enhancer_intersect_all.bed", sep = "\t", header = F)[,1:3]


promoters <- merge(instc_id, promoters, by.x = c("V1", "V2", "V3"), by.y = c("V1", "V2", "V3")) 
enhancers <- merge(instc_id, enhancers, by.x = c("V1", "V2", "V3"), by.y = c("V1", "V2", "V3"))


## annot heatmap
ann <- data.frame(row.names(mat), 
ifelse(row.names(mat) %in% promoters$V4, "promoter", 
    ifelse(row.names(mat) %in% enhancers$V4, "enhancer", "other")))

colnames(ann) <- c("peakid", "type")
ann$color <- ifelse(ann$type == "promoter", "#c9ddd7", ifelse(ann$type == "enhancer", "#D95F02", "lightgrey"))

annocol <- ann$color
names(annocol) <- ann$type

ha <- rowAnnotation(`Types` = ann$type,
                    border = F,
                    width = unit(.5,"mm"),
                    col = list(`Types` = annocol))


htcol2 <- colorRamp2(c(-4,-1, 0.5), c("#2166AC", "#FFFFFF", "#B2182B"))

ht2 <- Heatmap(mat, name = "log2(NSD2i/DMSO)", row_km = 4, cluster_columns = F, show_row_names=F,
                clustering_distance_rows = "pearson", col = htcol2,
                right_annotation = ha)
pdf("heatmap_v2.pdf", width = 5, height = 6)
draw(ht2)
dev.off()


## enhancer only
enhancermat <- mat[row.names(mat) %in% enhancers$V4,]


ht3 <- Heatmap(mat, name = "log2(NSD2i/DMSO)", row_km = 4, cluster_columns = F, show_row_names=F,
                clustering_distance_rows = "pearson", col = htcol2)
pdf("heatmap_enhancer.pdf", width = 5, height = 6)
draw(ht3)
dev.off()