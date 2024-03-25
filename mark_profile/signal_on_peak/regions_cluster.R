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


htcol <- colorRamp2(c(min(mat), 0, max(mat)),  c("#2E5A87", "white", "#A90C38"))


ht1 <- Heatmap(mat, name = "open chromatin", row_km = 3, cluster_columns = F, show_row_names=F,
clustering_distance_rows = "spearman", col = htcol)
pdf("heatmap.pdf", width = 3, height = 6)
draw(ht1)
dev.off()



pdf("veh_heatmap.pdf", width = 3, height = 6)
ht2 <- Heatmap(veh, name = "open chromatin", row_km = 3, cluster_columns = F)
dev.off()



## extract peak id
row.names(trt) <- peakid
nsd2i_openup_peaks <- row.names(trt[row_order(ht1)[[3]],])

write.table(nsd2i_openup_peaks, "nsd2i_openup_peaks.txt", quote=F, col.names=F, row.names=F)
