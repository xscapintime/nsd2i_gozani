rm(list = ls())
library(tidyverse)
library(ComplexHeatmap)


## load
files <- list.files(path="data", "*.ave.intsc.bed")

sig_data <- lapply(files, function(file) {
               read.table(paste0("data/", file), sep = "\t", header = F)}[,5])

spnames <- sub(".ave.intsc.bed", "", files) %>% sub("MiaPaCa2.", "", .)

m_mat <- Reduce(cbind, sig_data)
colnames(m_mat) <- spnames

## load rownames
peakid <- read.csv(paste0("data/", files[1]), sep = "\t", header = F)[,1]


## split
trt <- m_mat[,1:5]
veh <- m_mat[,6:11]

## density heatmap
pdf("nsd2i_density.pdf")
densityHeatmap(trt, ylim= c(0,2.5), title = " ", ylab= "NSD2i")
dev.off()

pdf("veh_density.pdf")
densityHeatmap(veh, ylim= c(0,2.5), title = " ", ylab= "Vehicle")
dev.off()


ht1 <- Heatmap(trt, name = "open chromatin", row_km = 3, cluster_columns = F)
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
