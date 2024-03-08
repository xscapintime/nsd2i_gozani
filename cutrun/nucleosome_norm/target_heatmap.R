library(tidyverse)
library(ComplexHeatmap)
library(stringr)
library(RColorBrewer)
library(circlize)

## load on-target normaliztion table
ontar_norm <- read.csv('ontar_norm.txt', sep = '\t', header = T, row.names =1 )
ontar_norm <- as.matrix(ontar_norm)

htcol <- colorRamp2(seq(0, 1, length=6), c("#FFC125", "#87CEFF", "#7EC0EE", "#6CA6CD", "#6CA6CD","#4A708B"))


## per mark
marks <- (colnames(ontar_norm) %>% str_split("_", simplify = T))[,4] %>% unique()

for (i in 1:length(marks)) {

    dat <- t(ontar_norm[,grep(marks[i], colnames(ontar_norm))])
    pdf(paste0(marks[i] ,"_ontarget_norm.pdf"), width = 10, height = 4)

    ht <- Heatmap(dat, name = "On-target Specificity",
        col = htcol,
        # left_annotation = ha,
        cluster_rows = F,
        cluster_columns = F,
        width = ncol(dat)*unit(3.7, "mm"),
        height = nrow(dat)*unit(3.7, "mm"))
    draw(ht)
    dev.off()

}
