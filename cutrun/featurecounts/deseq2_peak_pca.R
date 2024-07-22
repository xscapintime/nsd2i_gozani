rm(list = ls())

library(tidyverse)
library(DESeq2)
library(stringr)
library(ggplot2)

# load table
dat <- read.table("H3K27Ac.counts.txt", header = T, row.names = 1, skip = 1, sep = "\t")
dat <- dat[,6:17]

colnames(dat) <- str_split(colnames(dat), "\\.", simplify = T)[,13]
colnames(dat) <- str_replace(colnames(dat), "H3K27Ac", "H3K27ac")

# meta
meta <- as.data.frame(colnames(dat))
meta$timepoint <- str_split(colnames(dat),"_", simplify = T)[,2]
meta$treatment <- str_split(colnames(dat),"_", simplify = T)[,1]
meta$treatment <- str_replace(meta$treatment, "Vehicle", "DMSO")
# meta <- meta[,-1]

meta$treatment <- factor(meta$treatment, levels = c("NSD2i", "DMSO"))
meta$timepoint <- factor(meta$timepoint, levels = c("D1", "D5", "D9"))


## deseq2
dds <- DESeqDataSetFromMatrix(countData=dat,colData=meta, design= ~ treatment + timepoint)

keep <- rowSums(counts(dds)) >= 10

dds.filt <- dds[keep,]

vsd <- vst(dds.filt)

pca_dt <- DESeq2::plotPCA(vsd, intgroup = "treatment", ntop = 500,returnData = T)
pca_dt$timepoint <- meta$timepoint

percentVar <- round(100 * attr(pca_dt, "percentVar"))


# Plot
pca_dt %>%
  ggplot(aes(x = pca_dt[,1], y = pca_dt[,2], color = timepoint,shape=treatment)) +
  geom_point(size = 3, show.legend = T,position=position_jitter(width=10)) +
  geom_hline(yintercept = 0,alpha=0.1) +
  geom_vline(xintercept = 0,alpha=0.1) +
  xlab(label = paste(colnames(pca_dt[2]))) +
  ylab(label = paste(colnames(pca_dt[3]))) +
  labs(title = "Top 500 most variable peaks",colour="") +
  # scale_shape_manual(values=c(16,18,17),name="Day",labels=c("Day 1","Day 5","Day 9")) +
  scale_shape_manual(values=c(16,17),name="Treatment") +
  # scale_color_manual(values = c("blue","navyblue","purple","orange","red","maroon"),name="Condition") +
  scale_color_manual(values = c("blue","forestgreen","red"),name="Timepoint",labels=c("Day 1","Day 5","Day 9")) +
  theme(panel.border = element_rect(colour = "black", fill=NA, linewidth=1),
        panel.background = element_blank(),
        panel.grid.major = element_blank(),
        panel.grid.minor = element_blank(),
        axis.title = element_text(size=12,family="Helvetica",colour = "black"),
        axis.text.y=element_text(size=12,family="Helvetica",colour = "black"),
        axis.text.x=element_text(size=12,family="Helvetica",colour = "black"),
        axis.line = element_line(linewidth = 0.5, linetype = "solid",colour = "black"),
        plot.title = element_text(hjust = 0.5,color = "black",size=12,family="Helvetica"),
        strip.background =element_rect(fill="white"),
        strip.text = element_text(
        size = 12, color = "black",family = "Helvetica"),
        legend.text=element_text(size=12,family = "Helvetica",color = "black"),
        legend.title=element_text(size=12,family="Helvetica",color="black"),
        legend.background = element_rect(fill="white"),
        legend.key=element_rect(fill="white", colour = NA),
        legend.key.size = unit(0.1, 'cm'),
        plot.margin = unit(c(1.2,1.2,1.2,1.2), "mm"),
        panel.spacing = unit(0.1,'cm'),
        panel.spacing.y = unit(0.1,'cm'),
        panel.spacing.x = unit(0.1,'cm')) +
    xlab(paste0("PC1: ",percentVar[1],"% variance")) +
    ylab(paste0("PC2: ",percentVar[2],"% variance"))
# legend.position="right",
#     legend.justification="right",
#     legend.box.spacing = unit(-0.001, "cm"))
#
ggsave(filename = "top500variablepeak.H3K27ac.PCA.pdf",device = "pdf",units = "cm",width = 15,height=12,dpi = 600,bg="white")
    