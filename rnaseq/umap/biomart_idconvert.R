library(biomaRt)
library(tidyverse)

options(digits = 22)

# if SSL error
# https://github.com/grimbough/biomaRt/issues/31#issuecomment-718210434
httr::set_config(httr::config(ssl_verifypeer = FALSE))

## ensembl id and symble
mart <- biomaRt::useEnsembl(biomart = "genes", dataset = "hsapiens_gene_ensembl")


## convert
human_id <- read.csv("../tpm/tximport-tpm.csv", header = T, stringsAsFactors = F)[,1]
human_symbol_df <- getBM(attributes = c("ensembl_gene_id", "hgnc_symbol"),
                      filters = "ensembl_gene_id",
                      values = human_id, mart = mart)
human_symbol_df$hgnc_symbol <- ifelse(human_symbol_df$hgnc_symbol == "", NA, human_symbol_df$hgnc_symbol)


unq_df  <- human_symbol_df %>% na.omit() %>% distinct()

write.table(unq_df, file = "human_ensembl_syb.tsv", quote = F, sep = "\t", row.names = F)
