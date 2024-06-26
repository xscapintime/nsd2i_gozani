require(Organism.dplyr)
require(GenomicRanges)
require(BSgenome.Hsapiens.UCSC.hg38)

rush_promo <- function(db, upstream, downstream, genes, sequence) {
    prom <- promoters(db, upstream = upstream, downstream = downstream, filter = SymbolFilter(genes), use.names=TRUE)
    prom <- trim(prom)
    prom <- split(prom, prom$symbol)

    if (sequence == T) {
        prom_seq <- getSeq(Hsapiens, prom)
        #names(prom_seq) <- prom$symbol
        return(prom_seq)
    }
    else {
        return(prom)
    }
}
