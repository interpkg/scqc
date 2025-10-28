
#' scDblFinder for ATAC
#'
#' @param obj seurat object
#' @param fragfile fragments.tsv.gz file
#' @export
#'
Run_scDblFinder_ATAC <- function(obj=NULL, fragfile=NULL, repeats=NULL){
    # amulet
    #repeats =  import('resources/blacklist_repeats_segdups_rmsk_hg38.bed')
    otherChroms <- GRanges(c("chrM","chrX","chrY","MT","M","X","Y"),IRanges(1L,width=10^8)) # check which chromosome notation you are using c("M", "X", "Y", "MT")
    toExclude <- suppressWarnings(c(repeats, otherChroms))
    res <- scDblFinder::amulet(fragfile, regionsToExclude=toExclude)

    # scDblFinder
    # https://github.com/plger/scDblFinder
    # The scDblFinder method can be to single-cell ATACseq (on peak-level counts), however when doing so we recommend using the aggregateFeatures=TRUE parameter (see vignette).
    sce <- as.SingleCellExperiment(obj)
    sce <- scDblFinder(sce, artificialDoublets=1, aggregateFeatures=TRUE, nfeatures=25, processing="normFeatures")

    # Simple p-value combination
    res$scDblFinder.p <- 1-colData(sce)[row.names(res), "scDblFinder.score"]
    res$combined <- apply(res[,c("scDblFinder.p", "p.value")], 1, FUN=function(x){
            x[x<0.001] <- 0.001 # prevent too much skew from very small or 0 p-values
            suppressWarnings(aggregation::fisher(x))
        })

    return(res)
}





