
#' scDblFinder for ATAC
#'
#' @param obj seurat object
#' @param fragfile fragments.tsv.gz file
#' @param repeats blacklist GRanges form
#' @export
#'
Run_scDblFinder_ATAC <- function(obj=NULL, fragfile=NULL, repeats=NULL){
    # amulet
    #repeats =  import('resources/blacklist_repeats_segdups_rmsk_hg38.bed')
    otherChroms <- GRanges(c("chrM","chrX","chrY","MT","M","X","Y"),IRanges(1L,width=10^8))
    toExclude <- suppressWarnings(c(repeats, otherChroms))
    df <- scDblFinder::amulet(fragfile)

    # scDblFinder
    # https://github.com/plger/scDblFinder
    # The scDblFinder method can be to single-cell ATACseq (on peak-level counts), however when doing so we recommend using the aggregateFeatures=TRUE parameter (see vignette).
    sce <- as.SingleCellExperiment(obj)
    sce <- scDblFinder::scDblFinder(sce, artificialDoublets=1, aggregateFeatures=TRUE, nfeatures=25, processing="normFeatures")

    # Simple p-value combination
    df$scDblFinder.p <- 1 - SingleCellExperiment::colData(sce)[row.names(df), "scDblFinder.score"]
    df$combined <- apply(df[,c("scDblFinder.p", "p.value")], 1, FUN=function(x){
            x[x<0.001] <- 0.001
            suppressWarnings(aggregation::fisher(x))
        })

    # Conservative (fewer false positives)
    df$is_doublet <- df$combined < 0.01
    # data.frame
    # <cell>  nFrags  uniqFrags  nAbove2  total.nAbove2  p.value  q.value  scDblFinder.p  combinded  is_doublet

    # single cells (filtered doublet cells)
    #cells_singl <- rownames(df[!df$is_doublet,])

    return(df)
}









