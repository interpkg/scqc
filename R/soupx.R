
#' SoupX
#'
#' @param dir cellranger out dir
#' @export
#'
RunSoupX <- function(dir=NULL, adj=FALSE){
    sc <- SoupX::load10X(dir)
    pdf(paste0(dir, '/soupx.raw.pdf'), w=6, h=6)
    sc <- SoupX::autoEstCont(sc)
    dev.off()

    # adjust counts
    if (adj){
        adj.matrix <- SoupX::adjustCounts(sc, roundToInt = TRUE)
        DropletUtils:::write10xCounts(paste0(dir,'/soupx_adjustCounts.h5'), adj.matrix, type="HDF5", overwrite=TRUE)
    }
}




