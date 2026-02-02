
#' SoupX
#'
#' @param dir cellranger out dir
#' @export
#'
RunSoupX <- function(dir=NULL, adj=FALSE){
    out_dir <- paste0(dir, '/soupx')
    dir.create(out_dir)

    # read
    sc <- SoupX::load10X(dir)

    # out
    pdf(paste0(out_dir, '/soupx.raw.pdf'), w=6, h=6)
    sc <- SoupX::autoEstCont(sc)
    dev.off()

    # soupx- Estimated global rho
    write(unique(sc$metaData$rho), file=paste0(out_dir, '/soupx.raw.global_rho.out'))

    # adjust counts
    if (adj){
        adj.matrix <- SoupX::adjustCounts(sc, roundToInt = TRUE)
        DropletUtils:::write10xCounts(paste0(out_dir,'/soupx_adjustCounts.h5'), adj.matrix, type="HDF5", overwrite=TRUE)
    }
}




