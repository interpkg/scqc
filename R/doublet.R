
#' DoubletFinder
#'
#' @param obj object
#' @param singlet TRUE means remove doublet, and keep only singlet
#' @return object
#' @export
#'
RunDoubletFinder <- function(obj=NULL, singlet=TRUE, outdir='.')
{   
    DefaultAssay(obj) <- 'RNA'
    
    obj <- obj %>%
            NormalizeData() %>%
            FindVariableFeatures(selection.method = "vst", nfeatures = 2000) %>%
            ScaleData() %>%
            RunPCA()

    # Determine percent of variation associated with each PC
    pct <- obj[["pca"]]@stdev / sum(obj[["pca"]]@stdev) * 100
    # Calculate cumulative percents for each PC
    cumu <- cumsum(pct)
    # Determine which PC exhibits cumulative percent greater than 90% and % variation associated with the PC as less than 5
    co1 <- which(cumu > 90 & pct < 5)[1]
    # Determine the difference between variation of PC and subsequent PC
    co2 <- sort(which((pct[1:length(pct) - 1] - pct[2:length(pct)]) > 0.1), decreasing = T)[1] + 1
    # Minimum of the two calculation
    min.pc <- min(co1, co2)

    obj <- obj %>%
            FindNeighbors(reduction="pca", dims = 1:min.pc) %>%
            FindClusters(resolution = 0.1) %>%
            RunUMAP(dims=1:min.pc)
          
    # 10x Single Cell 3' Gene Expression v3.1 assay
    # https://kb.10xgenomics.com/hc/en-us/articles/360001378811-What-is-the-maximum-number-of-cells-that-can-be-profiled-
    n_total_cell = nrow(obj@meta.data)
    doublet_rate = n_total_cell/1000*0.0075
    print(paste0('[INFO] Total cells: ', n_total_cell))
    print(paste0('[INFO] Estimated doublet rate: ', doublet_rate))

    # v2.0.4 'paramSweep' not 'paramSweep_v3'
    sweep.res <- DoubletFinder::paramSweep(obj, PCs = 1:min.pc, sct = FALSE)
    sweep.stats <- DoubletFinder::summarizeSweep(sweep.res, GT = FALSE)
    bcmvn <- DoubletFinder::find.pK(sweep.stats)
    pK <- bcmvn %>% dplyr::filter(BCmetric == max(BCmetric)) %>% dplyr::select(pK)
    pK <- as.numeric(as.character(pK[[1]]))

    # Homotypic Doublet Proportion Estimate 
    homotypic.prop <- DoubletFinder::modelHomotypic(obj@meta.data$seurat_clusters)
    nExp_poi <- round(doublet_rate*nrow(obj@meta.data))
    nExp_poi.adj <- round(nExp_poi*(1-homotypic.prop))

    # Run DoubletFinder with varying classification stringencies 
    # v2.0.4 'doubletFinder' not 'doubletFinder_v3'
    obj <- DoubletFinder::doubletFinder(obj, PCs = 1:min.pc, pN = 0.25, pK = pK, nExp = nExp_poi.adj, reuse.pANN = FALSE, sct = FALSE)
    colnames(obj@meta.data) <- sub("^DF.classifications.*", "doubletfinder", colnames(obj@meta.data))

    # 
    sample_name <- unique(obj$orig.ident)
    pdf(paste0(outdir, '/doubletfinder.', sample_name, '.umap.pdf'), width = 6, height = 5, useDingbats=FALSE)
    p <- DimPlot(obj, group.by = 'doubletfinder')
    print(p)
    dev.off()

    # output for metadata with doublet results
    write.table(obj@meta.data, file = paste0(outdir, '/doubletfinder.', sample_name, '.metadata.xls'), sep = "\t", quote=F, col.names = NA)


    # filter doublet from 
    if (singlet){
        cell.use <- colnames(subset(x = obj, subset = doubletfinder == 'Singlet'))
        obj <- subset(obj, cells=cell.use)
    }
    

    return(obj)
}




#' DoubletFinder 2
#'
#' @param obj object
#' @param singlet TRUE means remove doublet, and keep only singlet
#' @param rate set rate assume total cells=10000
#' @return object
#' @export
#'
RunDoubletFinder2 <- function(
    obj=NULL, 
    singlet=TRUE, 
    rate=0.075, 
    max_dim=30, 
    outdir='.'
){   
    # 2026-01-26
    # Can I perform Cell Hashing in a GEM-X Universal 3’ Gene Expression v4 workflow?
    # https://kb.10xgenomics.com/s/article/33451184917389-Can-I-perform-Cell-Hashing-in-a-GEM-X-Universal-3-Gene-Expression-v4-workflow

    # > SingleCell3_v3
    # Chromium Single Cell 3' Reagent Kits User Guide (v3.1 Chemistry Dual Index)
    # https://cdn.10xgenomics.com/image/upload/v1722285481/support-documents/CG000315_ChromiumNextGEMSingleCell3__GeneExpression_v3.1_DualIndex__RevF.pdf
    # set rate=0.08 is for scRNA & assume total pre-filtered cells are 10000
    # rate <- 0.08
    
    # > SingleCell3_v4
    # Chromium GEM-X Single Cell 3' v4 Gene Expression User Guide
    # https://cdn.10xgenomics.com/image/upload/v1725314293/support-documents/CG000731_ChromiumGEM-X_SingleCell3v4_UserGuide_RevB.pdf
    # rate <- 0.04


    # 2026-01
    # https://github.com/chris-mcginnis-ucsf/DoubletFinder
    DefaultAssay(obj) <- 'RNA'
    obj <- NormalizeData(obj)
    obj <- FindVariableFeatures(obj, selection.method = "vst", nfeatures = 2000)
    obj <- ScaleData(obj)
    obj <- RunPCA(obj)
    obj <- RunUMAP(obj, dims = 1:max_dim)


    ## pK Identification (no ground-truth) 
    sweep.res <- DoubletFinder::paramSweep(obj, PCs = 1:max_dim, sct = FALSE)
    sweep.stats <- DoubletFinder::summarizeSweep(sweep.res, GT = FALSE)
    bcmvn <- DoubletFinder::find.pK(sweep.stats)
    optimal.pk <- as.numeric(as.character(bcmvn$pK[which.max(bcmvn$BCmetric)]))

    # 10x Single Cell 3' Gene Expression v3.1 assay
    # https://kb.10xgenomics.com/s/article/33451184917389-Can-I-perform-Cell-Hashing-in-a-GEM-X-Universal-3-Gene-Expression-v4-workflow
    # 
    n_total_cell = nrow(obj@meta.data)
    doublet_rate = n_total_cell/10000*rate
    print(paste0('[INFO] Total cells: ', n_total_cell))
    print(paste0('[INFO] Estimated doublet rate: ', doublet_rate))
    print(paste0('[INFO] pK: ', optimal.pk))

    ## Homotypic Doublet Proportion Estimate 
    homotypic.prop <- DoubletFinder::modelHomotypic(obj@meta.data$seurat_clusters)
    nExp_poi <- round(doublet_rate*nrow(obj@meta.data))
    nExp_poi.adj <- round(nExp_poi*(1-homotypic.prop))

    ## Run DoubletFinder with varying classification stringencies 
    obj <- DoubletFinder::doubletFinder(obj, PCs = 1:max_dim, pN = 0.25, pK = optimal.pk, nExp = nExp_poi.adj, reuse.pANN = FALSE, sct = FALSE)
    colnames(obj@meta.data) <- sub("^DF.classifications.*", "doubletfinder", colnames(obj@meta.data))


    # plot
    sample_name <- unique(obj$orig.ident)
    p <- Seurat::DimPlot(obj, group.by = 'doubletfinder')
    ggplot2::ggsave(paste0(outdir, '/doubletfinder.', sample_name, '.umap.pdf'), width = 6, height = 5)

    # output for metadata with doublet results
    write.table(obj@meta.data, file = paste0(outdir, '/doubletfinder.', sample_name, '.metadata.xls'), sep = "\t", quote=F, col.names = NA)

    # filter doublet from 
    if (singlet){
        cell.use <- colnames(subset(x = obj, subset = doubletfinder == 'Singlet'))
        obj <- subset(obj, cells=cell.use)
    }
    

    return(obj)
}





