
#' DropletQC
#'
#' @param dir cellranger out dir
#' @export
#'
RunDropletQC <- function(dir=NULL, form='cellranger'){
	if (form == 'cellranger'){
		dnf <- DropletQC::nuclear_fraction_tags(outs = dir, tiles = 1, cores = 1, verbose = FALSE)
	} 

	if (form == 'cellranger-arc'){
		bam <- paste0(dir, '/gex_possorted_bam.bam')
		barcodes <- paste0(dir, '/filtered_feature_bc_matrix/barcodes.tsv.gz')
		dnf <- DropletQC::nuclear_fraction_tags(bam = bam, barcodes=barcodes, tiles = 1, cores = 1, verbose = FALSE)
	}
    
	write.table(dnf, file=paste0(dir, '/nuclear_fraction.tsv'), sep='\t', quote=FALSE, col.names=FALSE)
}




#' DropletQC empty drops
#'
#' @param fnf dropletqc data
#' @param fmeta metadata
#' @param outdir out dir
#' @export
#'
RunEmptyDrop <- function(
	nf=NULL, meta=NULL, 
	nf_rescue=0.05, 
	umi_rescue=1000, 
	include_plot=FALSE,
	col_ct='', 
	outdir='.'
){
	droplet <- read.table(nf, sep='\t', header=F)
	metadata <- read.table(meta, sep='\t', header=T, row.names=1)
	dir.create(outdir)

	metadata$nf <- droplet$V2[match(rownames(metadata), droplet$V1)]
	nf_umi <- metadata[, c('nf', 'nCount_RNA')]
	colnames(nf_umi) <- c('nf', 'umi_count')
	write.table(nf_umi, file=paste0(outdir,'/nf_umi.tsv'), sep='\t', quote=FALSE, col.names=NA)


	# Identify drop cells
	# default: nf_rescue = 0.05,umi_rescue = 1000
	nfc.ed <- DropletQC::identify_empty_drops(nf_umi=nf_umi,
	                           		nf_rescue = nf_rescue, umi_rescue = umi_rescue,
	                                include_plot = include_plot)
	write.table(nfc.ed, file=paste0(outdir,'/nf_umi.empty_drops.tsv'), sep='\t', quote=FALSE, col.names=NA)


	# Identify damaged cells
	if (nchar(col_ct) > 0){
		nfc.ed$cell_type <- metadata[[col_ct]]

		nfc.ed_dc <- DropletQC::identify_damaged_cells(nfc.ed, verbose = FALSE)
		write.table(nfc.ed_dc, file=paste0(outdir,'/nf_umi.damaged_cells.tsv'), sep='\t', quote=FALSE, col.names=NA)
	}
}




#' DropletQC - define and remove empty-drop
#'
#' @param dir cellranger out dir
#' @param obj seurat object
#' @export
#'
RunDelEmptyDrop <- function(dir=NULL, obj=NULL, outdir='.'){
    dnf <- DropletQC::nuclear_fraction_tags(outs = dir, tiles = 1, cores = 1, verbose = FALSE)
    #   nuclear_fraction
	write.table(dnf, file=paste0(dir, '/nuclear_fraction.tsv'), sep='\t', quote=FALSE, col.names=FALSE)

	meta <- obj@meta.data

	meta$nf <- dnf$nuclear_fraction[match(rownames(meta), rownames(dnf))]
	nf_umi <- meta[, c('nf', 'nCount_RNA')]
	colnames(nf_umi) <- c('nf', 'umi_count')

	nf_ed <- DropletQC::identify_empty_drops(nf_umi=nf_umi)
	write.table(nf_ed, file=paste0(outdir,'/nf_umi.empty_drops.tsv'), sep='\t', quote=FALSE, col.names=NA)

	# real cell (filtered emptydrop)
	real_cell <- rownames(nf_ed)[nf_ed$cell_status == 'cell']
	obj <- subset(obj, cells=real_cell)

	return(obj)
}



