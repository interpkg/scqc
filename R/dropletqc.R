
#' DropletQC
#'
#' @param dir cellranger out dir
#' @export
#'
RunDropletQC <- function(dir=NULL, form='cellranger', outdir=''){
	dnf <- NULL

	if (form == 'cellranger'){
		dnf <- DropletQC::nuclear_fraction_tags(outs = dir, tiles = 1, cores = 1, verbose = FALSE)
	} 

	if (form == 'cellranger-arc'){
		bam <- paste0(dir, '/gex_possorted_bam.bam')
		barcodes <- paste0(dir, '/filtered_feature_bc_matrix/barcodes.tsv.gz')
		dnf <- DropletQC::nuclear_fraction_tags(bam = bam, barcodes=barcodes, tiles = 1, cores = 1, verbose = FALSE)
	}
    

	if (nchar(outdir) == 0){
		outdir <- dir
	} else {
		dir.create(outdir)
	}
	
	write.table(dnf, file=paste0(outdir, '/nuclear_fraction.tsv'), sep='\t', quote=FALSE, col.names=FALSE)
	
}




#' DropletQC empty drops
#'
#' @param fnf dropletqc data file
#' @param meta metadata data-frame
#' @param outdir out dir
#' @export
#'
EstimateEmptyDrops <- function(
	fnf=NULL, meta=NULL, 
	nf_rescue=0.05, 
	umi_rescue=1000, 
	include_plot=FALSE,
	col_ct='', 
	outdir='.'
){
	dir.create(outdir)

	dnf <- read.table(fnf, sep='\t', header=F)
	#meta <- read.table(fmeta, sep='\t', header=T, row.names=1)

	colnames(dnf) <- c('cell', 'nf')

	meta$nf <- dnf$nf[match(rownames(meta), dnf$cell)]
	nf_umi <- meta[, c('nf', 'nCount_RNA')]
	colnames(nf_umi) <- c('nf', 'umi_count')
	write.table(nf_umi, file=paste0(outdir,'/nf_umi.tsv'), sep='\t', quote=FALSE, col.names=NA)


	# Identify drop cells
	# default: nf_rescue = 0.05,umi_rescue = 1000
	nfc.ed <- DropletQC::identify_empty_drops(nf_umi=nf_umi,
	                           		nf_rescue = nf_rescue, umi_rescue = umi_rescue,
	                                include_plot = include_plot)
	write.table(nfc.ed, file=paste0(outdir,'/nf_umi.empty_drops.tsv'), sep='\t', quote=FALSE, col.names=NA)

	d_log <- as.data.frame(table(nfc.ed$cell_status))
	write.table(d_log, file=paste0(outdir,'/nf_umi.empty_drops.log'), sep='\t', quote=FALSE, row.names=F, col.names=F)


	# Identify damaged cells
	if (nchar(col_ct) > 0){
		nfc.ed$cell_type <- meta[[col_ct]]

		nfc.ed_dc <- DropletQC::identify_damaged_cells(nfc.ed, verbose = FALSE)
		write.table(nfc.ed_dc, file=paste0(outdir,'/nf_umi.damaged_cells.tsv'), sep='\t', quote=FALSE, col.names=NA)

		d_log <- as.data.frame(table(nfc.ed_dc$cell_status))
		write.table(d_log, file=paste0(outdir,'/nf_umi.damaged_cells.log'), sep='\t', quote=FALSE, row.names=F, col.names=F)
	}
}




#' Filter empty droplet using DropletQC 
#'
#' @param dir cellranger out dir
#' @param obj seurat object
#' @export
#'
RunFilterEmptyDrops <- function(dir=NULL, obj=NULL, form='cellranger', outdir='.'){
	meta <- obj@meta.data

	RunDropletQC(dir=dir, form=form, outdir=outdir)
	EstimateEmptyDrops(fnf=paste0(outdir, '/nuclear_fraction.tsv'), meta=meta, outdir=outdir)

	empty_drop <- read.table(paste0(outdir, '/nf_umi.empty_drops.tsv'), sep='\t', header=T, row.names=1)

	# real cell (filtered emptydrop)
	real_cell <- rownames(empty_drop)[empty_drop$cell_status == 'cell']
	obj <- subset(obj, cells=real_cell)

	return(obj)
}



