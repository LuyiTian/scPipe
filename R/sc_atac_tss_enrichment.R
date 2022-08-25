
# bamfile <- "/Users/voogd.o/Documents/pythonTesting/scPipeATAC/demux_completematch_small_chr21_R1_aligned_tagged_sorted.bam"
# outPath <- "~/Documents/pythonTesting/splitedCustom"
# # what file does this need to be ?
# gff <- "/Users/voogd.o/Documents/pythonTesting/scPipeATAC/ens_tiny_anno.gff3"
#' @title SC ATAC Plot TSS Enrichment
#' @description Generate a TSS enrichment plot displaying TSS scores aggregated over reads centered on the TSSs.
#' 
#' @param bam string path to input bam file
#' @param outdir string path to output directory
#' @param txs GRanges transcript object. Either one of `txs` and `seqinfo`, `gff` or `TxDb` file must be non null
#' @param seqinfo SeqInfo object containing sequence lengths. Either one of `txs` and `seqinfo`, `gff` or `TxDb` file must be non null
#' @param gff optional string to gff file path. Either one of `txs` and `seqinfo`, `gff` or `TxDb` file must be non null
#' @param TxDb TxDb type object containing sequence length data. Either one of `txs` and `seqinfo`, `gff` or `TxDb` file must be non null
#' @param seqnames optional character vector of sequence names to subset the TxDb or GFF file
#' @details The TSS enrichment calculation is a signal to noise calculation. 
#' The reads around a reference set of TSSs are collected to form an aggregate distribution of reads centered on the 
#' TSSs and extending to 2000 bp in either direction (for a total of 4000bp). This distribution is then normalized by
#' taking the average read depth in the 100 bps at each of the end flanks of the distribution (for a total of 200bp of 
#' averaged data) and calculating a fold change at each position over that average read depth. This means that the flanks 
#' should start at 1, and if there is high read signal at transcription start sites (highly open regions of the genome) 
#' there should be an increase in signal up to a peak in the middle. We take the signal value at the center of the 
#' distribution after this normalization as our TSS enrichment metric.
#' 
#' @importFrom ATACseqQC readBamFile shiftGAlignmentsList TSSEscore
#' @importFrom GenomicFeatures makeTxDbFromGFF transcripts seqinfo
#' @importFrom Rsamtools scanBam BamFile ScanBamParam
sc_atac_plot_TSSE <- function(bam, outdir, txs=NULL, seqinfo=NULL, gff=NULL, TxDb=NULL, seqnames=NULL, pairedReads=TRUE) {
	if (!dir.exists(outdir)) {dir.create(outdir)}

	if (is.null(gff) && is.null(TxDb) && is.null(txs)) {
		stop("Must provide GRanges transcripts object, gff file or a TxDb object")
	}
	if (is.null(txs)) {
		if (is.null(TxDb)) {
			# use the gff file
			TxDb <- GenomicFeatures::makeTxDbFromGFF(file=gff)
		}
		seqinfo <- GenomicFeatures::seqinfo(TxDb)
		txs <- GenomicFeatures::transcripts(TxDb)
	} else {
		# if transcripts is not null, we must also have seqinfo specified
		if (is.null(seqinfo)) {
			stop("If txs is given, must also provide seqinfo argument.")
		}
	}
	
	bamfile.labels <- gsub(".bam", "", basename(bam))

	possibleTag <- list("integer"=c("AM", "AS", "CM", "CP", "FI", "H0", "H1", "H2", 
									"HI", "IH", "MQ", "NH", "NM", "OP", "PQ", "SM",
									"TC", "UQ"), 
						"character"=c("BC", "BQ", "BZ", "CB", "CC", "CO", "CQ", "CR",
									"CS", "CT", "CY", "E2", "FS", "LB", "MC", "MD",
									"MI", "OA", "OC", "OQ", "OX", "PG", "PT", "PU",
									"Q2", "QT", "QX", "R2", "RG", "RX", "SA", "TS",
									"U2"))

	bamTop100 <- Rsamtools::scanBam(Rsamtools::BamFile(bam, yieldSize = 100), 
									param=Rsamtools::ScanBamParam(tag=unlist(possibleTag)))[[1]]$tag
	tags <- names(bamTop100)[lengths(bamTop100)>0]
	# tags
	
	which <- as(if (is.null(seqnames)) seqinfo else seqinfo[seqnames], "GRanges")

	gal <- ATACseqQC::readBamFile(bam, tag=tags, which=which, asMates=pairedReads, bigFile=TRUE)
	shiftedBam <- file.path(outdir, "shifted.bam")
	if (pairedReads) {
		gal1 <- ATACseqQC::shiftGAlignmentsList(gal, outbam=shiftedBam)
	} else {
		gal1 <- ATACseqQC::shiftGAlignments(gal, outbam=shiftedBam)
	}

	
	tsse <- ATACseqQC::TSSEscore(gal1, txs)
	
	# can we make this look a lot nicer?
	plot(1000 * (-9.5:9.5), tsse$values, type="b", xlab="distance to TSS", ylab="aggregate TSS score")
	
	# return (tsse$TSSEscore) # should we return anything?
}